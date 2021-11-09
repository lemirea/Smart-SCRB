#!/usr/bin/env python3

# script to run and manage STARsolo pipeline for Smart-SCRB DGE (3' digital gene expression) and WGB (whole gene body) methods
# 7/27/20 by Andrew Lemire
# Quantitative Genomics
# Janelia Research Campus
# HHMI
# released under BDS 3-clause license assigned to Howard Hughes Medical Institute


try:
	import os
	import sys 
	import time
	# import datetime
	import argparse
	import gzip
	from collections import defaultdict
	import warnings
	import subprocess
	import string
	import random


except ImportError as e:
    print('Required module not found: %s' % (e), file = sys.stderr)
    sys.exit(2)

print('\nAll required python modules found.')


######## Functions

def determine_method(pmap):
	## determine the assay type (spia or scrb) from the input platemap or sample.key.  this needs work
	with open(pmap, 'r') as fh:
		for rec in fh:
			rec = rec.rstrip()
			bits = rec.split('\t')
			if (len(bits) == 6) and (bits[0] in ['NB551104', 'MN00399', 'SN651', 'VH00613']):	# instrument serial numbers
				return('scrb')
			else:
				return('spia')



def get_single_record_from_fq (fq):
	## extract 4 lines at a time and return a list
	readid = fq.readline().rstrip()
	seq = fq.readline().rstrip()
	spacer = fq.readline().rstrip()
	qual = fq.readline().rstrip()
	if '+' not in spacer:   # this needs improvement: if not spacer.startswith('+'):
		raise TypeError("Problem with fastq format: \n%s\n%s\n%s\n%s" % (readid, seq, spacer, qual))
	return (readid, seq, spacer, qual)



def determine_assay_format(fname):
	## give a sample name from args.samples, determine if it is PE or SR, and DGE or WGB
	## assumes fastqs are named with Janelia Quantitative Genomics file naming convention
	assay = ''
	format = 'sr'
	read1 = gzip.open(fname + '_1.fastq.gz', 'rt')
	(rheader, rseq, rspacer, rqual) = get_single_record_from_fq(read1)
	if os.path.isfile(fname + '_2.fastq.gz'): 
		format = 'pe'
		read2 = gzip.open(fname + '_2.fastq.gz', 'rt')
		(mheader, mseq, mspacer, mqual) = get_single_record_from_fq(read2)
	if format == 'pe':
		if (len(rseq) < 30) and (len(rseq) != len(mseq)):
			assay = 'dge'
		else:
			assay = 'wgb'
	else:
		assay = 'wgb'
	read1.close()
	if os.path.isfile(fname + '_2.fastq.gz'): read2.close()
	return(assay, format)


def elapsed_time(start):
	## calculates elapsed time for a segment of code and returns tuple of (h, m, s)
	# use with print 'Elapsed time: %d:%02d:%02d' % (h,m,s) 
	end = time.perf_counter()
	elapsed = end - start
	m, s = divmod(elapsed, 60)
	h, m = divmod(m, 60)
	return (h,m,s)



#################################################################
########	Main body            ################################
#################################################################



in_the_beginning = time.perf_counter()

if not sys.argv[1:]:
	print('\n', '*' * 48)
	print(sys.argv[0], ' requires at least TWO arguments (the species index, and the platemap.tab or sample_key.txt).  Exiting.', '\n', '*' * 48)
	sys.exit()

parser = argparse.ArgumentParser(description='Configure %s script' % sys.argv[0])
# parser.add_argument('method', nargs='?', type=str, help='The assay used (spia or smrtscrb).')
# parser.add_argument('--Qindex', '-q', nargs='?', type=argparse.FileType('r'), help='The Q.index barcode file showing seq,plate,well.')
parser.add_argument('species', nargs='?', type=str, help='The species index to use. Must be one of [\'dmel\', \'mmus\', \'rnor\', \'drer\', \'hsap\', \'rtr\']')
parser.add_argument('samples', nargs='?', type=argparse.FileType('r'), help='The platemap.tab or sample.key file showing fastq base names and sample info.')
parser.add_argument('--method', '-m', nargs='?', type=str, help='Declare method used ("spia" or "scrb") or the program will attempt to determine if from the \
	samples file; the program will automatically determine DGE or WGB. This flag helps the program determine how to read the samples argument.')
# parser.add_argument('--kindex', '-k', nargs='?', type=argparse.FileType('r'), help='The path to the kallisto transcriptome index to use.')
parser.add_argument('--notrim', '-n', action='store_true', help='Use this flag to skip cutadapt trimming step (why would you want to do that?).')
parser.add_argument('--cores', '-j', nargs='?', type=int, help='The number of cores to request for the job (default = 4).')
parser.add_argument('--geneFull', '-g', action='store_true', help='Use this flag to count STARsolo in GeneFull mode (default = Gene mode)', default=False)
parser.add_argument('--echo', '-e', action='store_true', help='Turn this on to echo jobs instead of submitting them (default = submit jobs)', default=False)
# parser.add_argument('--gbam', '-g', nargs='?', type=argparse.FileType('r'), help='Use this flag to genearate a genome bam and pass the GTF to use')

args = parser.parse_args()

#################################################################
########	Parse args and prepare files    #####################
#################################################################

# validate input options and report them or die
# check if species index exists
baserefpath = '/groups/quantitativegenomics/sequencers/seq1/analysis/'	# path to genome references with STAR & bowtie1 indices, by convention
species_lut = {
	'dmel' : 'Dmelanogaster',
	'mmus' : 'Mus_musculus',
	'rnor' : 'Rattus_norvegicus',
	'drer' : 'Drerio',
	'hsap' : 'Homo_sapiens',
	'rtr' : 'RTR_rat'
}

rsem_lut = {
	'dmel' : 'dmel_flybase.r6.34.ERCCtg',
	'mmus' : 'mm10_plus.fa',
	'rnor' : 'rnor6_plus.fa',
	'drer' : 'danRer11_plus.fa',
	'hsap' : 'hg38_plus.fa',
	'rtr' : 'rnor6_plus.fa'
}
gtf_lut = {
	'dmel' : 'dmel_flybase.r6.34.ERCCtg.gtf',
	'mmus' : 'mm10_plus.gtf',
	'rnor' : 'rnor6_plus.gtf',
	'drer' : 'danRer11_plus.gtf',
	'hsap' : 'hg38_plusboth.gtf',
	'rtr' : 'rtr_plus.gtf'
}
refpath = baserefpath + 'reference_genomes/' + species_lut[args.species] + '/star_idx/'
rsempath = baserefpath + 'reference_genomes/' + species_lut[args.species] + '/bowtie1_idx/' + rsem_lut[args.species]
gtfpath = refpath + gtf_lut[args.species]
gid2gn = refpath + args.species + '.gid2gn'
internaldatapath = baserefpath + 'shared_software/'

# validate arguments
if not args.species:
	print('Species index is required! Exiting.', '\n', '*' * 48)
	sys.exit()

if not os.path.isdir(refpath):
	print('Species index %s does not exist!  Exiting.' % (refpath))
	print('*' * 48)
	sys.exit()

if args.geneFull:
	countMethod = 'GeneFull'
else:
	countMethod = 'Gene'

if args.echo:
	echo_cmd = 'echo '
else:
	echo_cmd = ''

# check method and set the platemap type. spia typically has 2 or 3 columns in sample.key, scrb (dge) always has 6 columns in platemap.tab. scrb-WGB is treated as spia
if (args.method) and (args.method not in ['spia', 'scrb']):
	print('--method argument must be one of [\'spia\', \'scrb\'].  Exiting.\n', '*' * 48)
	sys.exit()
elif args.method:
	method_type = args.method
else:
	method_type = determine_method(args.samples.name)

# find the uid for /scratch/ purposes:
myuser = subprocess.check_output('whoami', shell=True).strip().decode('utf-8')
scratchout = '/scratch/' + myuser + '/'

# check if the core request is reasonable (<= 8)
if args.cores:
	if args.cores > 8:
		print('You have requested more than 8 cores for the STARsolo pipline.  Code does not support this many cores (not needed!). Exiting.\n', '*' * 48)
		sys.exit()
else:
	args.cores = 4

# define adapter sequences for trimming Smart-SCRB (DGE and WGB) and SPIA libraries (SPIA is TruSeq, Smart-SCRB is Nextera):
revcompTruSeq_Univ_Adap = 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT'
TruSeq_Idx_Adap = 'GATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG'
pMENTS = 'CTGTCTCTTATACACATCT'
common_Tn5_end = 'AGATGTGTATAAGAGACAG'

# these files are kept in this location
qfile = internaldatapath + 'Q.index'
whitelist = internaldatapath + 'smartscrb_whitelist.txt'

qidx = defaultdict(lambda: defaultdict())		# {qidx[plate]['well'] : well, qidx[plate]['seq'] : seq}
with open(qfile, 'r') as fh:
	for rec in fh:
		rec = rec.rstrip()
		bits = rec.split('\t')
		bitstring = '_'.join(bits[1:])
		qidx[bitstring]['plate'] = bits[1]
		qidx[bitstring]['well'] = bits[2]
		qidx[bitstring]['seq'] = bits[0]

# this method of getting a list of fastqs depends on reading it in from platemap, which could be platemap.tab from smrtscrb2 or sample.key from run_bulk pipelines
pmap = defaultdict(lambda: defaultdict())
fqlist = []
bclut = defaultdict(lambda: defaultdict())

with args.samples as fh:
	if method_type == "scrb":	
		# now we know how to parse the args.platemap file
		for rec in fh:
			rec = rec.rstrip()
			bits = rec.split('\t')
			recstring = '_'.join(bits[2:5])
			bitstring = '_'.join([bits[4],bits[3]])
			bc = qidx[bitstring]['seq']
			pmap[recstring]['instrument'] = bits[0]
			pmap[recstring]['run'] = bits[1]
			fqlist.append(bits[2])
			pmap[recstring]['well'] = bits[3]
			pmap[recstring]['qplate'] = bits[4]
			pmap[recstring]['barcode'] = bc
			pmap[recstring]['sample'] = bits[5]
			bclut[bc]['bitstring'] = bitstring
			bclut[bc]['sample'] = bits[5]
	else:
		for rec in fh:
			rec = rec.rstrip()
			bits = rec.split('\t')
			pmap['basename'] = bits[0]
			if bits[1]:
				pmap['sample'] = bits[1]
			else:
				pmap['sample'] = ''
			try:	# sometimes there is a third column, but never a fourth, by convention
				if bits[2]:
					pmap['meta'] = bits[2]
				else:
					pmap['meta'] = ''
			except IndexError:
				pass
			fqlist.append(bits[0])

fqlist = sorted(list(set(fqlist)))
fq_scratch = defaultdict()	# for creating scratch directories later

# this means bclut is only as big as platemap.tab. or another way: this program only runs the fastqs declared in samples argument file
# regardless of what else is in the directory
# now that we have the method (scrb or spia) and the samples (platemap.tab or sample.key or similar) we need to determine the assay type (DGE or WGB)
# based on the fastqs.  if only 1 fastq, it must be WGB but single-ended in which case we will assum 600+/-200 bp insert size always
# known data types:
# spia wgb SR
# spia wgb PE
# scrb wgb PE (assume SR is fine too)
# scrb dge PE
# scrb dge SR becomes scrb wgb if SR (assumes upstream rename *_2.fastq.gz as *_1.fastq.gz by hand)

if os.path.isfile(fqlist[0] + '_1.fastq.gz'):
	(assay, format) = determine_assay_format(fqlist[0])
	fileSize = os.path.getsize(fqlist[0] + '_1.fastq.gz')
else:
	print('First fastq (%s) listed in %s is not recognized. Exiting.\n%s\n' % (fqlist[0] + '_1.fastq.gz', args.samples.name, fqlist), '*' * 48)
	sys.exit()

print(time.ctime())
print('Initializing %s with the following settings:\n' % (sys.argv[0]))

# for k in vars(args).keys():
# 	if (k == "method") and not args.method:
# 		print('%s\t::\t%s' % ('method', method_type))
# 	elif k in ['samples']:
# 		try:
# 			print('%s\t::\t%s' % (k, vars(args)[k].name))
# 		except: pass
# 	else: 
# 		print('%s\t::\t%s' % (k, vars(args)[k]))

print('%s\t::\t%s' % ('species', args.species))
print('%s\t::\t%s' % ('samples', args.samples.name))
print('%s\t::\t%s' % ('method', method_type.upper()))
print('%s\t::\t%s' % ('format', format.upper()))
print('%s\t::\t%s' % ('assay', assay.upper()))
print('%s\t::\t%s' % ('notrim', args.notrim))
print('%s\t::\t%s' % ('cores', args.cores))
print('%s\t::\t%s' % ('geneFull', args.geneFull))
print('%s\t::\t%s' % ('echo', args.echo))
print('%s\t::\t%s' % ('scratch', scratchout))



# trim the reads:

trim_assays = {
	'dge' : ['cutadapt', '-j', str(args.cores), '-m', '20', '--trim-n', '-n', '1', '-a', 'pMENTS='+pMENTS, '-A', 'revcompTruSeq_Univ_Adap='+revcompTruSeq_Univ_Adap],
	'wgb' : ['cutadapt', '-j', str(args.cores), '-m', '36', '--trim-n', '-n', '3', '-a', 'revcompTruSeq_Univ_Adap='+revcompTruSeq_Univ_Adap, \
	'-a', 'TruSeq_Idx_Adap='+TruSeq_Idx_Adap, '-a', 'pMENTS='+pMENTS, '-a', 'common_Tn5_end='+common_Tn5_end]
}

trimstring = trim_assays[assay]
if (assay == 'wgb') and (format == 'pe'):
	trimstring.append('-A')
	trimstring.append('revcompTruSeq_Univ_Adap='+revcompTruSeq_Univ_Adap)
	trimstring.append('-A')
	trimstring.append('TruSeq_Idx_Adap='+TruSeq_Idx_Adap)
	trimstring.append('-A')
	trimstring.append('pMENTS='+pMENTS)
	trimstring.append('-A')
	trimstring.append('common_Tn5_end='+common_Tn5_end)

# define the software parameters by method type
star_methods = {
	'spia' : 'STAR --genomeDir %s --runThreadN %s \
--sjdbGTFfile %s \
--alignSJoverhangMin 8 \
--alignSJDBoverhangMin 3 \
--outFilterMismatchNmax 999 \
--alignIntronMin 20 \
--alignIntronMax 150000 \
--readFilesCommand zcat \
--outMultimapperOrder Random \
--outSAMmultNmax 20 \
--outSAMattributes NH HI NM MD \
--outSAMtype BAM Unsorted \
--quantMode TranscriptomeSAM GeneCounts \
--quantTranscriptomeBan IndelSoftclipSingleend \
--outSAMunmapped Within ' % (refpath, args.cores, gtfpath),
	'scrb' : 'STAR --genomeDir %s --runThreadN %s \
--alignSJoverhangMin 8 \
--alignSJDBoverhangMin 8 \
--alignIntronMin 20 \
--alignIntronMax 150000 \
--readFilesCommand zcat \
--outMultimapperOrder Random \
--outSAMmultNmax 20 \
--outSAMattributes NH HI NM MD CR CY UR UY \
--outSAMtype BAM Unsorted \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--soloType CB_UMI_Simple \
--soloCBwhitelist %s \
--soloCBstart 2 \
--soloCBlen 8 \
--soloUMIstart 10 \
--soloUMIlen 10 \
--soloBarcodeReadLength 0 \
--soloCBmatchWLtype 1MM_multi_pseudocounts \
--soloStrand Forward \
--soloFeatures %s \
--soloUMIdedup 1MM_All \
--soloUMIfiltering MultiGeneUMI \
--soloCellFilter None' % (refpath, args.cores, whitelist, countMethod)
}
r_assays = {
	'dge' : '~/bin/STARsolo_DGE_downstream.R',
	'wgb' : '~/bin/STAR_WGB_downstream.R'
}

rsem_command = ['rsem-calculate-expression', '-p', str(args.cores)]
# if PE, let RSEM discover fragment length distribution, otherwise declare it
# if format == 'sr':
# 	rsem_command.append('--fragment-length-mean')
# 	rsem_command.append('600')
# 	rsem_command.append('--fragment-length-sd')
# 	rsem_command.append('200')

rsem_command.append('--fragment-length-mean')
rsem_command.append('600')
rsem_command.append('--fragment-length-sd')
rsem_command.append('200')

rsem_command.append('--no-bam-output')
rsem_command.append('--alignments')
# create a list of all job commands to be executed
bsub_outs = open('bsub_commands.txt', 'w')
if not os.path.isdir('./logs'):
	print('Creating logs directory.')
	os.system('mkdir logs')


# loop through all fastqs in args.samples and create bash scripts for each
for i in range(0,len(fqlist)):
	fq_scratch[fqlist[i]] = ''.join(random.choices(string.ascii_uppercase + string.ascii_lowercase + string.digits, k=10))
	read1 = fqlist[i] + '_1.fastq.gz'
	read2 = fqlist[i] + '_2.fastq.gz'
	star_cmd = [star_methods[method_type]]
	r_cmd = ['Rscript', r_assays[assay]]
	rsem_string = [' '.join(rsem_command)]
	if os.path.isfile(read1):
		print('Generating bash script for %s...' % fqlist[i], end = '')
		tmpout = scratchout + fqlist[i] + '_' + fq_scratch[fqlist[i]] + '/'
		if method_type == 'dge':
			trimout_for = tmpout + 'catrim.m20.' + read1
			trimout_rev = tmpout + 'catrim.m20.' + read2
		else:
			trimout_for = tmpout + 'catrim.m36.' + read1
			trimout_rev = tmpout + 'catrim.m36.' + read2

		# create bashfile and header, basic tests and housekeeping
		bashfile = open(fqlist[i] + '.sh', 'w')
		bashfile.write('#!/bin/bash\n')
		bashfile.write('# bash script to analyze %s from %s method with %s assay in %s mode with %s as the reference, created: %s\n' % (fqlist[i], method_type.upper(), assay.upper(), format.upper(), args.species, time.ctime()))
		bashfile.write('(test -d %s && echo "Results directory %s already exists.") || (mkdir %s && echo "Created %s directory for results." && echo)\n\n' % ('./star_analysis', './star_analysis', './star_analysis', './star_analysis'))
		bashfile.write('echo \'# Starting contents of scratch directory:\'\n')
		bashfile.write('ls -lah %s\n\n' % scratchout)
		bashfile.write('mkdir %s\n' % tmpout)
		# create cutadapt command
		if not args.notrim:
			bashfile.write('\necho \'# Begin cutadapt log:\'\n\n')
			bashfile.write(' '.join(trimstring))
			bashfile.write(' -o %s' % trimout_for)
			if format == 'sr':
				bashfile.write(' %s\n\n' % read1)
			else:
				bashfile.write(' -p %s' % trimout_rev)
				bashfile.write(' %s' % read1)
				bashfile.write(' %s\n\n' % read2)
		else:
			bashfile.write('\necho \'# CUTADAPT STEP WAS SKIPPED BY USER %s\'\n\n' % myuser)
			trimout_for = read1
			trimout_rev = read2
		# this will record users who skip adapter trimming (don't skip it...it's there for testing purposes or re-analysis from trimmed reads)
		# create STAR command
		bashfile.write('\necho \'# Begin STAR log:\'\n\n')		
		star_cmd.append('--readFilesIn')
		if method_type == 'scrb':	# STARsolo requires cDNA read as 1st input, BC+UMI read as 2nd input which is the reverse of the actual read order (read 1 is BC+UMI, 2 is cDNA)
			star_cmd.append(trimout_rev)
			star_cmd.append(trimout_for)
		else:
			star_cmd.append(trimout_for)
			if format == 'pe':
				star_cmd.append(trimout_rev)
		star_cmd.append('--outFileNamePrefix')
		star_cmd.append(tmpout + fqlist[i] + '.')
		bashfile.write(' '.join(star_cmd) + '\n\n')
		bashfile.write('echo \'# Scratch directory contents after STAR:\'\n')
		bashfile.write('ls -lah %s\n\n' % tmpout)

		if method_type == 'spia':
			if format == 'pe':
				rsem_string.append('--paired-end')
			rsem_string.append(tmpout + fqlist[i] + '.Aligned.toTranscriptome.out.bam')
			rsem_string.append(rsempath)
			rsem_string.append(tmpout + fqlist[i] + '.star.rsem')
			bashfile.write('\necho \'# Begin RSEM log:\'\n\n')
			bashfile.write(' '.join(rsem_string) + '\n\n')
			bashfile.write('echo \'# Scratch directory contents after RSEM:\'\n')
			bashfile.write('ls -lah %s\n\n' % tmpout)
					
# add downstream steps here

		# create samtools command
		bashfile.write('\necho \'# Begin samtools log:\'\n\n')
		samtools_cmd = ['samtools', 'sort', '-T', tmpout + fqlist[i], '-m', '6G', '-o', tmpout + fqlist[i] + '.sorted.bam', tmpout + fqlist[i] + '.Aligned.out.bam']
		bashfile.write(' '.join(samtools_cmd) + '\n\n')
		samtools_cmd = ['samtools', 'index', tmpout + fqlist[i] + '.sorted.bam']
		bashfile.write(' '.join(samtools_cmd) + '\n\n')
		samtools_cmd = ['samtools', 'flagstat', tmpout + fqlist[i] + '.sorted.bam', '>', tmpout + fqlist[i] + '.sorted.bam.flagstat']
		bashfile.write(' '.join(samtools_cmd) + '\n\n')
		# create cleanup command: delete trimmed fastqs
		cleanup_string = []
		if not args.notrim:
			cleanup_string.append('rm')
			cleanup_string.append('-v')
			cleanup_string.append(trimout_for)
			if format == 'pe':
				cleanup_string.append(trimout_rev + '\n')
		cleanup_string.append('rm')
		cleanup_string.append('-v')
		cleanup_string.append(tmpout + '*.Aligned*.out.bam;\n')
		cleanup_string.append('mv')
		cleanup_string.append('-v')		
		# here we check if the results directory (./star_analysis/fqlist[i]) already exists and clobber if so.
		if os.path.isfile('./star_analysis/%s' % (fqlist[i] + '.sorted.bam')):
			bashfile.write('echo \'Clobbering previous results!\'\n\n')
			print('Clobbering previous results!')
		cleanup_string.append(tmpout + '*')
		cleanup_string.append('./star_analysis/;\n')
		bashfile.write('\necho \'# Begin cleanup log:\'\n\n')
		bashfile.write(' '.join(cleanup_string) + '\n\n')
		# housekeeping
		bashfile.write('echo \'# Temp directory contents after cleanup:\'\n')
		bashfile.write('ls -lah %s\n\n' % tmpout)
		bashfile.write('rmdir %s\n\n' % tmpout)
		bashfile.write('echo \'# Scratch directory contents after cleanup:\'\n')
		bashfile.write('ls -lah %s\n\n' % scratchout)
		# create Rscript command
		if method_type == 'spia':
			bashfile.write('echo \'# Copying %s gene_id to gene_name LUT:\'\n' % (args.species))
			bashfile.write('scp -v %s ./' % (gid2gn))
			# do spia workflow
			# there is no downstream bustools for spia, but here we build the downstream Rscript command.
			# r_cmd = ['Rscript', '~/bin/kallisto_WGB_downstream.R', args.samples.name, t2g, tmpout + 'abundance.tsv', fqlist[i]]
			# r_cmd.append(args.samples.name)
			# r_cmd.append(t2g)
			# r_cmd.append(tmpout + 'abundance.tsv')
			# r_cmd.append(fqlist[i])

		else:   # scrb, so do STARsolo workflow
			# do scrb workflow
			# r_cmd = ['Rscript', '~/bin/STARsolo_DGE_downstream.R', qfile, args.samples.name, rscript_gid2gn, tmpout + fqlist[i] + '.genes.txt', tmpout + fqlist[i] + '.barcodes.txt', tmpout + fqlist[i] + '.mtx', fqlist[i]]
			bashfile.write('chmod -R 755 ./star_analysis/' + fqlist[i] + '.Solo.out/\n')
			r_cmd.append(qfile)
			r_cmd.append(args.samples.name)
			r_cmd.append('./star_analysis/' + fqlist[i] + '.Solo.out/' + countMethod + '/raw/' + 'features.tsv')
			r_cmd.append('./star_analysis/' + fqlist[i] + '.Solo.out/' + countMethod + '/raw/' + 'barcodes.tsv')
			r_cmd.append('./star_analysis/' + fqlist[i] + '.Solo.out/' + countMethod + '/raw/' + 'matrix.mtx')
			r_cmd.append(fqlist[i])
			r_cmd.append('./star_analysis/')
			bashfile.write('echo \'# Begin Rscript log:\'\n\n')
			bashfile.write(' '.join(r_cmd) + '\n\n')
			bashfile.write('echo \'# Results directory contents after Rscript:\'\n')
			bashfile.write('ls -lah %s\n\n' % './star_analysis/')

		print('done!')
		bashfile.close()
		os.system('chmod a+x %s' % bashfile.name)
		os.system(echo_cmd + 'bsub -n %s -J %s -o ./logs/%s "./%s"' % (args.cores, fqlist[i], fqlist[i] + '.\%J.out', bashfile.name))
		bsub_outs.write('bsub -n %s -J %s -o ./logs/%s "./%s" \n\n' % (args.cores, fqlist[i], fqlist[i] + '.\%J.out', bashfile.name))
		bashfile.close()
	else:
		print('Skipping sample file not found: %s' % (read1))
	

print ('\nFinished creating bash scripts and submitting jobs...')
bsub_outs.close()


print('Done!  Elapsed time for entire program (h:m:s): %d:%02d:%02d' % (elapsed_time(in_the_beginning)))




