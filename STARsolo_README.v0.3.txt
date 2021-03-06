###  README for output files from STARsolo pipeline
###  v.03

A brief outline of the STARsolo pipeline for Smart-SCRB DGE (3' digital gene expression):
1. declare species and platemap.tab (sample info) to python wrapper (run_STARsolo_pipeline.py)
2. python creates bash scripts to manage all steps for each fastq
	a. trim sequencing adapters with cutadapt
	b. STARsolo + whitelisted barcodes + species reference = counts matrix for all barcodes in whitelist
	c. sort and index genome bam for genome browsing
	d. run Rscript STARsolo_DGE_downstream.R to extract barcodes declared in platemap.tab and summarize results
3. final Rscript combine_STARsolo_DGE_results.R to merge results once all jobs are completed

Sample information is declared in platemap.tab (by QGSR staff) using a standard convention in QGSR. There are no headers in platemap.tab.
Tnere are six columns in platemap.tab:

1	sequencing instrument serial number
2	sequencing instrument run number
3	base fastq name (created by QGSR convention) [the * in *_1.fastq.gz]
4	Q.plate well position used
5	Q.plate used
6	sample identifier (created by QGSR staff with input from requestor)

There are four 96 well "Q.plates" in QGSR containing the set of 384 barcodes designed for our Smart-SCRB assay.  A copy of this file
is appended to the bottom of this README.  The first column of Q.index is used as the barcode whitelist passed to STARsolo.

The platemap.tab shows which well and Q.plate was used for each sample and the name of the fastq containing the sample reads.  
For Smart-SCRB DGE each fastq generally contains between 8 (low-cell) and 96 (single-cell) barcoded samples.  This file is located with the fastqs.

The STARsolo pipeline is managed by run_STARsolo_pipeline.py (QGSR).  This wrapper takes the species (dmel, mmus, hsap, rnor, drer, etc.) and the
platemap.tab file as inputs and creates a bash script to process each fastq separately.  It also creates a list of the commands
used to launch each bash script job.

Example files created by python wrapper, where fastq1 and fastq2 are two different libraries (fastqs):

 - bsub_commands.txt
 - <fastq1>.sh
 - <fastq2>.sh

There will be one bash script (and one entry in bsub_commands.txt) per fastq declared in platemap.tab.  These can also be used
as logs to show the software and parameters used.  Each bash script calls an R script (STARsolo_DGE_downstream.R) to summarize the 
STARsolo output, described below.

run_STARsolo_pipeline.py creates 'star_analysis' and 'logs' directories in the same location as the fastqs.  star_analysis contains the 
output from STARsolo (bam, counts matrix, logs, Rscript outputs, etc) and logs contains the wrapper logs.

These three files are created by the STARsolo_DGE_downstream.R script (QGSR) and are described along with the remaining STARsolo files below:
 - <fastq>.contam.txt
 - <fastq>.counts.txt
 - <fastq>.pstats.txt

The remaining files from STARsolo look like this, where <fastq> is the base fastq name:

<fastq>.contam.txt		barcodes from other Q.plates and wells that are not declared in the platemap.tab (typically lab contamination)
<fastq>.counts.txt		per-fastq counts for all samples in platemap.tab
<fastq>.Log.final.out		STAR final log
<fastq>.Log.out			STAR input parameters
<fastq>.Log.progress.out	STAR progress log
<fastq>.pstats.txt		per-fastq extension of platemap.tab; a subset of merged_sample_info.txt for this fastq only
<fastq>.ReadsPerGene.out.tab	STAR output file
<fastq>.SJ.out.tab		STAR output file
<fastq>.sorted.bam		gapped-alignment to full genome, ready for viewing in your favorite genome browser (IGV, UCSC, etc)
<fastq>.sorted.bam.bai		bam index, required for genome browsing
<fastq>.sorted.bam.flagstat	samtools flagstat output for this bam

<fastq>.Solo.out:		STARsolo output directory
Barcodes.stats			STARsolo output file
Gene (or GeneFull)		STARsolo output directory
    Features.stats		STARsolo output file; info is copied to merged_fastq_info.txt
    Summary.csv			STARsolo output file; some info is copied to merged_fastq_info
    raw				STARsolo output directory
        barcodes.tsv		STARsolo output file; these are column headers for unpacked matrix.mtx
        features.tsv		STARsolo output file; these are row names for unpacked matrix.mtx
        matrix.mtx		STARsolo output file; this is the counts matrix in MatrixMarket format (https://math.nist.gov/MatrixMarket/formats.html#MMformat)

The matrix.mtx generally contains all 384 barcodes from the whitelist (unless the QGSR parameters have changed).  The STARsolo_DGE_downstream.R script. 
(QGSR) uses the platemap.tab file to extract the declared barcodes from the /raw output from STARsolo.  The remainder are reported in <fastq>.contam.txt
if there are at least 100 counts.  Contamination could arise from aerosols when opening/closing barcode oligo plates, and/or from re-using Nextera i7 barcode
combinations with Q.plate combinations.  QGSR has 4 Q.plates of 96 barcodes each (384 total) and there are 24 Illumina Nextera i7 barcodes.  Over time, the 
same Q.plate-i7 combinations arise and may show up as contamination.  QGSR uses the metrics in these reports to discard lab water, barcode dilution plates,
and other lab reagents as necessary to minimize barcode background.

Once all jobs are complete, if desired, another R script (combine_STARsolo_DGE_results.R) is run to merge all results for the platemap.tab input.  These 
three files are created by the combine_STARsolo_DGE_results.R script (QGSR):

 - merged_counts.txt is the raw counts imported from STARsolo outputs (features.tsv, barcodes.tsv, matrix.mtx), 
	for all samples declared in platemap.tab.
 - merged_sample_info.txt is an extension of platemap.tab with post-processing stats, and is the sample key for merged_counts.txt.
 - merged_fastq_info.txt contains summary stats for all fastqs in the platemap.tab file.

The column headers in merged_counts.txt are the same as "qp_well_fq" in merged_sample_info.txt.
These are the columns in merged_counts.txt:

1	gene_id			gene_id from GTF file used in alignment and counting
2	gene_name		gene_name from GTF file used in alignment and counting
...	<qp_well_fq>	the remaining columns are counts data labeled by qp_well_fq column in merged_sample_info.txt

These are the columns in merged_sample_info.txt:

1	qp_well_fq	Q.plate + well + fastq name for this sample, used to link merged_sample_info.txt with merged_counts.txt
2	fastq		base fastq name for this sample
3	well		Q.plate well position for this sample
4	q_plate		Q.plate used for this sample
5	barcode		barcode for this sample
6	sample		sample identifier
7	inst		sequencing instrument serial number
8	run_num		sequencing instrument run number
9	gene_det	number of genes detected (count > 0) for this sample
10	sum_count	sum of all counts for this sample
11	pct_ercc	percent of counts that are ERCC controls for this sample
12	num_ercc	number of ERCC transcripts detected out of 92 for this sample
13	yield_ratio	ratio of sum_count to total_reads for this fastq (from merged_fastq_info.txt) = the useable fraction of reads

These are the columns in merged_fastq_info.txt

1	fastq			base fastq name (<fastq>_1.fastq.gz and <fastq>_2.fastq.gz)
2	total_reads		total starting reads in this fastq before any processing
3	trimmed_reads		total reads after cutadapt trimming
4	q30_bc_umi		ratio of BC+UMI reads to trimmed_reads if read quality is Q30+
5	q30_cDNA		ratio of cDNA reads to trimmed_reads if read quality is Q30+
6	uniq_aln_read_num	uniquely aligned reads
7	nUnmapped		extracted from STARsolo Features.stats
8	nNoFeature		extracted from STARsolo Features.stats
9	nAmbigFeature		extracted from STARsolo Features.stats
10	nAmbigFeatureMultimap	extracted from STARsolo Features.stats
11	nTooMany		extracted from STARsolo Features.stats
12	nNoExactMatch		extracted from STARsolo Features.stats
13	nExactMatch		extracted from STARsolo Features.stats
14	nMatch			extracted from STARsolo Features.stats
15	nCellBarcodes		extracted from STARsolo Features.stats = number of whitelisted barcodes reported
16	nUMIs			extracted from STARsolo Features.stats
17	map2txome		extracted from STARsolo Summary.csv = equal to ratio of nMatch to trimmed_reads
18	saturation		extracted from STARsolo Summary.csv
19	sum_all_count		this is the per-fastq sum of (sum_count) in merged_sample_info.txt
20	yield_ratio		this is the ratio of sum_all_count to total_reads for this fastq, the fraction of reads containing information

For Janelians: Please contact QGSR staff if you have any questions about these files, the R scripts, or the python wrapper used to run the STARsolo pipeline.

A copy of this README is provided with every Smart-SCRB analysis run performed by QGSR.

### INFO

cutadapt version:
# packages in environment at /groups/quantitativegenomics/home/qguser/bin/anaconda3:
#
# Name                    Version                   Build  Channel
cutadapt                  2.10             py37hf01694f_1    bioconda


STAR version:
## versions
versionGenome           2.7.4a
STAR --version
2.7.5c



python3 libraries used:
os, sys, time, subprocess, argparse, gzip, collections, string, random

R packages used:
pacman, ggplot2, dplyr, psych, Matrix

The R sessionInfo() is logged in logs/<fastq>.%J.out where %J is LSF job ID.

Three QGSR wrapper scripts:
run_STARsolo_pipeline.py
STARsolo_DGE_downstream.R
combine_STARsolo_DGE_results.R


The helpfile for run_STARsolo_pipeline.py:

Configure
/groups/quantitativegenomics/home/lemirea/bin/run_STARsolo_pipeline.py script

positional arguments:
  species               The species index to use. Must be one of ['dmel',
                        'mmus', 'rnor', 'drer', 'hsap']
  samples               The platemap.tab or sample.key file showing fastq base
                        names and sample info.

optional arguments:
  -h, --help            show this help message and exit
  --method [METHOD], -m [METHOD]
                        Declare method used ("spia" or "scrb") or the program
                        will attempt to determine if from the samples file;
                        the program will automatically determine DGE or WGB.
                        This flag helps the program determine how to read the
                        samples argument.
  --notrim, -n          Use this flag to skip cutadapt trimming step (why
                        would you want to do that?).
  --cores [CORES], -j [CORES]
                        The number of cores to request for the job (default =
                        4).
  --geneFull, -g        Use this flag to count STARsolo in GeneFull mode
                        (default = Gene mode)
  --echo, -e            Turn this on to echo jobs instead of submitting them
                        (default = submit jobs)


### Q.index
AACGTCTT	Q5	A01
GGCGTCAT	Q5	A02
AATGGCGT	Q5	A03
ATTGCAGT	Q5	A04
CACGGAGT	Q5	A05
CACGATGC	Q5	A06
ATTGGCTA	Q5	A07
CATGCATT	Q5	A08
CTCGCTTT	Q5	A09
GGGATCTT	Q5	A10
CGGATATA	Q5	A11
CACGCAAG	Q5	A12
CACCTATT	Q5	B01
CTTCTCTT	Q5	B02
CGAGTAGA	Q5	B03
TACGCTAC	Q5	B04
TTAGTCGA	Q5	B05
CGGAGGAA	Q5	B06
GGCATGTA	Q5	B07
AATGTGGA	Q5	B08
TTATCCGT	Q5	B09
CGTCCGAA	Q5	B10
ATAGCTCT	Q5	B11
ATCGGTAC	Q5	B12
AGTGCTTT	Q5	C01
TGCACGTC	Q5	C02
GATGCTTC	Q5	C03
TCCACTAA	Q5	C04
TGAGAACA	Q5	C05
TAGCTAGT	Q5	C06
CGGACTTC	Q5	C07
CAGAAGGC	Q5	C08
ATCGTCGC	Q5	C09
CTTAGCTG	Q5	C10
CTCACACT	Q5	C11
ACTGCTAC	Q5	C12
CTGAAGTA	Q5	D01
GTCCTAGC	Q5	D02
TGTCACAA	Q5	D03
TTCTGGGT	Q5	D04
GCGGTATA	Q5	D05
TGAATGGA	Q5	D06
ACTTCGAA	Q5	D07
ACCGCGTT	Q5	D08
TGGTGGCT	Q5	D09
CGAAGTCC	Q5	D10
TGTTCGCC	Q5	D11
TTACCGGA	Q5	D12
TCGTGGAA	Q5	E01
TTAACTGG	Q5	E02
CCGTGAGA	Q5	E03
TGGAGAGA	Q5	E04
AACGAAAC	Q5	E05
TAAGACGT	Q5	E06
GACAGATA	Q5	E07
ACGATTTC	Q5	E08
AACCATTC	Q5	E09
CAACAGGA	Q5	E10
ACCGATAG	Q5	E11
TCGCAGGT	Q5	E12
GACACCAT	Q5	F01
TGCGGAAT	Q5	F02
CTCTCTGA	Q5	F03
CGTCGAAG	Q5	F04
ATTTCCGC	Q5	F05
CAAGGTCA	Q5	F06
TCCATGCA	Q5	F07
CTTTTGTG	Q5	F08
AACTGCTC	Q5	F09
TTCGAGCT	Q5	F10
CTCCACGA	Q5	F11
TTTTGCGG	Q5	F12
TTGATGCT	Q5	G01
GTCGGTTG	Q5	G02
GTTTCCTT	Q5	G03
GTCGATGT	Q5	G04
TCCAGAAG	Q5	G05
TGAGTCAG	Q5	G06
GCAACTCA	Q5	G07
AACCGCCT	Q5	G08
TTACGGTC	Q5	G09
ACGATAGA	Q5	G10
CAAGTCGG	Q5	G11
GAATCGAA	Q5	G12
ATCTAGCC	Q5	H01
GGTAGAGC	Q5	H02
TTTGCCCA	Q5	H03
TCATGTGT	Q5	H04
TTCCCTAT	Q5	H05
GTGGAACT	Q5	H06
CTCTAGAA	Q5	H07
ACCTTTCC	Q5	H08
CCTAACCG	Q5	H09
CCGACAAG	Q5	H10
GCAATGAA	Q5	H11
GCAATTGC	Q5	H12
ATTTGTCC	Q6	A01
AGTACGGT	Q6	A02
GGACGATC	Q6	A03
CTAGCCAT	Q6	A04
TGCTCTCT	Q6	A05
GTGAATGA	Q6	A06
GTCACATG	Q6	A07
CGGTTCAT	Q6	A08
TGGATTGG	Q6	A09
GCGCAATT	Q6	A10
GGACGTCT	Q6	A11
GGTTCATG	Q6	A12
CCCAATTT	Q6	B01
CGAAGATG	Q6	B02
AGGGAGCA	Q6	B03
AGAAGCGG	Q6	B04
TGGCGTTT	Q6	B05
GCTGACAG	Q6	B06
ACAACCTC	Q6	B07
TGTTAGGT	Q6	B08
GTGGTAGG	Q6	B09
CTCTGTAT	Q6	B10
TCTGGCAC	Q6	B11
GCAGAACG	Q6	B12
GATCTGTT	Q6	C01
ACGACGCA	Q6	C02
AAAGGCAC	Q6	C03
CGTATAAC	Q6	C04
TACGACCG	Q6	C05
GTGTAGGG	Q6	C06
GAGTGTGG	Q6	C07
AAATTCGC	Q6	C08
GTTTCTCA	Q6	C09
AACTCCCG	Q6	C10
TAGTCCGG	Q6	C11
GAACTCAG	Q6	C12
AGACAAGA	Q6	D01
GTGCCAAC	Q6	D02
GTCATCTC	Q6	D03
GATTTGGC	Q6	D04
GCACATTA	Q6	D05
GCGTCCTA	Q6	D06
GTGAGCAG	Q6	D07
GAGACCCA	Q6	D08
AGCTATGC	Q6	D09
TTTGAGAC	Q6	D10
GGCAATCC	Q6	D11
CGGCTACT	Q6	D12
TCAGGTCG	Q6	E01
GTGTTTCT	Q6	E02
AGCTTGGG	Q6	E03
AGTAAGAG	Q6	E04
CATAGCAA	Q6	E05
TATCGGAC	Q6	E06
AAGCCCAT	Q6	E07
GTCGACTA	Q6	E08
CCATTGAT	Q6	E09
ATCACAGA	Q6	E10
TGAGGCTA	Q6	E11
AGCCTCAA	Q6	E12
GTAGTTCG	Q6	F01
CCAAAAGG	Q6	F02
AGGCGAAT	Q6	F03
CGGTCACA	Q6	F04
CACAGGGA	Q6	F05
GAAGCGTG	Q6	F06
TGACTCCC	Q6	F07
AAACAGAC	Q6	F08
GGCAACGA	Q6	F09
TGGCCCTA	Q6	F10
ACCCTATC	Q6	F11
TCGCACAG	Q6	F12
GACTCACA	Q6	G01
ACATAGCG	Q6	G02
TACCAATG	Q6	G03
AGGGTCTA	Q6	G04
GAAGTGAC	Q6	G05
TATTCGAG	Q6	G06
GTAGGGTA	Q6	G07
GGAGAGGT	Q6	G08
ACACCTCC	Q6	G09
GGAAGCAA	Q6	G10
GACCTTCC	Q6	G11
GTACGCCA	Q6	G12
CGCTTGAC	Q6	H01
TCTATGGG	Q6	H02
TGGTTTTC	Q6	H03
ATACTGCA	Q6	H04
AGGTCATC	Q6	H05
GATGATCG	Q6	H06
CAGTATTG	Q6	H07
ATGTGTAG	Q6	H08
TTTCCTCG	Q6	H09
TCAACAAC	Q6	H10
AGACCACT	Q6	H11
ACGAGCTG	Q6	H12
TCCAACAC	Q7	A01
ATAACGCG	Q7	A02
CTTGAAAG	Q7	A03
ATCCAGAG	Q7	A04
GGTAATGG	Q7	A05
ACAGCATA	Q7	A06
GTCCTTAA	Q7	A07
CCCCTACA	Q7	A08
AGCATACG	Q7	A09
ACTGACTT	Q7	A10
CATTGGCG	Q7	A11
TTGCTGAA	Q7	A12
CATGTACA	Q7	B01
GATTGTAC	Q7	B02
GGCAAAAT	Q7	B03
TACCTGCG	Q7	B04
ATTGACCC	Q7	B05
ATACGAGT	Q7	B06
ATACCATG	Q7	B07
GCTGTTTT	Q7	B08
CTGCTGTC	Q7	B09
CCTGATCT	Q7	B10
GCCAACCT	Q7	B11
CTGTTAAC	Q7	B12
CACCGTAG	Q7	C01
TTTCCATC	Q7	C02
ATCCGACA	Q7	C03
AGTCTACA	Q7	C04
GAGCCTAG	Q7	C05
GAAGACAA	Q7	C06
CGTTACAG	Q7	C07
GTGAGGTC	Q7	C08
GCACTCTC	Q7	C09
TCTCTTCC	Q7	C10
TAGGGCAA	Q7	C11
CTTCCAGG	Q7	C12
CAACCTAC	Q7	D01
GCCATTAG	Q7	D02
TTCAGCTT	Q7	D03
CCTCGTTT	Q7	D04
CTTATGGA	Q7	D05
CCTGTGTC	Q7	D06
GCAGCGAT	Q7	D07
GCTTGACG	Q7	D08
CCGTTTTT	Q7	D09
TTAGCACG	Q7	D10
ACTCGGTC	Q7	D11
ACAGGGCA	Q7	D12
GGTTTGCA	Q7	E01
CCCTAAAC	Q7	E02
CTCCGTTA	Q7	E03
GCTATGCC	Q7	E04
GCATAGGC	Q7	E05
AGATGGGA	Q7	E06
TGAGTGCT	Q7	E07
CCCCCAAT	Q7	E08
AACCGGTA	Q7	E09
TACTGTCC	Q7	E10
CAAACACG	Q7	E11
CGATATTC	Q7	E12
ACCCACAT	Q7	F01
ATGGTGCG	Q7	F02
GGATCCGA	Q7	F03
AAGGACTG	Q7	F04
CGAATCTC	Q7	F05
TAGGATGA	Q7	F06
GTCTCGAG	Q7	F07
TGTTGTCG	Q7	F08
TCGCTCCT	Q7	F09
AGAGCCTG	Q7	F10
CTCATGAT	Q7	F11
CAGATGCG	Q7	F12
CACTCCAC	Q7	G01
AAAAGCCA	Q7	G02
TCACGAAA	Q7	G03
GTGCCGTT	Q7	G04
GAATCAGC	Q7	G05
ATCAGTGG	Q7	G06
AGACTTCG	Q7	G07
AGTGTCAC	Q7	G08
ATCTCCTA	Q7	G09
CGTTCCTC	Q7	G10
GCCTAGAT	Q7	G11
TCGGGACA	Q7	G12
CGTGAATC	Q7	H01
TATGGGTT	Q7	H02
TAATACCC	Q7	H03
TCACCTTT	Q7	H04
CATACTCT	Q7	H05
CTGGTTGC	Q7	H06
TAGCGACG	Q7	H07
GCTACCAC	Q7	H08
TGCCGTAC	Q7	H09
AATCACGG	Q7	H10
TCTAACGT	Q7	H11
ACTAAACC	Q7	H12
AGGATGCC	Q8	A01
CGTTAAGA	Q8	A02
TGTACCTG	Q8	A03
TCAGAGGG	Q8	A04
CATCTGAG	Q8	A05
CAGGTCCT	Q8	A06
CGAGAAAT	Q8	A07
CTTACCGT	Q8	A08
TAATCGTC	Q8	A09
ACGCAAAA	Q8	A10
TGGAACAT	Q8	A11
CAAAAGCT	Q8	A12
TGACCGAG	Q8	B01
TATGACTC	Q8	B02
GCTCGAAT	Q8	B03
GATCGTCA	Q8	B04
GTGGTCCA	Q8	B05
CCGCATTC	Q8	B06
CAGGCGAT	Q8	B07
CGAACTAA	Q8	B08
GAAAACGC	Q8	B09
TGCGTCCA	Q8	B10
AATCCTGA	Q8	B11
CATGTTTG	Q8	B12
CAGCATGT	Q8	C01
ACTGTATG	Q8	C02
TGCTGCAA	Q8	C03
ACCTCCGT	Q8	C04
AGACACTT	Q8	C05
GTTCGTGG	Q8	C06
CAGGGATA	Q8	C07
TTTGGATG	Q8	C08
CATTCTGG	Q8	C09
GCATTCGG	Q8	C10
GCTAGGAG	Q8	C11
TCACTAGG	Q8	C12
TACAGAGC	Q8	D01
TCGGCTAG	Q8	D02
GCGATGTG	Q8	D03
GGCCATTT	Q8	D04
ACATGGTT	Q8	D05
ACACTGAG	Q8	D06
CTGGCATC	Q8	D07
TGCGATTC	Q8	D08
CTACACCT	Q8	D09
TGATAAGC	Q8	D10
CCTCAGTA	Q8	D11
GCAAGACC	Q8	D12
GTTTTCGA	Q8	E01
CGTTGCCA	Q8	E02
ACCGGTGT	Q8	E03
GCATATAG	Q8	E04
ATGTACGA	Q8	E05
TTCTTTGC	Q8	E06
CCAGCACT	Q8	E07
GTATGGAT	Q8	E08
AGTAACTC	Q8	E09
TGTGTCGT	Q8	E10
TAGCTTCA	Q8	E11
CATGAGGG	Q8	E12
ACTACTGG	Q8	F01
CAACTCCA	Q8	F02
TAGAAGTG	Q8	F03
ATTGATGG	Q8	F04
CTGACCAC	Q8	F05
TGTGCAGC	Q8	F06
CACCTTGA	Q8	F07
ATCGAATG	Q8	F08
CAGTACAA	Q8	F09
TACCCAAA	Q8	F10
TATCAGCA	Q8	F11
ATGGGCAT	Q8	F12
GAGTAGCA	Q8	G01
TGAACAGT	Q8	G02
GGTGCGTA	Q8	G03
GGGAAACG	Q8	G04
CTCCGAAC	Q8	G05
GTCAGGAA	Q8	G06
ATGGGAGA	Q8	G07
GCGGATAT	Q8	G08
GCATTAAC	Q8	G09
TGGAAATC	Q8	G10
GAAGGCTT	Q8	G11
GTTGTCTG	Q8	G12
GAACTACT	Q8	H01
GGTGGACA	Q8	H02
CCGGTGAA	Q8	H03
GTAACCCT	Q8	H04
ATCCTTCT	Q8	H05
GGGTTAAA	Q8	H06
AGGGTTAG	Q8	H07
TCTAGGCT	Q8	H08
CTCAAATC	Q8	H09
CTGCAAAT	Q8	H10
ACGGATTA	Q8	H11
AATGGGAG	Q8	H12
