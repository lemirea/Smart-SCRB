###  README for output files from STAR-RSEM pipeline for NuGEN SPIA or Smart-SCRB bulk RNA-seq (whole gene body, WGB)
###  v.02

A brief outline of the STAR-RSEM pipeline:
1. declare species and sample_key.txt (sample info) to python wrapper (run_STARsolo_pipeline.py)
2. python creates bash scripts to manage all steps for each fastq
	a. trim sequencing adapters with cutadapt
	b. STAR + species reference = transcriptome bam for RSEM + genome BAM
	c. sort and index genome bam for genome browsing
	d. RSEM + transcriptome bam = gene and isoform counts matrices
3. final Rscript combine_STAR-RSEM_WGB_results.R to merge results once all jobs are completed

Sample information is declared in sample_key.txt (by QGSR staff) using a standard convention in QGSR. There are no headers in sample_key.txt.
There are two, sometimes three, columns in sample_key.txt:

1	base fastq name (created by QGSR convention)
2	sample identifier (created by QGSR staff with input from requestor)
3	optional: additional identifying info

The sample_key.txt shows the sample name and the name of the fastq containing the sample reads.  
For Smart-SCRB WGB and NuGEN SPIA each fastq is a single sample.  The sample_key.txt file is located with the fastqs.

The STAR-RSEM pipeline is managed by run_STARsolo_pipeline.py (QGSR).  This wrapper takes the species (dmel, mmus, hsap, rnor, drer, etc.) and the
sample_key.txt file as inputs and creates a bash script to process each fastq separately.  It also creates a list of the commands
used to launch each bash script job.

Example files created by python wrapper, where fastq1 and fastq2 are two different libraries (sample fastqs):

 - bsub_commands.txt
 - <fastq1>.sh
 - <fastq2>.sh

There will be one bash script (and one entry in bsub_commands.txt) per fastq declared in sample_key.txt.  These can also be used
as logs to show the software and parameters used.  

run_STARsolo_pipeline.py creates 'star_analysis' and 'logs' directories in the same location as the fastqs.  star_analysis/ contains the 
output from STAR + RSEM (bam, gene and isoform counts matrix, logs, etc) and logs/ contains the wrapper logs.

The remaining files from STAR + RSEM look like this, where <fastq> is the base fastq name:

<fastq>.Log.final.out		STAR final log
<fastq>.Log.out			STAR input parameters
<fastq>.Log.progress.out	STAR progress log
<fastq>.ReadsPerGene.out.tab	STAR output file
<fastq>.SJ.out.tab		STAR output file
<fastq>.sorted.bam		gapped-alignment to full genome, ready for viewing in your favorite genome browser (IGV, UCSC, etc)
<fastq>.sorted.bam.bai		bam index, required for genome browsing
<fastq>.sorted.bam.flagstat		samtools flagstat output for this bam
<fastq>.star.rsem.genes.results		RSEM expected_count, TPM, and FPKM value for genes
<fastq>.star.rsem.isoforms.results	RSEM expected_count, TPM, and FPKM value for isoforms

<fastq>._STARgenome:	this directory contains information from STAR about the genome reference
exonGeTrInfo.tab
exonInfo.tab
geneInfo.tab
sjdbInfo.txt
sjdbList.fromGTF.out.tab
sjdbList.out.tab
transcriptInfo.tab

<fastq>.star.rsem.stat:		this directory contains information from RSEM about the expectation-maximization model
<fastq>.star.rsem.cnt
<fastq>.star.rsem.model
<fastq>.star.rsem.theta

Once all jobs are complete, if desired, another R script (combine_STAR-RSEM_WGB_results.R) is run to merge all results for the sample_key.txt input.  These 
three files are created by the combine_STAR-RSEM_WGB_results.R (QGSR) for gene-level data:

 - merged_counts.genes.txt is the expected_count imported from RSEM outputs (<fastq>.star.rsem.genes.results), 
	for all samples declared in sample_key.txt.
 - merged_TPM.genes.txt is the TPM imported from RSEM outputs (<fastq>.star.rsem.genes.results), 
	for all samples declared in sample_key.txt.
 - merged_FPKM.genes.txt is the FPKM imported from RSEM outputs (<fastq>.star.rsem.genes.results), 
	for all samples declared in sample_key.txt.

There are also three isoform-level merged files for each count type (expected_count, TPM, FPKM).  There is one additional file:
	
 - merged_sample_info.txt is an extension of sample_key.txt with post-processing stats, and is the sample key for merged_counts.<gene/isoform>.txt.

These are the columns in merged_counts.genes.txt:

1	gene_id			gene_id from GTF file used in alignment and counting
2	gene_name		gene_name from GTF file used in alignment and counting
...	<fastq>			the remaining columns are counts data labeled by fastq column in merged_sample_info.txt

Similary, these are the columns in merged_counts.isoforms.txt:

1	trasncript_id		transcript_id from GTF file used in alignment and counting
2	gene_id			gene_id from GTF file used in alignment and counting
3	gene_name		gene_name from GTF file used in alignment and counting
...	<fastq>			the remaining columns are counts data labeled by fastq column in merged_sample_info.txt

The column names must be unique.  The fastq field is the only guaranteed unique field in sample_key.txt, therefore the columns are named by fastq.
Rownames must also be unique.  The "transcript_id" and "gene_id" features from the genome annotation are the guaranteed unique value for each transcript/gene, 
and are used as rownames.  A copy of the gene_id to gene_name lookup table is also provided in the fastq directory (e.g. dmel.gid2gn).

These are the columns in merged_sample_info.txt:

1	fastq		base fastq name for this sample
2	sample		sample identifier (if sample_key has a 3rd column it is inserted here before total_reads)
3	total_reads		total starting reads in this fastq before any processing
4	trimmed_reads		total reads after cutadapt trimming
5	uniq_aln_read_num	uniquely aligned reads
6	gene_det	number of genes detected (expected_count > 0) for this sample
7	isoform_det	number of isoforms detected (expected_count > 0) for this sample
8	sum_gene_count	integerized sum of expected_counts (genes) for this sample
9	pct_ercc	percent of expected_counts (genes) that are ERCC controls for this sample
10	num_ercc	number of ERCC transcripts detected (genes) out of 92 for this sample
11	yield_ratio	ratio of sum_gene_count to total_reads for this sample = the useable fraction of reads

For Janelians: Please contact QGSR staff if you have any questions about these files, the R scripts, or the python wrapper used to run the STARsolo pipeline.

A copy of this README is provided with every bulk RNA-seq (WGB) analysis run performed by QGSR.

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

Two QGSR wrapper scripts:
run_STARsolo_pipeline.py
combine_STAR-RSEM_WGB_results.R


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
