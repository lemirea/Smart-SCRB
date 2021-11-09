# downstream R script for (entire platemap of samples) Smart-SCRB analysis with STARsolo
# 7/27/20 by Andrew Lemire
# Quantitative Genomics
# Janelia Research Campus
# HHMI
# released under BDS 3-clause license assigned to Howard Hughes Medical Institute

if (!require("pacman", quietly = TRUE)) install.packages("pacman", repos="http://cran.r-project.org", quiet = TRUE)
pacman::p_load(ggplot2, dplyr, psych, Matrix, grid, gridExtra)

myfiles <- commandArgs(TRUE)[1]
# 1 = /path/to/platemap.tab (just the path, not the file itself)
# this script will extract information from ./star_analysis and ./logs using the samples in myfiles[1]

# define order of rows/columns in plate format (column-wise, 8-channel pipet)
mycols <- sprintf("%02d",c(rep(1,8),rep(2,8),rep(3,8),rep(4,8),rep(5,8),rep(6,8),rep(7,8),rep(8,8),rep(9,8),rep(10,8),rep(11,8),rep(12,8)))
myrows <- rep(c("A","B","C","D","E","F","G","H"),12)
myplate <- paste0(myrows,mycols)
rm(mycols, myrows)

# read platemap.tab file (must be in the directory passed to this script, not as a symlink)
pmap <- read.delim(paste(myfiles[1], 'platemap.tab', sep="/"), header=F, sep="\t", col.names = c("inst","run_num","fastq","well","q_plate","sample"), stringsAsFactors = F)
pmap$qp_well_fq <- paste(pmap$q_plate, pmap$well, pmap$fastq, sep="_")
rownames(pmap) <- pmap$qp_well_fq
pmap$well <- factor(pmap$well,levels=myplate, ordered = TRUE)
pmap$q_plate <- factor(pmap$q_plate,levels=c("Q5","Q6","Q7","Q8"), ordered = TRUE)
pmap$sample <- factor(pmap$sample)
pmap$fastq <- as.character(pmap$fastq)

# get all fastq names and define data locations
allsamples <- unique(pmap$fastq)
outs <- paste(myfiles[1], '/star_analysis/', sep = '')
logs <- paste(myfiles[1], '/logs/', sep = '')

# loop through all samples to extract per-sample pstats (from STARsolo_DGE_downstream.R)
# and to extract info from logs to complete platemap metadata (reads, trimmed reads, etc)
merged_count <- data.frame()
merged_pstat <- data.frame()
per_fastq_stats <- data.frame(row.names = allsamples)

for (samp in allsamples){
  if (dir.exists(paste(outs, samp, '.Solo.out/Gene/raw', sep = ''))) {
    countMethod <- 'Gene'
  } else if (dir.exists(paste(outs, samp, '.Solo.out/GeneFull/raw', sep = ''))) {
    countMethod <- 'GeneFull'
  }
  rawData <- paste(outs, samp, '.Solo.out/', countMethod, '/raw/', sep = '')
  pmap.sample <- read.delim(paste(outs, samp, '.pstats.txt', sep = ''), header = T, sep = "\t", stringsAsFactors = F)
  rownames(pmap.sample) <- pmap.sample$qp_well_fq
  stargene <- read.delim(paste(rawData, 'features.tsv', sep = ''), header = F, sep = "\t", col.names = c("gene_id", "gene_name", "countMethod"), stringsAsFactors = F, row.names = 1)
  stargene$gene_id <- rownames(stargene)
  starbcs <- read.delim(paste(rawData, 'barcodes.tsv', sep = ''), header=F, sep="\t", col.names = c("barcode"), stringsAsFactors = F)
  starcnt <- readMM(paste(rawData, 'matrix.mtx', sep = ''))
  cnts <- as.matrix(starcnt)
  cnts <- as.data.frame(cnts)
  rownames(cnts) <- stargene$gene_id
  colnames(cnts) <- starbcs$barcode
  cnts <- as.data.frame(cnts[,pmap.sample$barcode], row.names = stargene$gene_id)
  colnames(cnts) <- pmap.sample$qp_well_fq
  # extract fastq metadata from various logs
  logfile <- list.files(path = logs, pattern = samp, full.names = T)[[1]]
  tot_reads <- strsplit(grep('Total read', readLines(logfile), value = T), ' ') %>%  unlist()
  tot_reads <- as.integer(gsub(',', '', tot_reads[length(tot_reads)]))
  pmap.sample$yield_ratio <- lapply(pmap.sample$sum_count, function(x) x/tot_reads) %>% unlist()
  per_fastq_stats[samp, "total_reads"] <- tot_reads
  trim_reads <- strsplit(grep('Pairs written', readLines(logfile), value = T), ' ') %>%  unlist()
  trim_reads <- as.integer(gsub(',', '', trim_reads[length(trim_reads)-1]))
  per_fastq_stats[samp, "trimmed_reads"] <- trim_reads
  data_summary_csv <- paste(outs, samp, '.Solo.out/', countMethod, '/Summary.csv',  sep='')
  sumStats <- read.csv(data_summary_csv, header = F, col.names = c('variable', 'value'))
  sumStats$variable <- gsub(pattern = ' ', '_', sumStats$variable)
  rownames(sumStats) <-  sumStats$variable
  per_fastq_stats[samp, "q30_bc_umi"] <- sumStats["Q30_Bases_in_CB+UMI", 2]
  per_fastq_stats[samp, "q30_cDNA"] <- sumStats["Q30_Bases_in_RNA_read", 2]
  logfile <- list.files(path = outs, pattern = paste(samp, '.Log.final.out', sep = ''), full.names = T)[[1]]
  uniq_reads <- strsplit(grep('Uniquely mapped reads number', readLines(logfile), value = T), '\t') %>%  unlist()
  uniq_reads <- as.integer(uniq_reads[length(uniq_reads)])
  per_fastq_stats[samp, "uniq_aln_read_num"] <- uniq_reads
  feats <-  paste(outs, samp, '.Solo.out/', countMethod, '/Features.stats', sep='')
  feat_data <- read.table(feats, header = F, col.names = c("variable", "value"), stringsAsFactors = F)
  rownames(feat_data) <-  feat_data$variable
  for (n in feat_data$variable) {
    per_fastq_stats[samp, n] <- feat_data[n, 2]
  }
  if (countMethod == 'Gene') {
    per_fastq_stats[samp, "map2txome"] <- sumStats["Reads_Mapped_to_Transcriptome:_Unique_Genes", 2]
  } else {
    per_fastq_stats[samp, "map2txome"] <- sumStats["Reads_Mapped_to_Transcriptome:_Unique_GeneFulls", 2]
  }
  per_fastq_stats[samp, "saturation"] <- sumStats["Sequencing_Saturation",2]
  per_fastq_stats[samp, 'sum_all_count'] <- sum(pmap.sample$sum_count)
  per_fastq_stats[samp, 'yield_ratio'] <- round(per_fastq_stats[samp, 'sum_all_count']/per_fastq_stats[samp, 'total_reads'], 6)
  # merge counts and metadata
  if (sum(colSums(cnts) != pmap.sample$sum_count) != 0) {print('colSum QC check failed')}
  if (length(merged_count) == 0) {
    merged_count <- cnts
  } else {
    if (sum(rownames(cnts) != rownames(merged_count)) != 0) {
      print('rownames QC check failed')
    } else{
      merged_count <- cbind(merged_count, cnts)
    }
  }
  if (length(merged_pstat) == 0) {
    merged_pstat <- pmap.sample
  } else {
    if (sum(colnames(pmap.sample) != colnames(merged_pstat)) != 0) {
      print('platemap stats colnames QC check failed')
    } else {
      merged_pstat <- rbind(merged_pstat, pmap.sample)
    }
  }
}

merged_count$gene_name <- stargene$gene_name
merged_count$gene_id <- rownames(merged_count)
gn_col <- length(colnames(merged_count))
data_col <- length(colnames(merged_count)) - 2
merged_count <- merged_count[,c(gn_col, gn_col-1, 1:data_col)]

merged_pstat <- merged_pstat[,c(8, 3:5, 7, 6, 1:2, 9:13)]

per_fastq_stats$fastq <- rownames(per_fastq_stats)
sn_col <- length(colnames(per_fastq_stats))
per_fastq_stats <- per_fastq_stats[,c(sn_col, 1:(sn_col-1))]

fc_okay <- file.copy('/groups/quantitativegenomics/sequencers/seq1/analysis/shared_software/STARsolo_README.txt', paste(myfiles[1], '/STARsolo_README.txt', sep = ''), overwrite = T, copy.mode = F)
if (fc_okay) {
  rm(fc_okay)
} else {
  print('README file copy failed.')
  rm(fc_okay)
}

# write files and cleanup
write.table(per_fastq_stats, paste(myfiles[1], "/merged_fastq_info.txt", sep=""), col.names = T, row.names = F, quote=F, sep="\t")
write.table(merged_count, paste(myfiles[1], "/merged_counts.txt", sep=""), col.names = T, row.names = F, quote=F, sep="\t")
write.table(merged_pstat, paste(myfiles[1], "/merged_sample_info.txt", sep=""), col.names = T, row.names = F, quote=F, sep="\t")

# trim_ratio is a rough approximation of the use-able fraction of reads per fastq: the sum of all gene counts divided by the total number of reads in the fastq
per_fastq_stats$trim_ratio <- per_fastq_stats$trimmed_reads/per_fastq_stats$total_reads
# make plots for QC
p1 <- ggplot(merged_pstat, aes(x=fastq, y=pct_ercc)) + geom_jitter(alpha = 0.35, width = 0.1, size = 1) + geom_violin(alpha = 0.7, aes(fill = fastq)) + xlab('') + theme_bw() + theme(axis.text.x = element_text(angle = 90, size = 6)) + theme(legend.position = 'none') + geom_hline(yintercept = 5, color = 'red') + ggtitle('Percent ERCC', subtitle = 'per sample')
p2 <- ggplot(merged_pstat, aes(x=fastq, y=gene_det)) + geom_jitter(alpha = 0.35, width = 0.1, size = 1) + geom_violin(alpha = 0.7, aes(fill = fastq)) + xlab('') + theme_bw() + theme(axis.text.x = element_text(angle = 90, size = 6)) + theme(legend.position = 'none') + ggtitle('Gene detection', subtitle = 'per sample')
p3 <- ggplot(per_fastq_stats, aes(x=fastq, y=trim_ratio)) + geom_segment(aes(x = fastq, xend = fastq, y=0, yend = trim_ratio), color = 'gray', size=1 ) + geom_point(size = 4, color = '#69b3a2') + theme_bw() + theme(legend.position = 'none') + xlab('') + geom_hline(yintercept = 0.85, color = 'red') + ylim(0, 1) + coord_flip() + theme(panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank(), legend.position="none") + ggtitle('Ratio of trimmed : total reads', subtitle = 'per fastq')
p4 <- ggplot(per_fastq_stats, aes(x=fastq, y=yield_ratio)) + geom_segment(aes(x = fastq, xend = fastq, y=0, yend = yield_ratio), color = 'gray', size=1) + geom_point(size = 4, color = '#69b3a2') + theme_bw() + theme(legend.position = 'none') + xlab('') + ylim(0, 1) + coord_flip() + theme(panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank(), legend.position="none") + ggtitle('Yield ratio (fraction of useable reads)', subtitle = 'per fastq')
pdf(paste(myfiles[1], '/QC_plots.pdf', sep=''), width=12, height=8)
grid.arrange(p1, p3, p2, p4)
dev.off()

rm(pmap.sample, tot_reads, trim_reads, uniq_reads, samp, rawData, logfile, data_summary_csv, feats, feat_data, sumStats, starbcs, starcnt, stargene, cnts, gn_col, data_col, sn_col)
rm(allsamples, logs, outs, myfiles, myplate, n, pmap, merged_count, merged_pstat, per_fastq_stats)

# sessionInfo()


