# downstream R script for (entire sample.key of samples) Smart-SCRB WGB analysis with STAR-RSEM
# 7/27/20 by Andrew Lemire
# Quantitative Genomics
# Janelia Research Campus
# HHMI
# released under BDS 3-clause license assigned to Howard Hughes Medical Institute

if (!require("pacman", quietly = TRUE)) install.packages("pacman", repos="http://cran.r-project.org", quiet = TRUE)
pacman::p_load(ggplot2, dplyr, psych, Matrix, grid, gridExtra)

myfiles <- commandArgs(TRUE)[1]
# 1 = /path/to/sample_key.txt (just the path, not the file itself)
# this script will extract information from ./star_analysis and ./logs using the samples in myfiles[1]

skey <- read.delim(paste(myfiles[1], 'sample_key.txt', sep="/"), header=F, sep="\t", stringsAsFactors = F)
colnames(skey)[1:2] <- c('fastq', 'sample')
if (length(colnames(skey)) > 3) {print('sample_key.txt has too many columns')}
if (length(colnames(skey)) == 3) {
  colnames(skey)[3] <- 'info'
}

gid2gn <- list.files(path=myfiles[1], pattern = ".*\\.gid2gn", full.names = T)
gid2gn <- read.delim(gid2gn, header = F, sep = "\t", stringsAsFactors = F, col.names = c("gene_id", "gene_name"))
rownames(gid2gn) <- gid2gn$gene_id

merged.gene.cnt <- data.frame()
merged.gene.tpm <- data.frame()
merged.gene.fpkm <- data.frame()
merged.isoform.cnt <- data.frame()
merged.isoform.tpm <- data.frame()
merged.isoform.fpkm <-  data.frame()
merged.isoform.pct <- data.frame()
allsamples <- unique(skey$fastq)
outs <- paste(myfiles[1], '/star_analysis/', sep = '')
logs <- paste(myfiles[1], '/logs/', sep = '')
per_fastq_stats <- data.frame(skey, row.names = allsamples, stringsAsFactors = F)
for (samp in allsamples) {
  logfile <- list.files(path = logs, pattern = samp, full.names = T)[[1]]
  tot_reads <- strsplit(grep('Total read', readLines(logfile), value = T), ' ') %>%  unlist()
  tot_reads <- as.integer(gsub(',', '', tot_reads[length(tot_reads)]))
  per_fastq_stats[samp, "total_reads"] <- tot_reads
  trim_reads <- strsplit(grep('Pairs written|Reads written', readLines(logfile), value = T), ' ') %>%  unlist()
  trim_reads <- as.integer(gsub(',', '', trim_reads[length(trim_reads)-1]))
  per_fastq_stats[samp, "trimmed_reads"] <- trim_reads
  logfile <- list.files(path = outs, pattern = paste(samp, '.Log.final.out', sep = ''), full.names = T)[[1]]
  uniq_reads <- strsplit(grep('Uniquely mapped reads number', readLines(logfile), value = T), '\t') %>%  unlist()
  uniq_reads <- as.integer(uniq_reads[length(uniq_reads)])
  per_fastq_stats[samp, "uniq_aln_read_num"] <- uniq_reads
  gene.cnt <- read.delim(paste(outs, samp, '.star.rsem.genes.results', sep=''), stringsAsFactors = F, sep = "\t")
  rownames(gene.cnt) <- gene.cnt$gene_id
  per_fastq_stats[samp, "gene_det"] <- sum(gene.cnt$expected_count>0)
  iso.cnt <- read.delim(paste(outs, samp, '.star.rsem.isoforms.results', sep=''), stringsAsFactors = F, sep = "\t")
  rownames(iso.cnt) <- iso.cnt$transcript_id
  per_fastq_stats[samp, "isoform_det"] <- sum(iso.cnt$expected_count>0)
  per_fastq_stats[samp, "sum_gene_count"] <- round(sum(gene.cnt$expected_count),0)
  ercc.idx <- grepl('ERCC-0', rownames(gene.cnt))
  per_fastq_stats[samp, "pct_ercc"] <- round((sum(gene.cnt[ercc.idx,'expected_count'])/per_fastq_stats[samp, "sum_gene_count"])*100,4)
  per_fastq_stats[samp, "num_ercc"] <- sum(gene.cnt[ercc.idx,'expected_count']>0)
  per_fastq_stats[samp, "yield_ratio"] <- per_fastq_stats[samp, "sum_gene_count"]/per_fastq_stats[samp, "total_reads"]
  # now build the 6 output files

  if (length(merged.gene.cnt) == 0) {
    merged.gene.cnt <- data.frame(gene.cnt$expected_count)
    colnames(merged.gene.cnt) <- samp
    rownames(merged.gene.cnt) <- rownames(gene.cnt)
  } else {
    if (sum(rownames(gene.cnt) != rownames(merged.gene.cnt)) != 0) {
      print('rownames QC check failed')
    } else{
      merged.gene.cnt[, samp] <- gene.cnt$expected_count
    }
  }
  if (length(merged.gene.tpm) == 0) {
    merged.gene.tpm <- data.frame(gene.cnt$TPM)
    colnames(merged.gene.tpm) <- samp
    rownames(merged.gene.tpm) <- rownames(gene.cnt)
  } else {
    if (sum(rownames(gene.cnt) != rownames(merged.gene.tpm)) != 0) {
      print('rownames QC check failed')
    } else{
      merged.gene.tpm[, samp] <- gene.cnt$TPM
    }
  }
  if (length(merged.gene.fpkm) == 0) {
    merged.gene.fpkm <- data.frame(gene.cnt$FPKM)
    colnames(merged.gene.fpkm) <- samp
    rownames(merged.gene.fpkm) <- rownames(gene.cnt)
  } else {
    if (sum(rownames(gene.cnt) != rownames(merged.gene.fpkm)) != 0) {
      print('rownames QC check failed')
    } else{
      merged.gene.fpkm[, samp] <- gene.cnt$FPKM
    }
  }
  if (length(merged.isoform.cnt) == 0) {
    merged.isoform.cnt <- data.frame(iso.cnt$expected_count)
    colnames(merged.isoform.cnt) <- samp
    rownames(merged.isoform.cnt) <- rownames(iso.cnt)
  } else {
    if (sum(rownames(iso.cnt) != rownames(merged.isoform.cnt)) != 0) {
      print('rownames QC check failed')
    } else{
      merged.isoform.cnt[, samp] <- iso.cnt$expected_count
    }
  }
  if (length(merged.isoform.tpm) == 0) {
    merged.isoform.tpm <- data.frame(iso.cnt$TPM)
    colnames(merged.isoform.tpm) <- samp
    rownames(merged.isoform.tpm) <- rownames(iso.cnt)
  } else {
    if (sum(rownames(iso.cnt) != rownames(merged.isoform.cnt)) != 0) {
      print('rownames QC check failed')
    } else{
      merged.isoform.tpm[, samp] <- iso.cnt$TPM
    }
  }
  if (length(merged.isoform.fpkm) == 0) {
    merged.isoform.fpkm <- data.frame(iso.cnt$FPKM)
    colnames(merged.isoform.fpkm) <- samp
    rownames(merged.isoform.fpkm) <- rownames(iso.cnt)
  } else {
    if (sum(rownames(iso.cnt) != rownames(merged.isoform.cnt)) != 0) {
      print('rownames QC check failed')
    } else{
      merged.isoform.fpkm[, samp] <- iso.cnt$FPKM
    }
  }
  if (length(merged.isoform.pct) == 0) {
    merged.isoform.pct <- data.frame(iso.cnt$IsoPct)
    colnames(merged.isoform.pct) <- samp
    rownames(merged.isoform.pct) <- rownames(iso.cnt)
  } else {
    if (sum(rownames(iso.cnt) != rownames(merged.isoform.pct)) != 0) {
      print('rownames QC check failed')
    } else {
      merged.isoform.pct[, samp] <- iso.cnt$IsoPct
    }
  }
}

# copy WGB README
fc_okay <- file.copy('/groups/quantitativegenomics/sequencers/seq1/analysis/shared_software/STAR-RSEM_README.txt', paste(myfiles[1], '/STAR-RSEM_README.txt', sep = ''), overwrite = T, copy.mode = F)
if (fc_okay) {
  rm(fc_okay)
} else {
  print('README file copy failed.')
  rm(fc_okay)
}

# prep files
datacols <- length(allsamples)
merged.gene.cnt$gene_id <- rownames(merged.gene.cnt)
merged.gene.cnt$gene_name <- gid2gn[rownames(merged.gene.cnt), "gene_name"]
merged.gene.cnt <- merged.gene.cnt[,c(datacols+1:2, 1:datacols)]
merged.gene.tpm$gene_id <- rownames(merged.gene.tpm)
merged.gene.tpm$gene_name <- gid2gn[rownames(merged.gene.tpm), "gene_name"]
merged.gene.tpm <- merged.gene.tpm[,c(datacols+1:2, 1:datacols)]
merged.gene.fpkm$gene_id <- rownames(merged.gene.fpkm)
merged.gene.fpkm$gene_name <- gid2gn[rownames(merged.gene.fpkm), "gene_name"]
merged.gene.fpkm <- merged.gene.fpkm[,c(datacols+1:2, 1:datacols)]
merged.isoform.cnt$transcript_id <- rownames(merged.isoform.cnt)
merged.isoform.cnt$gene_id <- iso.cnt[rownames(merged.isoform.cnt), "gene_id"]
merged.isoform.cnt$gene_name <- gid2gn[merged.isoform.cnt$gene_id, "gene_name"]
merged.isoform.cnt <- merged.isoform.cnt[,c(datacols+1:3, 1:datacols)]
merged.isoform.tpm$transcript_id <- rownames(merged.isoform.tpm)
merged.isoform.tpm$gene_id <- iso.cnt[rownames(merged.isoform.tpm), "gene_id"]
merged.isoform.tpm$gene_name <- gid2gn[merged.isoform.tpm$gene_id, "gene_name"]
merged.isoform.tpm <- merged.isoform.tpm[,c(datacols+1:3, 1:datacols)]
merged.isoform.fpkm$transcript_id <- rownames(merged.isoform.fpkm)
merged.isoform.fpkm$gene_id <- iso.cnt[rownames(merged.isoform.fpkm), "gene_id"]
merged.isoform.fpkm$gene_name <- gid2gn[merged.isoform.fpkm$gene_id, "gene_name"]
merged.isoform.fpkm <- merged.isoform.fpkm[,c(datacols+1:3, 1:datacols)]
merged.isoform.pct$transcript_id <- rownames(merged.isoform.pct)
merged.isoform.pct$gene_id <- iso.cnt[rownames(merged.isoform.pct), "gene_id"]
merged.isoform.pct$gene_name <- gid2gn[merged.isoform.pct$gene_id, "gene_name"]
merged.isoform.pct <- merged.isoform.pct[,c(datacols+1:3, 1:datacols)]
# write files and cleanup
write.table(per_fastq_stats, paste(myfiles[1], "/merged_sample_info.txt", sep=""), col.names = T, row.names = F, quote=F, sep="\t")

write.table(merged.gene.cnt, paste(myfiles[1], "/merged_counts.genes.txt", sep=""), col.names = T, row.names = F, quote=F, sep="\t")
write.table(merged.gene.tpm, paste(myfiles[1], "/merged_TPM.genes.txt", sep=""), col.names = T, row.names = F, quote=F, sep="\t")
write.table(merged.gene.fpkm, paste(myfiles[1], "/merged_FPKM.genes.txt", sep=""), col.names = T, row.names = F, quote=F, sep="\t")

write.table(merged.isoform.cnt, paste(myfiles[1], "/merged_counts.isoforms.txt", sep=""), col.names = T, row.names = F, quote=F, sep="\t")
write.table(merged.isoform.tpm, paste(myfiles[1], "/merged_TPM.isoforms.txt", sep=""), col.names = T, row.names = F, quote=F, sep="\t")
write.table(merged.isoform.fpkm, paste(myfiles[1], "/merged_FPKM.isoforms.txt", sep=""), col.names = T, row.names = F, quote=F, sep="\t")
write.table(merged.isoform.pct, paste(myfiles[1], "/merged_IsoPct.isoforms.txt", sep=""), col.names = T, row.names = F, quote=F, sep="\t")
# QC plots
# make plots for QC
# trim_ratio is a rough approximation of the use-able fraction of reads per fastq: the sum of all gene counts divided by the total number of reads in the fastq
per_fastq_stats$trim_ratio <- per_fastq_stats$trimmed_reads/per_fastq_stats$total_reads

p1 <- ggplot(per_fastq_stats, aes(x=fastq, y=pct_ercc)) + geom_segment(aes(x = fastq, xend = fastq, y=0, yend = pct_ercc), color = 'gray', size=1) + geom_point(size = 4, color = '#69b3a2') + theme_bw() + theme(legend.position = 'none') + xlab('') + geom_hline(yintercept = 5, color = 'red') + ylim(0, 100) + coord_flip() + theme(panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank(), legend.position="none") + ggtitle('Percent ERCC', subtitle = 'per fastq')
p2 <- ggplot(per_fastq_stats, aes(x=fastq, y=gene_det)) + geom_segment(aes(x = fastq, xend = fastq, y=0, yend = gene_det), color = 'gray', size=1) + geom_point(size = 4, color = '#69b3a2') + theme_bw() + theme(legend.position = 'none') + xlab('') + coord_flip() + theme(panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank(), legend.position="none") + ggtitle('Gene detection', subtitle = 'per fastq')
p3 <- ggplot(per_fastq_stats, aes(x=fastq, y=trim_ratio)) + geom_segment(aes(x = fastq, xend = fastq, y=0, yend = trim_ratio), color = 'gray', size=1) + geom_point(size = 4, color = '#69b3a2') + theme_bw() + theme(legend.position = 'none') + xlab('') + geom_hline(yintercept = 0.85, color = 'red') + ylim(0, 1) + coord_flip() + theme(panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank(), legend.position="none") + ggtitle('Ratio of trimmed : total reads', subtitle = 'per fastq')
p4 <- ggplot(per_fastq_stats, aes(x=fastq, y=yield_ratio)) + geom_segment(aes(x = fastq, xend = fastq, y=0, yend = yield_ratio), color = 'gray', size=1) + geom_point(size = 4, color = '#69b3a2') + theme_bw() + theme(legend.position = 'none') + xlab('') + ylim(0, 1) + coord_flip() + theme(panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank(), legend.position="none") + ggtitle('Yield ratio (fraction of useable reads)', subtitle = 'per fastq')
pdf(paste(myfiles[1], '/QC_plots.pdf', sep=''), width=12, height=8)
grid.arrange(p1, p3, p2, p4)
dev.off()
