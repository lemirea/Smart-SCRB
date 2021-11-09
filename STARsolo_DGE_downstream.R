# downstream R script for (single sample) Smart-SCRB analysis with STARsolo
# 7/27/20 by Andrew Lemire
# Quantitative Genomics
# Janelia Research Campus
# HHMI
# released under BDS 3-clause license assigned to Howard Hughes Medical Institute

if (!require("pacman", quietly = TRUE)) install.packages("pacman", repos="http://cran.r-project.org", quiet = TRUE)
pacman::p_load(ggplot2, dplyr, psych, Matrix)


myfiles <- commandArgs(TRUE)[1:7]
# 1 = Q.index
# 2 = platemap.tab
# 3 = STARsolo features.tsv (list of gene_id (rownames for matrix), gene_name, countMethod)
# 4 = STARsolo barcodes.tsv (list of sample barcodes (colnames for matrix))
# 5 = STARsolo matrix.mtx results file
# 6 = root sample name (* in *_1.fastq.gz)
# 7 = ./star_analysis output dir

# qgsr files
qidx <- read.delim(myfiles[1], header=F, sep="\t", col.names = c("barcode","q_plate","well"), stringsAsFactors = F)
pmap <- read.delim(myfiles[2], header=F, sep="\t", col.names = c("inst","run_num","fastq","well","q_plate","sample"), stringsAsFactors = F)

# STARsolo files
stargene <- read.delim(myfiles[3], header = F, sep = "\t", col.names = c("gene_id", "gene_name", "countMethod"), stringsAsFactors = F, row.names = 1)
stargene$gene_id <- rownames(stargene)
starbcs <- read.delim(myfiles[4], header=F, sep="\t", col.names = c("barcode"), stringsAsFactors = F)
starcnt <- readMM(myfiles[5])

# process qgsr files
qidx <- qidx %>% group_by(q_plate, well) %>% mutate("qp_well" = paste(q_plate, well, sep="_"))
pmap <- pmap %>% group_by(q_plate, well, fastq) %>% mutate("qp_well" = paste(q_plate, well, sep="_"))
pmap <- left_join(x = pmap, y = qidx, by = "qp_well", suffix = c("", ".y"))
pmap[,9:10] <- NULL
pmap$qp_well <- NULL
pmap$qp_well_fq <- paste(pmap$q_plate, pmap$well, pmap$fastq, sep="_")

# process STARsolo results
g2 <- as.matrix(starcnt)
g2 <- as.data.frame(g2)
rownames(g2) <- stargene$gene_id
colnames(g2) <- starbcs$barcode
g2$gene_name <- stargene$gene_name

# split results by expected and found
bc_expected <- pmap$barcode[pmap$fastq==myfiles[6]]
qpfq_expected <- pmap$qp_well_fq[pmap$fastq==myfiles[6]]
bc_found <- bc_expected %in% starbcs$barcode

# i want to force R to output all barcodes expected from the platemap regardless of what STARsolo found
# because in the case of lazy downsampling some wells may be missing. report as NA
expect <- as.data.frame(matrix(,nrow = length(rownames(g2)), ncol = length(bc_expected)))
rownames(expect) <- rownames(g2)
colnames(expect) <- qpfq_expected
expect[,!bc_found] <- NaN
expect[,bc_found] <- g2[,bc_expected[bc_found]]

# calculate stats
pmap.expect <- pmap[pmap$fastq==myfiles[6],]
pmap.expect$gene_det <- colSums(expect>0)
pmap.expect$sum_count <- colSums(expect)
ercc.idx <- grepl('ERCC-0', rownames(expect))
pmap.expect$pct_ercc <- round((colSums(data.frame(expect[ercc.idx,]))/pmap.expect$sum_count)*100,4)
pmap.expect$num_ercc <- colSums(data.frame(expect[ercc.idx,]>0))

# write .pstats.txt file
write.table(pmap.expect, paste(myfiles[7], myfiles[6], ".pstats.txt", sep=""), col.names = T, row.names = F, quote=F, sep="\t")

# re-order columns in expect (counts matrix for expected barcodes)
expect$gene_id <- rownames(expect)
expect$gene_name <- g2$gene_name
gn_col <- length(colnames(expect))
gid_col <- gn_col - 1
dat_col <- gn_col - 2
expect <- expect[,c(gid_col, gn_col, 1:dat_col)]

# write the expected results to file.
write.table(expect, paste(myfiles[7], myfiles[6], ".counts.txt", sep=""), col.names = T, row.names = F, quote=F, sep="\t")


# what about the barcodes STARsolo reports that weren't in our platemap?
# this is highly likely to be lab contamination or a labeling error in platemap.tab
bc_other <- setdiff(starbcs$barcode, bc_expected)
if (length(bc_other) != 0) {
samples_other <- qidx$qp_well[match(bc_other, qidx$barcode)]
other <- as.data.frame(matrix(,nrow = length(rownames(g2)), ncol = length(bc_other)))
rownames(other) <- rownames(g2)
colnames(other) <- samples_other
other[,samples_other] <- g2[,bc_other]
# now re-order "other" so that the column order loosely represents qp_well
neworder_other <- qidx$qp_well[sort(match(samples_other, qidx$qp_well))]
other <- other[,neworder_other]
sumCounts <- colSums(other)
# filter to remove contam bcs with <=100 counts
other <- other[,sumCounts>100]
if (length(colnames(other)>0)) {
other$gene_id <- rownames(other)
other$gene_name <- g2$gene_name
gn_col <- length(colnames(other))
gid_col <- gn_col - 1
dat_col <- gn_col - 2
other <- other[,c(gid_col, gn_col, 1:dat_col)]
print(paste0("Found reads from other Q.index wells:", collapse=" "))
print(colSums(other[,3:gn_col]))
write.table(other, paste(myfiles[7], myfiles[6], ".contam.txt", sep=""), col.names = T, row.names = T, quote=F, sep="\t")
} else {
  print(paste0("No significant contamination detected.  No contam.txt file created.", collapse=""))
}
}

# log sessionInfo
sessionInfo()



