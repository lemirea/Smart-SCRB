# Smart-SCRB
Single-cell RNA-seq analysis pipeline for Smart-SCRB method developed by Quantitative Genomics shared resource at Janelia Research Campus, HHMI.

This repository contains the scripts necessary to run the Smart-SCRB analysis pipeline developed by Quantitative Genomics for either DGE or WGB methods.

Requires:
- STAR v2.7.5 or higher
- RSEM v1.3.0 (for WGB)
- cutadapt v2.10 or higher
- R v4.0 or higher

The Smart-SCRB assay has several primer sets (DGE, WGB, others not released) that allow us to use a unified chemistry and protocol to create either 
stranded, barcoded, UMI-containing reads for 3' digital gene expression, or unstranded non-barcoded non-UMI-containing reads for whole gene body RNA-seq 
profiling. The unified chemistry and protocol requires only a simple change of the primer sets to switch between WGB and DGE methods at the reverse transcription
step.  Smart-SCRB works with as little as 2 pg total RNA (Drosophila whole brain Trizol prep) and as much as 5 ng total RNA, and is suitable for single-cell RNA-seq
as well as de novo transcriptome assembly and genome annotation.

Smart-SCRB DGE is used for single-cell RNA-seq as well as low input bulk RNA-seq (such as a few dozen Drosophila neurons).  DGE is useful when the researcher
needs gene-level results with a stranded assay (read 2 is the RNA strand, read 1 is barcode + UMI).

Smart-SCRB WGB is classic bulk RNA-seq using an unstranded method with no barcodes or UMIs that produces reads along the entire gene body.  This is useful when
the researcher needs isoform-level results.  Smart-SCRB WGB has been optimized to be stranded at the 5' and 3' ends and unstranded in the middle (due to 
Nextera library prep).  The Nextera read 1 sequencing primer is fused to the 5'biotin-oligo-dT(30)VN RT primer and Nextera read 2 sequencing primer is fused
to the template switch oligo at the 5' end.  Smart-SCRB WGB has been used in single-cell RNA-seq where isoforms are believed to better define cell type differences,
as well as in de novo transcriptome assembly.

Smart-SCRB was first described (without the name) in Cembrowski et al, Cell. 2018a (https://www.cell.com/cell/pdf/S0092-8674(18)30311-8.pdf) and has been used in other
studies of Drosophila, mice, zebrafish, aphids, and cell lines, some of which are listed below:





