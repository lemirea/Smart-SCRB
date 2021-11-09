# Smart-SCRB 
Single-cell or bulk RNA-seq analysis pipeline for the Smart-SCRB method developed by the Quantitative Genomics shared resource at Janelia Research Campus, HHMI.

This repository contains the scripts necessary to run the Smart-SCRB analysis pipeline developed by Quantitative Genomics for either DGE or WGB methods, which are
described in detail in their respective STAR\*README\*.txt files.

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

de Rus Jacquet A, Tancredi JL, Lemire AL, DeSantis MC, Li WP, O'Shea EK. The LRRK2 G2019S mutation alters astrocyte-to-neuron communication via extracellular vesicles and induces neuron atrophy in a human iPSC-derived model of Parkinson's disease. Elife. 2021 Sep 30;10. doi: 10.7554/eLife.73062. PubMed PMID: 34590578; PubMed Central PMCID: PMC8514240. 

Erwin SR, Bristow BN, Sullivan KE, Kendrick RM, Marriott B, Wang L, Clements J, Lemire AL, Jackson J, Cembrowski MS. Spatially patterned excitatory neuron subtypes and projections of the claustrum. Elife. 2021 Aug 16;10. doi: 10.7554/eLife.68967. PubMed PMID: 34397382; PubMed Central PMCID: PMC8367382. 

Korgaonkar A, Han C, Lemire AL, Siwanowicz I, Bennouna D, Kopec RE, Andolfatto P, Shigenobu S, Stern DL. A novel family of secreted insect proteins linked to plant gall development. Curr Biol. 2021 May 10;31(9):2038. doi: 10.1016/j.cub.2021.03.001. PubMed PMID: 33974861. 

Xu S, Yang H, Menon V, Lemire AL, Wang L, Henry FE, Turaga SC, Sternson SM. Behavioral state coding by molecularly defined paraventricular hypothalamic cell type ensembles. Science. 2020 Oct 16;370(6514). doi: 10.1126/science.abb2494. PubMed PMID: 33060330. 

O'Leary TP, Sullivan KE, Wang L, Clements J, Lemire AL, Cembrowski MS. Extensive and spatially variable within-cell-type heterogeneity across the basolateral amygdala. Elife. 2020 Sep 1;9. doi: 10.7554/eLife.59003. PubMed PMID: 32869744; PubMed Central PMCID: PMC7486123. 

Aso Y, Ray RP, Long X, Bushey D, Cichewicz K, Ngo TT, Sharp B, Christoforou C, Hu A, Lemire AL, Tillberg P, Hirsh J, Litwin-Kumar A, Rubin GM. Nitric oxide acts as a cotransmitter in a subset of dopaminergic neurons to diversify memory dynamics. Elife. 2019 Nov 14;8. doi: 10.7554/eLife.49257. PubMed PMID: 31724947; PubMed Central PMCID: PMC6948953. 

Phillips JW, Schulmann A, Hara E, Winnubst J, Liu C, Valakh V, Wang L, Shields BC, Korff W, Chandrashekar J, Lemire AL, Mensh B, Dudman JT, Nelson SB, Hantman AW. A repeated molecular architecture across thalamic pathways. Nat Neurosci. 2019 Nov;22(11):1925-1935. doi: 10.1038/s41593-019-0483-3. Epub 2019 Sep 16. PubMed PMID: 31527803; PubMed Central PMCID: PMC6819258. 

Cembrowski MS, Wang L, Lemire AL, Copeland M, DiLisio SF, Clements J, Spruston N. The subiculum is a patchwork of discrete subregions. Elife. 2018 Oct 30;7. doi: 10.7554/eLife.37701. PubMed PMID: 30375971; PubMed Central PMCID: PMC6226292. 








