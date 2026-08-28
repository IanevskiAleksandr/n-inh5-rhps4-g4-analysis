# Computational methods

## G4 CUT&Tag analysis

The consensus BG4 CUT&Tag matrix contained 4,777 peaks measured in three replicates of DMSO, RHPS4, N-inh5 and combined treatment. Integer peak counts were analysed with DESeq2 1.50.2. Size factors were fixed from mapped reads remaining after filtering and scaled to a geometric mean of one. Direct treatment effects were estimated with a four-level design using DMSO as the reference. Non-additive effects were estimated with the factorial model `~ RHPS4 + N-inh5 + RHPS4:N-inh5`. P values were adjusted by the Benjamini–Hochberg method. Gain and loss denote positive and negative direct-treatment log2 fold changes, respectively; synergistic and antagonistic denote positive and negative interaction coefficients.

The factorial analysis identified 43 peaks at FDR < 0.05 (21 positive and 22 negative). Restriction to chromosomes 1–22, X and Y removed one negative peak on `chrUn_GL000224v1`, leaving the 42 peaks displayed in the genome-wide and motif panels (21 positive and 21 negative).

## TSS signal and interaction heatmap

BG4 signal was summarized in 100-bp bins across ±10 kb of 35,143 gene-level transcription start sites. Missing bins were assigned zero and regions without signal were omitted when the matrices were generated. Panel A shows condition means and condition-specific ranked heatmaps. Panel E shows log2 BG4 signal relative to the mean DMSO signal for nine treatment libraries at 40 finalized display peaks. The BAP1 track spans chr3:52,405,664–52,408,042 (hg38); the tested peak is chr3:52,406,264–52,407,442.

## G4 motif and promoter-proximity analysis

Sequences were retrieved from hg38 and screened with pqsfinder 2.26.0 at a minimum score of 40. Interaction peaks were compared with 2,000 random regions sampled across autosomes in proportion to chromosome length and matched to the interaction-peak width distribution (`set.seed(2024)`). Motifs were detected in 20 of 42 interaction peaks and 536 of 2,000 random regions (two-sided Fisher exact P = 0.004529).

For the displayed promoter-proximity analysis, the supplied gene-level TSS coordinates were expanded symmetrically by 3 kb and reduced before overlap testing. Six of the 42 standard-chromosome interaction peaks met this operational definition; all six had positive interaction coefficients (6 of 21 positive versus 0 of 21 negative; two-sided Fisher exact P = 0.02069). This frozen definition uses gene-level boundary TSS coordinates and is not equivalent to testing ±3 kb around every annotated transcript TSS.

## G4–RNA integration

Combination-versus-DMSO G4 effects were joined to matched processed RNA contrasts by gene symbol. When more than one G4 peak mapped to a gene, the peak with the smallest finite adjusted P value and then the largest absolute effect was retained. Panel C contains 1,577 genes with both measurements. The scatter is descriptive and does not imply that a local G4 change caused the RNA response. Panel F shows selected G4-linked, p21, DREAM/cell-cycle and senescence genes. RNA log2 fold changes for each treatment were combined with a DMSO reference of zero and standardized within each gene across the four conditions.

## BAP1 external validation

BAP1 shRNA results were reconstructed from Supplementary files 1 and 2 of Zyner et al. (eLife, 2019). BAP1 CRISPR–drug responses were reconstructed from the 31-screen normalized Z-score and negative-direction DrugZ FDR matrices of Olivieri et al. (Cell, 2020; Mendeley Data version 2). `R/validate_bap1.R` downloads checksum-pinned source files and verifies every displayed value against the compact repository tables.
