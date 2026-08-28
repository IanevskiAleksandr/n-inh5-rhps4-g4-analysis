# N-inh5–RHPS4 G4 analysis

Code and compact source data for the analysis in *G-Quadruplex stabilization via N-inh5 and RHPS4 induces cell cycle arrest and senescence in glioblastoma*.

The repository contains processed BG4 CUT&Tag and RNA contrast tables used in the study. Raw sequencing data are not duplicated here and should be cited through the accession reported with the manuscript. Repository authorship covers the code release; manuscript authorship is separate.

## Reproduction

From the repository root:

```bash
Rscript run_all.R
```

This regenerates the differential G4 analysis, motif analysis, Figure 7 panels, the assembled Figure 7 and Supplementary Figure 5. Numerical checks are applied at the end of the run. Individual steps are:

```bash
Rscript R/01_g4_analysis.R
Rscript R/02_motifs.R
Rscript R/03_figure7.R
Rscript R/04_bap1_figure.R
Rscript R/05_validate.R
```

The submitted raster is retained as `figures/Figure7_reference.jpg`; the scripted reconstruction is written to `figures/Figure7_reproduced.pdf` and `.png`. Panel source tables are written under `results/`.

The analysis distinguishes 43 factorial interaction peaks genome-wide from the 42 peaks on standard chromosomes displayed in panels B and D. Panel D uses the frozen figure definition of promoter proximity: overlap with a symmetric ±3-kb window around the gene-level TSS coordinates in `data/processed/hg38_TSS.bed.gz`. It is not an all-transcript promoter annotation.

## BAP1 external validation

Supplementary Figure 5 evaluates BAP1 using independent shRNA screens from Zyner *et al.* and pooled CRISPR screens from Olivieri *et al.* The validation script downloads the four public source files, verifies SHA-256 checksums, reconstructs the compact BAP1 tables, and asserts equality to the supplied source data.

```bash
Rscript R/validate_bap1.R
```

Downloaded files are cached under `data/raw_external/bap1/` and are excluded from version control. Set `BAP1_SOURCE_DIR` to use an existing local copy. A successful run writes `results/bap1_validation.tsv`.

Expected anchors are:

- PDS shRNA depletion: A375 genome-wide, 4/8; A375 focused, 7/7; HT1080 focused, 4/7.
- Phen-DC3 CRISPR screen: BAP1 normZ = −5.81, treatment rank 1/31, gene rank 7/17,298 and negative-direction DrugZ FDR estimate 7.87 × 10⁻⁶.
- Pyridostatin CRISPR screen: BAP1 normZ = −2.56, treatment rank 3/31, gene rank 163/17,374 and negative-direction DrugZ FDR estimate 0.552.

Negative normZ values indicate sensitization after BAP1 knockout. The Olivieri `FDRNeg` values are direction-specific DrugZ estimates supplied by the original study; they are preserved without truncation and can exceed 1. The Phen-DC3 result passes the source study's negative-direction FDR threshold, whereas the pyridostatin result is exploratory.

### Public sources

- Zyner, K.G. *et al.* Genetic interactions of G-quadruplexes in humans. *eLife* **8**, e46793 (2019), [doi:10.7554/eLife.46793](https://doi.org/10.7554/eLife.46793). Supplementary files 1 and 2 are downloaded from eLife and are licensed under CC BY 4.0 with the article.
- Olivieri, M. *et al.* A genetic map of the response to DNA damage in human cells. *Cell* **182**, 481–496.e21 (2020), [doi:10.1016/j.cell.2020.05.040](https://doi.org/10.1016/j.cell.2020.05.040). The analysis files are from [Mendeley Data version 2](https://doi.org/10.17632/gfcn2wmrpf.2), licensed under CC BY 4.0.
- DrugZ is available from the [Hart laboratory repository](https://github.com/hart-lab/drugz) under the MIT License.

Only compact derived tables are distributed in `data/external/`. The original files remain attributable to their authors and retain their source licenses.

## Requirements

R ≥ 4.5.0 with the packages recorded in `DESCRIPTION`. Bioconductor dependencies include DESeq2, GenomicRanges, BSgenome.Hsapiens.UCSC.hg38 and pqsfinder; BAP1 validation additionally uses readxl.

The computational methods and complete panel definitions are in [METHODS.md](METHODS.md) and [FIGURE_LEGENDS.md](FIGURE_LEGENDS.md).

## License

Repository code is released under the MIT License. External datasets retain their original licenses.
