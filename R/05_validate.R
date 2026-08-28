source("R/functions.R")

analysis <- readRDS("results/g4_analysis.rds")
motif <- readRDS("results/motif_analysis.rds")
summary_table <- read.delim("results/g4_summary.tsv")

stopifnot(
  nrow(analysis$coordinates) == 4777L,
  identical(summary_table$significant, c(176L, 195L, 228L, 43L)),
  identical(summary_table$gained, c(108L, 106L, 125L, 21L)),
  identical(summary_table$lost, c(68L, 89L, 103L, 22L)),
  nrow(analysis$interaction_standard) == 42L,
  sum(analysis$interaction_standard$log2FoldChange > 0) == 21L,
  sum(analysis$interaction_standard$log2FoldChange < 0) == 21L,
  identical(motif$motif_summary$n, c(6L, 42L, 4754L, 2000L)),
  identical(motif$motif_summary$motif_positive, c(5L, 20L, 1951L, 536L)),
  all(motif$promoter_table == matrix(c(15L, 6L, 21L, 0L), nrow = 2L)),
  abs(motif$motif_test$p.value - 0.004529359) < 1e-8,
  abs(motif$promoter_test$p.value - 0.02068861) < 1e-8,
  nrow(read.delim("results/figure7C_source_data.tsv")) == 1577L,
  nrow(read.delim("data/processed/interaction_peak_signal.tsv.gz")) == 40L
)

bap1 <- analysis$interaction[analysis$interaction$CHR == "chr3" & analysis$interaction$START == 52406264L & analysis$interaction$END == 52407442L, ]
stopifnot(
  nrow(bap1) == 1L,
  abs(bap1$log2FoldChange - 16.62547) < 1e-4,
  abs(bap1$padj - 0.007127501) < 1e-8
)

cat("Figure 7 analytical checks: PASS\n")
