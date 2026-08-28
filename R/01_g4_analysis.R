suppressPackageStartupMessages({
  library(DESeq2)
  library(GenomicRanges)
})

source("R/functions.R")

count_table <- read.delim(gzfile("data/processed/g4_consensus_counts.tsv.gz"), check.names = FALSE)
sample_table <- read.delim("data/processed/samples.tsv", check.names = FALSE)
coordinates <- count_table[, c("CHR", "START", "END")]
counts <- as.matrix(count_table[, sample_table$sample, drop = FALSE])
storage.mode(counts) <- "integer"
rownames(counts) <- peak_key(coordinates$CHR, coordinates$START, coordinates$END)
stopifnot(identical(colnames(counts), sample_table$sample), nrow(counts) == 4777L)

pair_metadata <- data.frame(
  condition = factor(sample_table$condition, levels = c("NT", "RHPS4", "Inh5", "Combo")),
  row.names = sample_table$sample
)
pair_dds <- DESeqDataSetFromMatrix(counts, pair_metadata, ~condition)
sizeFactors(pair_dds) <- sample_table$filtered_library_reads / geometric_mean(sample_table$filtered_library_reads)
pair_dds <- DESeq(pair_dds, quiet = TRUE)

extract_result <- function(dds, contrast) {
  x <- as.data.frame(results(dds, contrast = contrast))
  cbind(coordinates, x)
}

pairwise <- list(
  RHPS4 = extract_result(pair_dds, c("condition", "RHPS4", "NT")),
  `N-inh5` = extract_result(pair_dds, c("condition", "Inh5", "NT")),
  Combo = extract_result(pair_dds, c("condition", "Combo", "NT"))
)

factorial_metadata <- data.frame(
  RHPS4 = factor(sample_table$RHPS4, levels = c("no", "yes")),
  Inh5 = factor(sample_table$Inh5, levels = c("no", "yes")),
  row.names = sample_table$sample
)
factorial_dds <- DESeqDataSetFromMatrix(counts, factorial_metadata, ~RHPS4 + Inh5 + RHPS4:Inh5)
sizeFactors(factorial_dds) <- sample_table$filtered_library_reads / geometric_mean(sample_table$filtered_library_reads)
factorial_dds <- DESeq(factorial_dds, quiet = TRUE)
interaction_name <- grep("RHPS4yes.*Inh5yes", resultsNames(factorial_dds), value = TRUE)
stopifnot(length(interaction_name) == 1L)
interaction <- cbind(coordinates, as.data.frame(results(factorial_dds, name = interaction_name)))

annotations <- read.csv(gzfile("data/processed/g4_peak_annotations.csv.gz"), check.names = FALSE)
annotation_key <- peak_key(annotations$seqnames, annotations$start, annotations$end)
for (name in names(pairwise)) {
  k <- peak_key(pairwise[[name]]$CHR, pairwise[[name]]$START, pairwise[[name]]$END)
  pairwise[[name]]$SYMBOL <- annotations$SYMBOL[match(k, annotation_key)]
  pairwise[[name]]$location <- annotations$loc_cat[match(k, annotation_key)]
}
interaction_key <- peak_key(interaction$CHR, interaction$START, interaction$END)
interaction$SYMBOL <- annotations$SYMBOL[match(interaction_key, annotation_key)]
interaction$location <- annotations$loc_cat[match(interaction_key, annotation_key)]

interaction_fdr05 <- interaction[!is.na(interaction$padj) & interaction$padj < 0.05, ]
interaction_standard <- interaction_fdr05[interaction_fdr05$CHR %in% standard_chromosomes, ]
interaction_standard <- interaction_standard[!duplicated(peak_key(interaction_standard$CHR, interaction_standard$START, interaction_standard$END)), ]
interaction_standard$direction <- ifelse(interaction_standard$log2FoldChange > 0, "Synergistic", "Antagonistic")

saveRDS(
  list(
    coordinates = coordinates,
    counts = counts,
    sample_table = sample_table,
    pairwise = pairwise,
    interaction = interaction,
    interaction_standard = interaction_standard
  ),
  "results/g4_analysis.rds"
)

for (name in names(pairwise)) write_tsv(pairwise[[name]], file.path("results", paste0("pairwise_", gsub("-", "", name), "_vs_NT.tsv")))
write_tsv(interaction, "results/factorial_interaction.tsv")
write_tsv(interaction_standard, "results/factorial_interaction_FDR05_standard.tsv")

summary_table <- do.call(rbind, lapply(names(pairwise), function(name) {
  x <- pairwise[[name]]
  x <- x[!is.na(x$padj) & x$padj < 0.05, ]
  data.frame(contrast = paste(name, "vs NT"), significant = nrow(x), gained = sum(x$log2FoldChange > 0), lost = sum(x$log2FoldChange < 0))
}))
summary_table <- rbind(
  summary_table,
  data.frame(
    contrast = "RHPS4:N-inh5 interaction",
    significant = nrow(interaction_fdr05),
    gained = sum(interaction_fdr05$log2FoldChange > 0),
    lost = sum(interaction_fdr05$log2FoldChange < 0)
  )
)
write_tsv(summary_table, "results/g4_summary.tsv")
