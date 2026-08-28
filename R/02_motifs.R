suppressPackageStartupMessages({
  library(GenomicRanges)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(pqsfinder)
})

source("R/functions.R")

analysis <- readRDS("results/g4_analysis.rds")
genome <- BSgenome.Hsapiens.UCSC.hg38

interaction <- analysis$interaction_standard
interaction_gr <- GRanges(interaction$CHR, IRanges(interaction$START, interaction$END))
all_gr <- GRanges(analysis$coordinates$CHR, IRanges(analysis$coordinates$START, analysis$coordinates$END))
all_gr <- unique(keepStandardChromosomes(all_gr, pruning.mode = "coarse"))

tss <- read.delim(gzfile("data/processed/hg38_TSS.bed.gz"), header = FALSE)
colnames(tss) <- c("chr", "start", "end", "name", "score", "strand")
tss_gr <- GRanges(tss$chr, IRanges(tss$start + 1L, tss$end), strand = tss$strand)
promoters <- reduce(suppressWarnings(trim(GRanges(seqnames(tss_gr), IRanges(pmax(1L, start(tss_gr) - 3000L), end(tss_gr) + 3000L)))))
interaction$is_promoter_proximal <- overlapsAny(interaction_gr, promoters)
promoter_gr <- interaction_gr[interaction$is_promoter_proximal]

set.seed(2024)
chromosome_lengths <- seqlengths(genome)[paste0("chr", 1:22)]
random_widths <- sample(width(interaction_gr), 2000L, replace = TRUE)
random_chromosomes <- sample(names(chromosome_lengths), 2000L, replace = TRUE, prob = chromosome_lengths)
random_starts <- mapply(function(chr, width) sample.int(chromosome_lengths[[chr]] - width - 1L, 1L), random_chromosomes, random_widths)
random_gr <- GRanges(random_chromosomes, IRanges(random_starts, random_starts + random_widths - 1L))

has_motif <- function(gr) {
  sequences <- getSeq(genome, gr)
  sink(nullfile())
  on.exit(sink())
  vapply(seq_along(sequences), function(i) length(pqsfinder(sequences[[i]], min_score = 40)) > 0L, logical(1))
}

motifs <- list(
  `Promoter interaction peaks` = has_motif(promoter_gr),
  `All interaction peaks` = has_motif(interaction_gr),
  `All G4 peaks` = has_motif(all_gr),
  `Random regions` = has_motif(random_gr)
)

motif_summary <- do.call(rbind, lapply(names(motifs), function(name) {
  x <- motifs[[name]]
  data.frame(set = name, n = length(x), motif_positive = sum(x), percent = 100 * mean(x))
}))

motif_test <- fisher.test(matrix(c(sum(motifs[["All interaction peaks"]]), sum(!motifs[["All interaction peaks"]]), sum(motifs[["Random regions"]]), sum(!motifs[["Random regions"]])), nrow = 2, byrow = TRUE))
promoter_table <- table(
  location = factor(ifelse(interaction$is_promoter_proximal, "Promoter", "Non-promoter"), levels = c("Non-promoter", "Promoter")),
  direction = factor(interaction$direction, levels = c("Synergistic", "Antagonistic"))
)
promoter_test <- fisher.test(promoter_table)

interaction$motif_pqsfinder40 <- motifs[["All interaction peaks"]]
saveRDS(
  list(
    interaction = interaction,
    motif_summary = motif_summary,
    motif_test = motif_test,
    promoter_table = promoter_table,
    promoter_test = promoter_test
  ),
  "results/motif_analysis.rds"
)
write_tsv(motif_summary, "results/motif_summary.tsv")
write_tsv(as.data.frame.matrix(promoter_table), "results/promoter_direction_table.tsv")
write_tsv(interaction, "results/interaction_peak_motifs.tsv")
