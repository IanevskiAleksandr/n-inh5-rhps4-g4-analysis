suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(scales)
  library(circlize)
  library(GenomicRanges)
  library(BSgenome.Hsapiens.UCSC.hg38)
})

source("R/functions.R")
analysis <- readRDS("results/g4_analysis.rds")
motif <- readRDS("results/motif_analysis.rds")

profile <- read.delim("data/processed/tss_profile.tsv")
profile$condition <- sub("^Inh5$", "N-inh5", profile$condition)
profile$condition <- sub("^NT$", "DMSO", profile$condition)
profile_levels <- c("DMSO", "RHPS4", "N-inh5", "Combo")
profile_colors <- c(DMSO = condition_colors[["NT"]], condition_colors[c("RHPS4", "N-inh5", "Combo")])
profile$condition <- factor(profile$condition, levels = profile_levels)

p_a1 <- ggplot(profile, aes(position_bp, signal, colour = condition)) +
  annotate("rect", xmin = -500, xmax = 500, ymin = -Inf, ymax = Inf, fill = "grey90") +
  geom_vline(xintercept = 0, linetype = 3, colour = "grey55", linewidth = 0.25) +
  geom_line(linewidth = 0.55) +
  scale_colour_manual(values = profile_colors, name = NULL) +
  scale_x_continuous(breaks = c(-10000, -5000, 0, 5000, 10000), labels = c("−10 kb", "−5 kb", "0", "5 kb", "10 kb")) +
  labs(x = NULL, y = "Mean G4 signal (RPM)") +
  theme_figure(8) +
  theme(legend.position = c(0.86, 0.76), legend.background = element_blank())

p_a2 <- ggplot(profile, aes(position_bp, signal, colour = condition)) +
  geom_vline(xintercept = 0, linetype = 3, colour = "grey65", linewidth = 0.2) +
  geom_line(linewidth = 0.42) +
  facet_wrap(~condition, nrow = 1) +
  scale_colour_manual(values = profile_colors, guide = "none") +
  scale_x_continuous(breaks = c(-10000, 0, 10000), labels = c("−10", "0", "10")) +
  labs(x = NULL, y = NULL) +
  theme_figure(6) +
  theme(strip.background = element_blank(), strip.text = element_text(size = 7), axis.text.y = element_blank(), axis.ticks.y = element_blank())

tss_matrix <- read_deeptools_matrix("data/processed/tss_condition_matrix.tsv.gz", 4L)
heatmap_rows <- lapply(seq_along(condition_levels), function(i) {
  j <- ((i - 1L) * tss_matrix$bins + 1L):(i * tss_matrix$bins)
  x <- tss_matrix$values[, j, drop = FALSE]
  ranking <- order(rowMeans(x[, 96:105, drop = FALSE]), decreasing = TRUE)
  z <- x[ranking[unique(round(seq(1, length(ranking), length.out = 900)))], , drop = FALSE]
  data.frame(
    condition = profile_levels[i],
    rank = rep(seq_len(nrow(z)), times = ncol(z)),
    position_bp = rep(seq(-9950, 9950, length.out = ncol(z)), each = nrow(z)),
    signal = as.vector(z)
  )
})
heatmap_data <- do.call(rbind, heatmap_rows)
heatmap_data$condition <- factor(heatmap_data$condition, levels = profile_levels)

p_a3 <- ggplot(heatmap_data, aes(position_bp, rank, fill = signal)) +
  geom_raster() +
  facet_wrap(~condition, nrow = 1) +
  scale_y_reverse() +
  scale_x_continuous(breaks = c(-10000, 0, 10000), labels = c("−10", "0", "10")) +
  scale_fill_gradient(low = "white", high = "#b2182b", limits = c(0, 1.2), oob = squish, name = "G4 signal\n(RPM)") +
  labs(x = "Distance to TSS (kb)", y = "Genes (ranked)") +
  theme_minimal(base_size = 6) +
  theme(panel.grid = element_blank(), strip.text = element_blank(), axis.text.y = element_blank(), legend.position = "right")

panel_a <- p_a1 / p_a2 / p_a3 + plot_layout(heights = c(1.05, 0.42, 1.4))
save_ggplot(panel_a, "Figure7A", 5.0, 6.2)
rm(tss_matrix, heatmap_data, heatmap_rows)
gc()

pairwise_significant <- lapply(analysis$pairwise, function(x) x[!is.na(x$padj) & x$padj < 0.05 & x$CHR %in% standard_chromosomes, ])
consensus <- data.frame(chr = analysis$coordinates$CHR, start = analysis$coordinates$START, end = analysis$coordinates$END)
consensus <- consensus[consensus$chr %in% standard_chromosomes, ]

draw_circos <- function() {
  circos.clear()
  circos.par(start.degree = 90, gap.after = c(rep(1.2, 23), 18), cell.padding = c(0, 0, 0, 0), track.margin = c(0.003, 0.003))
  circos.initializeWithIdeogram(species = "hg38", chromosome.index = standard_chromosomes, plotType = c("ideogram", "labels"), track.height = 0.045, labels.cex = 0.5)
  circos.genomicDensity(consensus, col = "grey45", window.size = 1e6, track.height = 0.055)
  for (name in c("RHPS4", "N-inh5", "Combo")) {
    x <- pairwise_significant[[name]]
    gain <- data.frame(chr = x$CHR[x$log2FoldChange > 0], start = x$START[x$log2FoldChange > 0], end = x$END[x$log2FoldChange > 0], value = pmin(x$log2FoldChange[x$log2FoldChange > 0], 25))
    loss <- data.frame(chr = x$CHR[x$log2FoldChange < 0], start = x$START[x$log2FoldChange < 0], end = x$END[x$log2FoldChange < 0], value = pmax(x$log2FoldChange[x$log2FoldChange < 0], -25))
    circos.genomicTrackPlotRegion(list(gain, loss), ylim = c(-25, 25), track.height = 0.075, bg.border = NA, panel.fun = function(region, value, ...) {
      i <- getI(...)
      circos.genomicPoints(region, value, pch = 16, cex = 0.45, col = if (i == 1L) direction_colors[["Gain"]] else direction_colors[["Loss"]])
      circos.lines(CELL_META$cell.xlim, c(0, 0), col = "grey80", lty = 3, lwd = 0.3)
    })
  }
  x <- motif$interaction
  syn <- data.frame(chr = x$CHR[x$direction == "Synergistic"], start = x$START[x$direction == "Synergistic"], end = x$END[x$direction == "Synergistic"], value = 1)
  ant <- data.frame(chr = x$CHR[x$direction == "Antagonistic"], start = x$START[x$direction == "Antagonistic"], end = x$END[x$direction == "Antagonistic"], value = -1)
  circos.genomicTrackPlotRegion(list(syn, ant), ylim = c(-2, 2), track.height = 0.06, bg.border = NA, panel.fun = function(region, value, ...) {
    i <- getI(...)
    circos.genomicPoints(region, value, pch = 17, cex = 0.72, col = if (i == 1L) direction_colors[["Synergistic"]] else direction_colors[["Antagonistic"]])
  })
  legend("center", legend = c("Gain", "Loss", "Pairwise peak", "Interaction peak"), col = c(direction_colors[["Gain"]], direction_colors[["Loss"]], "black", "black"), pch = c(16, 16, 16, 17), bty = "n", cex = 0.65, title = "G4 change vs DMSO")
  circos.clear()
}

pdf("figures/Figure7B.pdf", width = 7, height = 7)
par(mar = c(0, 0, 0, 0))
draw_circos()
dev.off()
png("figures/Figure7B.png", width = 2100, height = 2100, res = 300)
par(mar = c(0, 0, 0, 0))
draw_circos()
dev.off()

combo <- analysis$pairwise$Combo
annotations <- read.csv(gzfile("data/processed/g4_peak_annotations.csv.gz"), check.names = FALSE)
annotation_key <- peak_key(annotations$seqnames, annotations$start, annotations$end)
k <- peak_key(combo$CHR, combo$START, combo$END)
combo$gene <- annotations$SYMBOL[match(k, annotation_key)]
combo$location <- annotations$loc_cat[match(k, annotation_key)]
combo <- combo[!is.na(combo$gene) & combo$gene != "", ]
combo <- combo[order(is.na(combo$padj), combo$padj, -abs(combo$log2FoldChange)), ]
combo <- combo[!duplicated(combo$gene), ]
rna <- read.csv(gzfile("data/processed/rna_contrasts.csv.gz"))
scatter <- merge(combo[, c("gene", "log2FoldChange", "padj", "location")], rna[, c("gene", "rna_combo_fc")], by = "gene")
names(scatter)[names(scatter) == "log2FoldChange"] <- "g4_fc"
names(scatter)[names(scatter) == "rna_combo_fc"] <- "rna_fc"
scatter <- scatter[is.finite(scatter$g4_fc) & is.finite(scatter$rna_fc), ]

lower_right <- scatter[scatter$g4_fc > 0 & scatter$rna_fc < 0, ]
lower_right$score <- rank(lower_right$g4_fc) + rank(-lower_right$rna_fc)
upper_left <- scatter[scatter$g4_fc < 0 & scatter$rna_fc > 0, ]
upper_left$score <- rank(-upper_left$g4_fc) + rank(upper_left$rna_fc)
labels <- unique(c("BAP1", "ATP13A3-DT", "RBM15-AS1", "SLC29A1", "NINJ2", "RAD51", "RPA2", head(lower_right$gene[order(-lower_right$score)], 20), head(upper_left$gene[order(-upper_left$score)], 20)))
scatter$label <- ifelse(scatter$gene %in% labels, scatter$gene, NA_character_)
scatter$location <- factor(scatter$location, levels = names(location_colors))

panel_c <- ggplot(scatter, aes(g4_fc, rna_fc)) +
  geom_hline(yintercept = 0, linetype = 3, colour = "grey65", linewidth = 0.25) +
  geom_vline(xintercept = 0, linetype = 3, colour = "grey65", linewidth = 0.25) +
  geom_point(colour = "grey78", size = 0.4, alpha = 0.35) +
  geom_point(data = scatter[!is.na(scatter$label), ], aes(fill = location), shape = 21, size = 2.1, stroke = 0.2) +
  geom_text_repel(data = scatter[!is.na(scatter$label), ], aes(label = label), size = 2.2, seed = 42, min.segment.length = 0, segment.size = 0.18, max.overlaps = Inf) +
  scale_fill_manual(values = location_colors, drop = FALSE, name = "Peak location") +
  labs(x = expression(G4~Combo~vs~DMSO~log[2]~FC), y = expression(RNA~Combo~vs~DMSO~log[2]~FC)) +
  theme_figure(8) +
  theme(legend.position = c(0.86, 0.82), legend.background = element_rect(fill = alpha("white", 0.8), colour = NA))
save_ggplot(panel_c, "Figure7C", 5.6, 4.3)
write_tsv(scatter, "results/figure7C_source_data.tsv")

bars <- motif$motif_summary
bars$set <- factor(bars$set, levels = c("Promoter interaction peaks", "All interaction peaks", "All G4 peaks", "Random regions"))
bar_colors <- c(`Random regions` = "#aab0b4", `All G4 peaks` = "#849794", `All interaction peaks` = "#8e44ad", `Promoter interaction peaks` = "#ef8a17")
panel_d1 <- ggplot(bars, aes(percent, set, fill = set)) +
  geom_col(width = 0.64) +
  geom_text(aes(label = paste0("n=", n)), hjust = -0.12, size = 2.4, fontface = "italic") +
  scale_fill_manual(values = bar_colors, guide = "none") +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 80, 20), expand = expansion(mult = c(0, 0))) +
  labs(x = "% peaks with predicted\nG4 motif (pqsfinder ≥ 40)", y = NULL) +
  theme_figure(7)

contingency <- as.data.frame(motif$promoter_table)
names(contingency) <- c("location", "direction", "count")
contingency$location <- factor(contingency$location, levels = c("Promoter", "Non-promoter"))
panel_d2 <- ggplot(contingency, aes(direction, location, fill = count)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = count), size = 4) +
  scale_fill_gradient(low = "#f8e9e4", high = "#e45745", guide = "none") +
  labs(x = NULL, y = NULL, title = sprintf("Fisher p = %.3f", motif$promoter_test$p.value)) +
  theme_minimal(base_size = 7) +
  theme(panel.grid = element_blank(), plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 35, hjust = 1))

panel_d <- panel_d1 + panel_d2 + plot_layout(widths = c(1.8, 1))
save_ggplot(panel_d, "Figure7D", 6.1, 2.7)

peak_signal <- read.delim(gzfile("data/processed/interaction_peak_signal.tsv.gz"), check.names = FALSE)
signal_columns <- c("RHPS4.1", "RHPS4.2", "RHPS4.3", "Inh5.1", "Inh5.2", "Inh5.3", "Combo.1", "Combo.2", "Combo.3")
heat <- data.frame(
  gene = rep(peak_signal$label, each = length(signal_columns)),
  sample = rep(signal_columns, times = nrow(peak_signal)),
  value = as.vector(t(as.matrix(peak_signal[, signal_columns])))
)
heat$gene <- factor(heat$gene, levels = peak_signal$label)
heat$sample <- factor(heat$sample, levels = rev(signal_columns))
levels(heat$sample) <- gsub("^(RHPS4|Inh5|Combo)\\.([123])$", "\\1 (\\2)", levels(heat$sample))

annotation_data <- rbind(
  data.frame(gene = peak_signal$label, strip = "Location", value = peak_signal$loc_cat),
  data.frame(gene = peak_signal$label, strip = "Direction", value = peak_signal$direction)
)
annotation_data$gene <- factor(annotation_data$gene, levels = peak_signal$label)
annotation_data$strip <- factor(annotation_data$strip, levels = c("Direction", "Location"))
annotation_palette <- c(location_colors, direction_colors[c("Synergistic", "Antagonistic")])

p_e_anno <- ggplot(annotation_data, aes(gene, strip, fill = value)) +
  geom_tile(colour = "white", linewidth = 0.15) +
  scale_fill_manual(values = annotation_palette, name = NULL) +
  theme_void(base_size = 6) +
  theme(axis.text.y = element_text(), axis.title.y = element_blank(), legend.position = "right", legend.text = element_text(size = 5), legend.key.height = grid::unit(2.4, "mm"))

p_e_heat <- ggplot(heat, aes(gene, sample, fill = value)) +
  geom_tile(colour = "white", linewidth = 0.12) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", limits = c(-5, 5), oob = squish, name = "Log2FC vs DMSO") +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 6) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5), legend.position = "right")

bap1 <- read.delim(gzfile("data/processed/bap1_locus_signal.tsv.gz"), check.names = FALSE)
bap1_samples <- setdiff(names(bap1), "position")
bap1_long <- data.frame(position = rep(bap1$position, times = length(bap1_samples)), sample = rep(bap1_samples, each = nrow(bap1)), signal = unlist(bap1[, bap1_samples], use.names = FALSE))
bap1_long$condition <- sub("rep[123]$", "", bap1_long$sample)
bap1_long$condition <- sub("^Inh5$", "N-inh5", bap1_long$condition)
bap1_long$condition <- sub("^Combo$", "Combo", bap1_long$condition)
bap1_long$condition <- factor(bap1_long$condition, levels = condition_levels)
bap1_mean <- aggregate(signal ~ position + condition, bap1_long, mean)

p_e_bap1 <- ggplot() +
  annotate("rect", xmin = 52406264, xmax = 52407442, ymin = -Inf, ymax = Inf, fill = "#fee8a6", alpha = 0.45) +
  geom_line(data = bap1_long, aes(position, signal, group = sample, colour = condition), linewidth = 0.18, alpha = 0.45) +
  geom_line(data = bap1_mean, aes(position, signal, colour = condition), linewidth = 0.55) +
  geom_vline(xintercept = c(52406264, 52407442), linetype = 3, colour = "#c9a227", linewidth = 0.2) +
  facet_grid(condition ~ .) +
  scale_colour_manual(values = condition_colors, guide = "none") +
  scale_x_continuous(breaks = c(52406264, 52407442), labels = c("52,406,264\n(peak start)", "52,407,442\n(peak end)")) +
  labs(x = NULL, y = "G4 signal (RPM)", title = "BAP1 promoter interaction peak region") +
  theme_figure(6) +
  theme(strip.background = element_blank(), strip.text.y = element_text(angle = 0), plot.title = element_text(hjust = 0.5, size = 7))

panel_e_left <- p_e_anno / p_e_heat + plot_layout(heights = c(0.18, 1))
panel_e <- panel_e_left | p_e_bap1
panel_e <- panel_e + plot_layout(widths = c(2.6, 1))
save_ggplot(panel_e, "Figure7E", 12.0, 4.4)

gene_groups <- list(
  `Promoter synergy` = c("BAP1", "CTDNEP1", "PEG3", "PARTICL"),
  `Direct G4/RNA loci` = c("ATP13A3-DT", "RBM15-AS1", "SLC29A1", "NINJ2"),
  `p21 response` = c("CDKN1A", "MIR34AHG"),
  `DREAM/cell-cycle` = c("MYBL2", "LIN9", "CDC20", "BIRC5", "CDK1", "TOP2A", "MKI67"),
  Senescence = c("LMNB1", "LMNB2")
)
genes <- unlist(gene_groups, use.names = FALSE)
groups <- rep(names(gene_groups), lengths(gene_groups))
rna_sub <- rna[match(genes, rna$gene), ]
rna_matrix <- cbind(NT = 0, RHPS4 = rna_sub$rna_rhps4_fc, `N-inh5` = rna_sub$rna_inh5_fc, Combo = rna_sub$rna_combo_fc)
rownames(rna_matrix) <- genes
z_matrix <- t(scale(t(rna_matrix)))
z_matrix[z_matrix > 2] <- 2
z_matrix[z_matrix < -2] <- -2
expression_data <- data.frame(
  gene = rep(genes, each = 4),
  group = rep(groups, each = 4),
  condition = rep(colnames(z_matrix), times = length(genes)),
  z = as.vector(t(z_matrix))
)
expression_data$gene <- factor(expression_data$gene, levels = rev(genes))
expression_data$group <- factor(expression_data$group, levels = names(gene_groups))
expression_data$condition <- factor(expression_data$condition, levels = colnames(z_matrix))

panel_f <- ggplot(expression_data, aes(condition, gene, fill = z)) +
  geom_tile(colour = "white", linewidth = 0.25) +
  facet_grid(group ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", limits = c(-2, 2), name = "z-score") +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 7) +
  theme(panel.grid = element_blank(), strip.placement = "outside", strip.background = element_blank(), strip.text.y.left = element_text(angle = 0, hjust = 1), axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "top")
save_ggplot(panel_f, "Figure7F", 4.0, 5.5)
write_tsv(expression_data, "results/figure7F_source_data.tsv")

panel_raster <- function(path, label) {
  image <- png::readPNG(path)
  ggplot() +
    annotation_raster(image, xmin = 0, xmax = 1, ymin = 0, ymax = 1) +
    annotate("text", x = 0, y = 1, label = label, hjust = -0.15, vjust = 1.1, fontface = "bold", size = 5) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
    theme_void()
}

composite <- (panel_raster("figures/Figure7A.png", "A") + panel_raster("figures/Figure7B.png", "B") + (panel_raster("figures/Figure7C.png", "C") / panel_raster("figures/Figure7D.png", "D")) + plot_layout(widths = c(1.0, 1.35, 1.05))) /
  (panel_raster("figures/Figure7E.png", "E") + panel_raster("figures/Figure7F.png", "F") + plot_layout(widths = c(2.8, 1))) +
  plot_layout(heights = c(1.05, 0.72))
save_ggplot(composite, "Figure7_reproduced", 14.8, 9.7)
