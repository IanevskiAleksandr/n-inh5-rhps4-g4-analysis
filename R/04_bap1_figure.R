suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

source("R/functions.R")

shrna <- read.delim("data/external/bap1_shrna.tsv")
shrna$screen <- factor(shrna$screen, levels = shrna$screen)

panel_a <- ggplot(shrna, aes(screen, fraction)) +
  geom_hline(yintercept = 0.5, linetype = 2, colour = "grey45", linewidth = 0.35) +
  geom_col(width = 0.62, fill = "#7455a6") +
  geom_text(aes(label = count_label), vjust = -0.45, size = 3.2, fontface = "bold") +
  scale_y_continuous(limits = c(0, 1.12), breaks = seq(0, 1, 0.25), labels = scales::percent_format()) +
  labs(x = NULL, y = "BAP1-targeting shRNAs depleted", title = "Pyridostatin sensitization") +
  theme_figure(9) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

crispr <- read.delim("data/external/bap1_crispr_31screens.tsv")
crispr$screen_rank <- factor(crispr$screen_rank, levels = crispr$screen_rank)
crispr$label <- ifelse(crispr$highlight == "Other treatments", NA_character_, crispr$screen)

panel_b <- ggplot(crispr, aes(screen_rank, normZ)) +
  geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.3) +
  geom_segment(aes(xend = screen_rank, y = 0, yend = normZ), colour = "grey80", linewidth = 0.35) +
  geom_point(aes(colour = highlight), size = 2.4) +
  geom_text(data = crispr[!is.na(crispr$label), ], aes(label = label), vjust = 1.4, hjust = 0, size = 3) +
  scale_colour_manual(values = c(`Other treatments` = "#a9afb4", `Phen-DC3` = "#d44e41", Pyridostatin = "#7455a6"), guide = "none") +
  scale_x_discrete(breaks = levels(crispr$screen_rank)[c(1, 5, 10, 15, 20, 25, 31)]) +
  labs(x = "Rank across 31 treatments", y = "BAP1 knockout response (normZ)", title = "Independent CRISPR–drug screens") +
  theme_figure(9)

supplementary_figure <- panel_a + panel_b + plot_annotation(tag_levels = "a") + plot_layout(widths = c(1, 1.45))
save_ggplot(supplementary_figure, "SupplementaryFigure5_BAP1", 10.2, 4.3)
