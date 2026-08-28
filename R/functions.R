options(stringsAsFactors = FALSE)

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

standard_chromosomes <- paste0("chr", c(1:22, "X", "Y"))
condition_levels <- c("NT", "RHPS4", "N-inh5", "Combo")
condition_colors <- c(NT = "#30343b", RHPS4 = "#c4473a", `N-inh5` = "#2b7bba", Combo = "#7e57a6")
direction_colors <- c(Gain = "#cf4c3c", Loss = "#2b83ba", Synergistic = "#e76f51", Antagonistic = "#4b9cd3")
location_colors <- c(Promoter = "#d95f59", UTR = "#e6a21a", Exon = "#2ca25f", Intron = "#4c91c6", Intergenic = "#7b61a8")

peak_key <- function(chr, start, end) paste(chr, start, end, sep = ":")

geometric_mean <- function(x) exp(mean(log(x)))

read_deeptools_matrix <- function(path, blocks) {
  z <- readLines(gzfile(path))[-1]
  x <- do.call(rbind, strsplit(z, "\t", fixed = TRUE))
  coords <- x[, 1:6, drop = FALSE]
  values <- matrix(as.numeric(x[, -(1:6), drop = FALSE]), nrow = nrow(x))
  stopifnot(ncol(values) %% blocks == 0)
  list(coords = coords, values = values, bins = ncol(values) / blocks)
}

write_tsv <- function(x, path) write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE)

save_ggplot <- function(plot, stem, width, height) {
  ggplot2::ggsave(file.path("figures", paste0(stem, ".pdf")), plot, width = width, height = height, device = grDevices::cairo_pdf)
  ggplot2::ggsave(file.path("figures", paste0(stem, ".png")), plot, width = width, height = height, dpi = 300)
}

theme_figure <- function(base_size = 8) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      axis.title = ggplot2::element_text(colour = "black"),
      axis.text = ggplot2::element_text(colour = "black"),
      plot.title = ggplot2::element_text(face = "bold"),
      legend.title = ggplot2::element_text(face = "bold")
    )
}
