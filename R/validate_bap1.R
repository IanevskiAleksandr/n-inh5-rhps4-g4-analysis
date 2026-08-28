#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(readxl))

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Run with Rscript R/validate_bap1.R", call. = FALSE)
root <- dirname(dirname(normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)))
raw_dir <- Sys.getenv("BAP1_SOURCE_DIR", file.path(root, "data", "raw_external", "bap1"))
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)

sources <- data.frame(
  file = c(
    "elife-46793-supp1-v1.xlsx",
    "elife-46793-supp2-v1.xlsx",
    "Zscores_31screens.csv",
    "FDRNeg_31screens.csv"
  ),
  url = c(
    "https://cdn.elifesciences.org/articles/46793/elife-46793-supp1-v1.xlsx",
    "https://cdn.elifesciences.org/articles/46793/elife-46793-supp2-v1.xlsx",
    "https://data.mendeley.com/public-files/datasets/gfcn2wmrpf/files/f27bf615-7747-46d4-8f78-39cb422895a0/file_downloaded",
    "https://data.mendeley.com/public-files/datasets/gfcn2wmrpf/files/454913cb-0223-4af6-9ea9-f591b7a96312/file_downloaded"
  ),
  sha256 = c(
    "80ac5b5836ee9ce955d09ed8f3c9a9c32f232f186af18a382891458a5c1fb183",
    "98e458a97835652d1d497818c8adb7800ff8376db6e30fdb61ea520c7234f486",
    "7f36e57221240ab3b65d54baabc648e684110c8e35ee0ac1b458317c39e48d33",
    "9a220861ae704635dae29a2bedc87acfc9c351b937cb3a06791a3f408e488c9e"
  ),
  stringsAsFactors = FALSE
)

fetch <- function(file, url, sha256) {
  path <- file.path(raw_dir, file)
  if (!file.exists(path)) {
    tmp <- tempfile(tmpdir = raw_dir)
    utils::download.file(url, tmp, mode = "wb", quiet = TRUE)
    observed <- unname(tools::sha256sum(tmp))
    if (!identical(observed, sha256)) stop("Checksum failed for ", file, call. = FALSE)
    if (!file.rename(tmp, path)) stop("Could not install ", file, call. = FALSE)
  }
  observed <- unname(tools::sha256sum(path))
  if (!identical(observed, sha256)) stop("Checksum failed for ", path, call. = FALSE)
  path
}

paths <- mapply(fetch, sources$file, sources$url, sources$sha256, USE.NAMES = FALSE)
names(paths) <- sources$file

assert <- function(ok, label) if (!isTRUE(ok)) stop("Validation failed: ", label, call. = FALSE)
near <- function(x, y, tolerance = 1e-12) length(x) == length(y) && all(abs(x - y) <= tolerance)

gene_row <- function(file, sheet) {
  x <- as.data.frame(read_excel(file, sheet = sheet), check.names = FALSE)
  out <- x[x$Gene == "BAP1", , drop = FALSE]
  assert(nrow(out) == 1L, paste(sheet, "contains one BAP1 row"))
  out
}

gw <- gene_row(paths[["elife-46793-supp1-v1.xlsx"]], "Fig 3B GW PDS_genes")
a375 <- gene_row(paths[["elife-46793-supp1-v1.xlsx"]], "Fig 6B A375 Focused_PDS")
ht1080 <- gene_row(paths[["elife-46793-supp2-v1.xlsx"]], "Fig7Sup1D_PDS_genes")
ht_shrna <- as.data.frame(
  read_excel(paths[["elife-46793-supp2-v1.xlsx"]], sheet = "Fig7Sup1C_PDS_shRNAs"),
  check.names = FALSE
)
ht_shrna <- ht_shrna[ht_shrna$gene == "BAP1", , drop = FALSE]

shrna <- data.frame(
  screen = c("A375 genome-wide", "A375 focused", "HT1080 focused"),
  depleted = as.integer(c(gw$`Depleted shRNAs`, a375$`Depleted shRNAs`, ht1080$`Depleted shRNAs`)),
  detected = as.integer(c(gw$`shRNAs in t0`, a375$`shRNAs in t0`, ht1080$`shRNAs in t0`)),
  stringsAsFactors = FALSE
)
shrna$fraction <- shrna$depleted / shrna$detected
shrna$count_label <- sprintf("%d/%d", shrna$depleted, shrna$detected)
shrna$percent_label <- sprintf("%.0f%%", 100 * shrna$fraction)

assert(identical(shrna$depleted, c(4L, 7L, 4L)), "Zyner depleted shRNA counts")
assert(identical(shrna$detected, c(8L, 7L, 7L)), "Zyner detected shRNA counts")
assert(nrow(ht_shrna) == 4L && all(ht_shrna$FDR < 0.05), "four significant HT1080 BAP1 shRNAs")

supplied_shrna <- read.delim(
  file.path(root, "data", "external", "bap1_shrna.tsv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
assert(identical(supplied_shrna$screen, shrna$screen), "supplied shRNA screen order")
assert(identical(as.integer(supplied_shrna$depleted), shrna$depleted), "supplied depleted counts")
assert(identical(as.integer(supplied_shrna$detected), shrna$detected), "supplied detected counts")
assert(near(supplied_shrna$fraction, shrna$fraction), "supplied depleted fractions")
assert(identical(supplied_shrna$count_label, shrna$count_label), "supplied count labels")
assert(identical(supplied_shrna$percent_label, shrna$percent_label), "supplied percentage labels")

z <- read.csv(paths[["Zscores_31screens.csv"]], check.names = FALSE, stringsAsFactors = FALSE)
fdr <- read.csv(paths[["FDRNeg_31screens.csv"]], check.names = FALSE, stringsAsFactors = FALSE)
z_cols <- names(z)[-1L]
fdr_cols <- names(fdr)[-1L]
fdr_key <- sub("^HU\\.1$", "HU", fdr_cols)
assert(length(z_cols) == 31L && identical(z_cols, fdr_key), "31 aligned Olivieri screens")

z_row <- z[z$Gene == "BAP1", , drop = FALSE]
fdr_row <- fdr[fdr$gene == "BAP1", , drop = FALSE]
assert(nrow(z_row) == 1L && nrow(fdr_row) == 1L, "one Olivieri BAP1 row")

gene_rank <- function(column) {
  values <- suppressWarnings(as.numeric(z[[column]]))
  target <- which(z$Gene == "BAP1")
  c(rank = rank(values, ties.method = "min", na.last = "keep")[target], total = sum(!is.na(values)))
}

crispr <- data.frame(
  screen = ifelse(z_cols == "PhenDC3", "Phen-DC3", z_cols),
  source_column = z_cols,
  normZ = as.numeric(z_row[1L, z_cols]),
  FDR = as.numeric(fdr_row[1L, fdr_cols]),
  stringsAsFactors = FALSE
)
crispr <- crispr[order(crispr$normZ, crispr$screen), , drop = FALSE]
crispr$screen_rank <- seq_len(nrow(crispr))
crispr$highlight <- ifelse(
  crispr$screen == "Phen-DC3", "Phen-DC3",
  ifelse(crispr$screen == "Pyridostatin", "Pyridostatin", "Other treatments")
)
ranks <- t(vapply(crispr$source_column, gene_rank, numeric(2L)))
crispr$gene_rank <- as.integer(ranks[, "rank"])
crispr$n_genes <- as.integer(ranks[, "total"])

supplied_crispr <- read.delim(
  file.path(root, "data", "external", "bap1_crispr_31screens.tsv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
text_cols <- c("screen", "source_column", "highlight")
number_cols <- c("normZ", "FDR", "screen_rank", "gene_rank", "n_genes")
assert(all(vapply(text_cols, function(column) identical(supplied_crispr[[column]], crispr[[column]]), logical(1L))), "supplied CRISPR labels")
for (column in number_cols) assert(near(supplied_crispr[[column]], crispr[[column]]), paste("supplied", column))

phen <- crispr[crispr$screen == "Phen-DC3", , drop = FALSE]
pds <- crispr[crispr$screen == "Pyridostatin", , drop = FALSE]
assert(near(phen$normZ, -5.81) && phen$screen_rank == 1L, "Phen-DC3 normZ and screen rank")
assert(phen$gene_rank == 7L && phen$n_genes == 17298L && near(phen$FDR, 7.87e-6), "Phen-DC3 gene rank and DrugZ FDR")
assert(near(pds$normZ, -2.56) && pds$screen_rank == 3L, "pyridostatin normZ and screen rank")
assert(pds$gene_rank == 163L && pds$n_genes == 17374L && near(pds$FDR, 0.552), "pyridostatin gene rank and DrugZ FDR")

report <- data.frame(
  check = c(
    "A375 genome-wide PDS shRNAs",
    "A375 focused PDS shRNAs",
    "HT1080 focused PDS shRNAs",
    "HT1080 significant shRNA rows",
    "Olivieri screens",
    "Phen-DC3 BAP1",
    "Pyridostatin BAP1",
    "Compact source tables"
  ),
  observed = c(
    "4/8", "7/7", "4/7", "4",
    "31",
    "normZ=-5.81; screen=1/31; gene=7/17298; DrugZ FDRneg=7.87e-06",
    "normZ=-2.56; screen=3/31; gene=163/17374; DrugZ FDRneg=0.552",
    "exact match"
  ),
  status = "PASS",
  stringsAsFactors = FALSE
)
dir.create(file.path(root, "results"), recursive = TRUE, showWarnings = FALSE)
write.table(
  report,
  file.path(root, "results", "bap1_validation.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
cat("BAP1 external validation: PASS\n")
