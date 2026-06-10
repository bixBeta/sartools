#!/usr/bin/env Rscript
# split_counts.R
# Splits a counts matrix into per-sample files.
# Each output file: two columns (feature-ID + sample counts), named <sample>.rawCounts
#
# Usage:
#   Rscript split_counts.R --input counts_matrix.tsv [--outdir ./split_counts] [--header]
#
# Flags:
#   --input    Path to counts matrix (tab or comma delimited; auto-detected)
#   --outdir   Output directory (default: current directory)
#   --header   If set, output files include a header line (default: no header)

suppressPackageStartupMessages(library(tidyverse))

# ── CLI args ──────────────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)

parse_arg <- function(flag, default = NULL) {
  idx <- which(args == flag)
  if (length(idx) && length(args) >= idx + 1) args[idx + 1] else default
}

has_flag <- function(flag) flag %in% args

input_file  <- parse_arg("--input")
out_dir     <- parse_arg("--outdir", default = ".")
keep_header <- has_flag("--header")

if (is.null(input_file)) {
  stop("Usage: Rscript split_counts.R --input <counts_matrix> [--outdir <output_dir>] [--header]")
}

if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}

message(sprintf("Header in output: %s", if (keep_header) "yes" else "no"))

# ── Auto-detect delimiter (tab or comma) ─────────────────────────────────────
detect_delim <- function(path, n = 5) {
  lines  <- readLines(path, n = n, warn = FALSE)
  tabs   <- sum(stringr::str_count(lines, "\t"))
  commas <- sum(stringr::str_count(lines, ","))
  if (tabs >= commas) "\t" else ","
}

delim <- detect_delim(input_file)
message(sprintf("Detected delimiter: %s", if (delim == "\t") "tab" else "comma"))

# ── Read counts matrix ────────────────────────────────────────────────────────
# First column is treated as the feature/gene ID regardless of its name.
counts <- read_delim(input_file, delim = delim, col_types = cols(.default = "c"),
                     show_col_types = FALSE, trim_ws = TRUE)

gene_col    <- colnames(counts)[1]   # whatever the first column is called
sample_cols <- colnames(counts)[-1]

message(sprintf("Read %d genes x %d samples from: %s", nrow(counts), length(sample_cols), input_file))
message(sprintf("Output directory: %s", out_dir))

# ── Create output directory if needed ────────────────────────────────────────
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
  message("Created output directory: ", out_dir)
}

# ── Split and write ───────────────────────────────────────────────────────────
walk(sample_cols, function(sample) {
  out_path <- file.path(out_dir, paste0(sample, ".rawCounts"))

  counts |>
    select(all_of(gene_col), all_of(sample)) |>
    write_tsv(out_path, col_names = keep_header)

  message("  Written: ", basename(out_path))
})

message(sprintf("\nDone. %d files written to: %s", length(sample_cols), out_dir))
