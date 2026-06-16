#!/usr/bin/env Rscript

# QC for PureTarget BAMs
# Goal: per-sample summaries of read quality (rq/RQ -> predicted QV),
# read length, MAPQ, and optional tandem repeat locus tag (TR).
#
# Notes:
# - PacBio rq/RQ tags are used when available.
# - Some older BAMs may not contain rq/RQ metadata; QV-derived metrics are
#   reported as NA for those samples.
# - BAM files are expected to be indexed.

suppressPackageStartupMessages({
  library(Rsamtools)
  library(dplyr)
  library(purrr)
  library(stringr)
  library(ggplot2)
  library(readr)
  library(tibble)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

# ----------------------------
# 0) Config
# ----------------------------
BAM_DIR <- "data/spanning_bams"

OUT_DIR_FIG <- "figures"
OUT_DIR_RES <- "results"

dir.create(OUT_DIR_FIG, showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_DIR_RES, showWarnings = FALSE, recursive = TRUE)

# ----------------------------
# 1) Helper functions
# ----------------------------
sample_name <- function(path) {
  x <- basename(path)
  x <- str_replace(x, "\\.bam$", "")
  x <- str_replace(x, "\\.sorted$", "")
  x
}

qv_from_rq <- function(rq) {
  case_when(
    is.na(rq) ~ NA_real_,
    rq >= 1 ~ 60,
    rq <= 0 ~ NA_real_,
    TRUE ~ -10 * log10(1 - rq)
  )
}

tag_to_char <- function(tag, n) {
  if (is.null(tag)) return(rep(NA_character_, n))
  if (is.list(tag)) {
    return(vapply(tag, function(z) {
      if (is.null(z) || length(z) == 0) NA_character_ else as.character(z[[1]])
    }, character(1)))
  }
  as.character(tag)
}

tag_to_int <- function(tag, n) {
  if (is.null(tag)) return(rep(NA_integer_, n))
  if (is.list(tag)) {
    return(vapply(tag, function(z) {
      if (is.null(z) || length(z) == 0) NA_integer_ else as.integer(z[[1]])
    }, integer(1)))
  }
  as.integer(tag)
}

tag_to_num <- function(tag, n) {
  if (is.null(tag)) return(rep(NA_real_, n))
  if (is.list(tag)) {
    return(vapply(tag, function(z) {
      if (is.null(z) || length(z) == 0) NA_real_ else as.numeric(z[[1]])
    }, numeric(1)))
  }
  as.numeric(tag)
}

extract_read_qc <- function(bam_path) {
  param <- ScanBamParam(
    what = c("qname", "flag", "rname", "pos", "mapq", "cigar", "seq", "qual"),
    tag = c("rq", "RQ", "RG", "HP", "TR")
  )
  
  x <- scanBam(bam_path, param = param)[[1]]
  n <- length(x$qname)
  
  rq_tag <- x$tag$rq %||% x$tag$RQ
  
  tibble(
    bam = basename(bam_path),
    sample = sample_name(bam_path),
    qname = x$qname,
    flag = x$flag,
    chr = as.character(x$rname),
    pos = x$pos,
    mapq = x$mapq,
    cigar = x$cigar,
    read_len = nchar(as.character(x$seq)),
    rq = tag_to_num(rq_tag, n),
    RG = tag_to_char(x$tag$RG, n),
    HP = tag_to_int(x$tag$HP, n),
    TR = tag_to_char(x$tag$TR, n)
  ) %>%
    mutate(
      qv_phred = qv_from_rq(rq),
      qv_cap = pmin(qv_phred, 60),
      has_rq = !is.na(rq)
    )
}

# ----------------------------
# 2) Load BAMs and check indexes
# ----------------------------
bams <- list.files(BAM_DIR, pattern = "\\.bam$", full.names = TRUE)
stopifnot(length(bams) > 0)

bai1 <- paste0(bams, ".bai")
bai2 <- str_replace(bams, "\\.bam$", ".bai")
has_index <- file.exists(bai1) | file.exists(bai2)

if (any(!has_index)) {
  stop(
    "Missing BAM index for:\n",
    paste(basename(bams[!has_index]), collapse = "\n"),
    "\nFix with: samtools index <bam>"
  )
}

message("All BAMs indexed.")

# ----------------------------
# 3) Extract read-level QC
# ----------------------------
read_qc <- map_dfr(bams, extract_read_qc)

write_csv(
  read_qc,
  file.path(OUT_DIR_RES, "read_level_quality_metrics.csv")
)

# ----------------------------
# 4) Summarize QC by sample
# ----------------------------
read_quality_summary <- read_qc %>%
  group_by(sample) %>%
  summarise(
    n_reads = n(),
    frac_has_rq = mean(has_rq, na.rm = TRUE),
    mean_qv_phred = mean(qv_phred, na.rm = TRUE),
    median_qv_phred = median(qv_phred, na.rm = TRUE),
    frac_qv60 = mean(qv_phred >= 60, na.rm = TRUE),
    median_read_len = median(read_len, na.rm = TRUE),
    p90_read_len = as.numeric(quantile(read_len, 0.90, na.rm = TRUE)),
    median_mapq = median(mapq, na.rm = TRUE),
    frac_mapq0 = mean(mapq == 0, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    across(
      c(mean_qv_phred, median_qv_phred, frac_qv60),
      ~ ifelse(frac_has_rq == 0, NA_real_, .x)
    )
  ) %>%
  arrange(desc(n_reads))

write_csv(
  read_quality_summary,
  file.path(OUT_DIR_RES, "read_quality_summary_by_sample.csv")
)

# ----------------------------
# 5) Optional locus-level read counts
# ----------------------------
locus_read_counts <- read_qc %>%
  count(sample, TR, name = "n_reads") %>%
  arrange(sample, desc(n_reads))

write_csv(
  locus_read_counts,
  file.path(OUT_DIR_RES, "read_counts_by_sample_and_locus.csv")
)

# ----------------------------
# 6) QC plots
# ----------------------------
p_len_qv <- read_qc %>%
  filter(!is.na(qv_cap)) %>%
  ggplot(aes(x = read_len, y = qv_cap)) +
  geom_hex(bins = 60) +
  scale_x_log10() +
  theme_bw() +
  labs(
    x = "Read length (bp, log10)",
    y = "Predicted QV (capped at 60)",
    title = "Read length vs predicted QV"
  )

ggsave(
  file.path(OUT_DIR_FIG, "read_length_vs_predicted_qv.pdf"),
  p_len_qv,
  width = 6,
  height = 4,
  units = "in",
  dpi = 300
)

p_locus_counts <- locus_read_counts %>%
  filter(!is.na(TR)) %>%
  ggplot(aes(x = TR, y = sample, fill = n_reads)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c(trans = "log10", name = "Spanning reads") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 8)
  ) +
  labs(
    x = "Repeat locus",
    y = "Sample",
    title = "Spanning read coverage per sample and locus"
  )

ggsave(
  file.path(OUT_DIR_FIG, "spanning_read_counts_by_sample_and_locus.pdf"),
  p_locus_counts,
  width = 8,
  height = 7,
  units = "in",
  dpi = 300
)

message("Done. Read-quality QC outputs saved to: ", OUT_DIR_RES, " and ", OUT_DIR_FIG)
