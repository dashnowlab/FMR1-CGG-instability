#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(ggridges)
})

# ----------------------------
# 0) Config (edit these)
# ----------------------------
READS_PATH <- "data/allSamples-readsdata.txt"
META_PATH  <- "data/sample-metadata.tsv"

OUT_DIR_FIG <- "figures"
OUT_DIR_RES <- "results"

PM_MIN <- 55
PM_MAX <- 199      
FM_MIN <- 200
FM_MAX <- 6000     

METH_THR <- 0.6
COR_METHOD <- "spearman"
BANDWIDTH <- 1.5

# ----------------------------
# 1) Decide cohort- or sample-level analysis
# ----------------------------
COHORT_FILTER <- function(meta) {
  meta %>%
    filter(Diagnosis == "FXTAS") %>%
    filter(Region %in% c("BA10","CBL","PVWM"))
}

# If you want a within-sample report, set ONE_SAMPLE_ID.
# If NULL, runs cohort-level only.
ONE_SAMPLE_ID <- NULL
# ONE_SAMPLE_ID <- "bc2011"

# Save filtered per-read data used for analysis?
WRITE_FILTERED_DATA <- TRUE

dir.create(OUT_DIR_FIG, showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_DIR_RES, showWarnings = FALSE, recursive = TRUE)

# ----------------------------
# 2) Utilities
# ----------------------------
assert_has_cols <- function(df, cols, df_name = "data.frame") {
  missing <- setdiff(cols, names(df))
  if (length(missing) > 0) {
    stop(df_name, " is missing required columns: ", paste(missing, collapse = ", "))
  }
}

safe_cor <- function(x, y, method = "spearman") {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]
  if (length(x) < 3 || length(unique(x)) < 2 || length(unique(y)) < 2) return(NA_real_)
  suppressWarnings(cor(x, y, method = method))
}

prep_repeat_meth <- function(df, repeat_min, repeat_max, max_inclusive = FALSE) {
  df %>%
    mutate(
      allele_length = as.numeric(allele_length),
      denom = methylated_bases + umethylated_bases,
      meth_frac = methylated_bases / denom
    ) %>%
    filter(
      !is.na(allele_length),
      is.finite(meth_frac),
      denom > 0
    ) %>%
    { if (max_inclusive) filter(., allele_length >= repeat_min, allele_length <= repeat_max)
      else               filter(., allele_length >= repeat_min, allele_length <  repeat_max)
    }
}

calc_bin_pct <- function(df, thr = 0.6) {
  df %>%
    group_by(bam_id) %>%
    summarise(
      pct_high = mean(meth_frac >= thr, na.rm = TRUE),
      pct_low  = mean(meth_frac <  thr, na.rm = TRUE),
      n_reads  = n(),
      .groups = "drop"
    ) %>%
    mutate(lab = sprintf("H=%0.1f%%", 100 * pct_high))
}

calc_coupling_rho <- function(df, method = "spearman") {
  df %>%
    group_by(bam_id) %>%
    summarise(
      rho = safe_cor(allele_length, meth_frac, method = method),
      n_reads = n(),
      .groups = "drop"
    )
}

order_ids_by_median_len <- function(df) {
  df %>%
    group_by(bam_id) %>%
    summarise(med_len = median(allele_length, na.rm = TRUE), .groups = "drop") %>%
    arrange(med_len) %>%
    pull(bam_id)
}

plot_ridges_meth_bins <- function(df, thr = 0.6, xlim = c(55, 200),
                                  title = "", bandwidth = 1.5) {
  id_levels <- order_ids_by_median_len(df)
  
  dfp <- df %>%
    mutate(
      bam_id = factor(bam_id, levels = id_levels),
      meth_bin = ifelse(meth_frac >= thr, "High", "Low"),
      meth_bin = factor(meth_bin, levels = c("Low", "High"))
    )
  
  ggplot(dfp, aes(x = allele_length, y = bam_id)) +
    geom_density_ridges(
      data = subset(dfp, meth_bin == "Low"),
      alpha = 0.55, scale = 1.0, bandwidth = bandwidth
    ) +
    geom_density_ridges(
      data = subset(dfp, meth_bin == "High"),
      alpha = 0.35, scale = 1.0, bandwidth = bandwidth, size = 1.2
    ) +
    coord_cartesian(xlim = xlim) +
    theme_void() +
    theme(
      axis.text.y  = element_text(size = 8, color = "black"),
      axis.text.x  = element_text(size = 8, color = "black"),
      plot.title   = element_text(size = 11, face = "bold")
    ) +
    labs(x = "Repeat length (reads)", y = NULL, title = title)
}

save_outputs <- function(prefix, rho_tbl, bin_tbl, plot_obj) {
  write_csv(rho_tbl, file.path(OUT_DIR_RES, paste0(prefix, "_rho.csv")))
  write_csv(bin_tbl, file.path(OUT_DIR_RES, paste0(prefix, "_binpct.csv")))
  ggsave(
    filename = file.path(OUT_DIR_FIG, paste0(prefix, "_ridges.pdf")),
    plot = plot_obj, width = 6.5, height = 3.2, useDingbats = FALSE
  )
}

# ----------------------------
# 3) Load data
# ----------------------------
if (!file.exists(READS_PATH)) stop("READS_PATH not found: ", READS_PATH)
if (!file.exists(META_PATH))  stop("META_PATH not found: ", META_PATH)

dat <- read.delim(
  READS_PATH, header = TRUE, sep = "\t", fill = TRUE,
  stringsAsFactors = FALSE
)
meta <- read_tsv(META_PATH, show_col_types = FALSE)

meta <- meta %>%
  mutate(Sex = tools::toTitleCase(tolower(Sex)))

# normalize the Batch/Run column if present
if ("Batch/Run" %in% names(meta)) {
  meta <- meta %>% rename(batch_run = `Batch/Run`)
}

assert_has_cols(dat,  c("bam_id", "allele_length", "methylated_bases", "umethylated_bases"), "reads data")
assert_has_cols(meta, c("bam_id"), "metadata")

meta_subset <- COHORT_FILTER(meta)
if (nrow(meta_subset) == 0) stop("COHORT_FILTER returned 0 rows. Check metadata columns/values.")

dat_cohort <- dat %>%
  semi_join(meta_subset, by = "bam_id") %>%
  filter(!is.na(allele_length))

if (!is.null(ONE_SAMPLE_ID)) {
  dat_cohort <- dat_cohort %>% filter(bam_id == ONE_SAMPLE_ID)
  if (nrow(dat_cohort) == 0) stop("ONE_SAMPLE_ID not found after filtering: ", ONE_SAMPLE_ID)
}


# ----------------------------
# 4) Run PM workflow
# ----------------------------
pm_df <- prep_repeat_meth(dat_cohort, PM_MIN, PM_MAX, max_inclusive = FALSE)

if (WRITE_FILTERED_DATA) {
  write_csv(pm_df, file.path(OUT_DIR_RES, "PM_filtered_per_read.csv"))
}

pm_rho <- calc_coupling_rho(pm_df, COR_METHOD)
pm_bin <- calc_bin_pct(pm_df, METH_THR)

pm_plot <- plot_ridges_meth_bins(
  pm_df, thr = METH_THR, xlim = c(PM_MIN, PM_MAX),
  title = ifelse(is.null(ONE_SAMPLE_ID), "Cohort: Premutation", paste0(ONE_SAMPLE_ID, ": Premutation")),
  bandwidth = BANDWIDTH
)

save_outputs(prefix = "PM", rho_tbl = pm_rho, bin_tbl = pm_bin, plot_obj = pm_plot)

# Cohort-level test (if cohort mode)
if (is.null(ONE_SAMPLE_ID)) {
  pm_wilcox <- wilcox.test(pm_rho$rho, mu = 0)
  capture.output(pm_wilcox, file = file.path(OUT_DIR_RES, "PM_wilcox_rho_vs0.txt"))
}

# ----------------------------
# 5) Run FM workflow (optional)
# ----------------------------
fm_df <- prep_repeat_meth(dat_cohort, FM_MIN, FM_MAX, max_inclusive = FALSE)

if (nrow(fm_df) > 0) {
  if (WRITE_FILTERED_DATA) {
    write_csv(fm_df, file.path(OUT_DIR_RES, "FM_filtered_per_read.csv"))
  }

  fm_rho <- calc_coupling_rho(fm_df, COR_METHOD)
  fm_bin <- calc_bin_pct(fm_df, METH_THR)
  
  fm_plot <- plot_ridges_meth_bins(
    fm_df, thr = METH_THR, xlim = c(FM_MIN, min(1300, FM_MAX)),
    title = ifelse(is.null(ONE_SAMPLE_ID), "Cohort: Full mutation", paste0(ONE_SAMPLE_ID, ": Full mutation")),
    bandwidth = BANDWIDTH
  )
  
  save_outputs(prefix = "FM", rho_tbl = fm_rho, bin_tbl = fm_bin, plot_obj = fm_plot)
  
  if (is.null(ONE_SAMPLE_ID)) {
    fm_wilcox <- wilcox.test(fm_rho$rho, mu = 0)
    capture.output(fm_wilcox, file = file.path(OUT_DIR_RES, "FM_wilcox_rho_vs0.txt"))
  }
} else {
  message("No FM-range reads found in this dataset selection.")
}

message("Done. Outputs saved to: ", OUT_DIR_RES, " and ", OUT_DIR_FIG)
