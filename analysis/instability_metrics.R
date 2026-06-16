#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
})

# ----------------------------
# 0) Config
# ----------------------------
READS_PATH <- "data/allSamples-readsdata.txt"
META_PATH  <- "data/sample-metadata.tsv"

OUT_DIR_FIG <- "figures"
OUT_DIR_RES <- "results"

PM_MIN <- 55
PM_MAX <- 199

dir.create(OUT_DIR_FIG, showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_DIR_RES, showWarnings = FALSE, recursive = TRUE)

# ----------------------------
# 1) Utilities
# ----------------------------
assert_has_cols <- function(df, cols, df_name = "data.frame") {
  missing <- setdiff(cols, names(df))
  if (length(missing) > 0) {
    stop(df_name, " is missing required columns: ", paste(missing, collapse = ", "))
  }
}

density_mode <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 3 || length(unique(x)) < 2) return(NA_real_)
  d <- density(x)
  d$x[which.max(d$y)]
}

calc_instability_metrics <- function(df, pm_min = 55, pm_max = 200) {
  df %>%
    mutate(allele_length = as.numeric(allele_length)) %>%
    filter(
      !is.na(allele_length),
      allele_length >= pm_min,
      allele_length <= pm_max
    ) %>%
    group_by(bam_id) %>%
    summarise(
      n_pm = n(),
      modal_pm = density_mode(allele_length),
      mean_length = mean(allele_length, na.rm = TRUE),
      median_length = median(allele_length, na.rm = TRUE),
      instability_abs_mean = mean(abs(allele_length - modal_pm), na.rm = TRUE),
      instability_right_mean = mean(pmax(allele_length - modal_pm, 0), na.rm = TRUE),
      instability_right_median = median(pmax(allele_length - modal_pm, 0), na.rm = TRUE),
      pct_right = mean(allele_length > modal_pm, na.rm = TRUE),
      .groups = "drop"
    )
}

save_model_summary <- function(model, path) {
  capture.output(summary(model), file = path)
}

save_cor_test <- function(x, y, path, method = "spearman") {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3 || length(unique(x[ok])) < 2 || length(unique(y[ok])) < 2) {
    writeLines("Correlation test not run: insufficient finite/nonconstant observations.", path)
  } else {
    capture.output(
      suppressWarnings(cor.test(x[ok], y[ok], method = method)),
      file = path
    )
  }
}

# ----------------------------
# 2) Load data
# ----------------------------
if (!file.exists(READS_PATH)) stop("READS_PATH not found: ", READS_PATH)
if (!file.exists(META_PATH))  stop("META_PATH not found: ", META_PATH)

dat <- read.delim(
  READS_PATH,
  header = TRUE,
  sep = "\t",
  fill = TRUE,
  stringsAsFactors = FALSE
)

meta <- read_tsv(META_PATH, show_col_types = FALSE)

if ("Batch/Run" %in% names(meta)) {
  meta <- meta %>% rename(batch_run = `Batch/Run`)
}

if ("Sex" %in% names(meta)) {
  meta <- meta %>%
    mutate(Sex = tools::toTitleCase(tolower(Sex)))
}

assert_has_cols(dat, c("bam_id", "allele_length"), "reads data")
assert_has_cols(meta, c("bam_id"), "metadata")

# ----------------------------
# 3) Calculate instability metrics
# ----------------------------
instab_all <- dat %>%
  semi_join(meta, by = "bam_id") %>%
  calc_instability_metrics(pm_min = PM_MIN, pm_max = PM_MAX) %>%
  left_join(
    meta %>%
      select(any_of(c("bam_id", "Sex", "Diagnosis", "Region", "batch_run"))),
    by = "bam_id"
  )

write_csv(instab_all, file.path(OUT_DIR_RES, "instability_metrics_by_sample.csv"))

# ----------------------------
# 4) Statistical tests
# ----------------------------
model_abs <- lm(instability_abs_mean ~ modal_pm, data = instab_all)
model_right <- lm(instability_right_mean ~ modal_pm, data = instab_all)

save_model_summary(
  model_abs,
  file.path(OUT_DIR_RES, "lm_abs_instability_vs_modal_repeat.txt")
)

save_model_summary(
  model_right,
  file.path(OUT_DIR_RES, "lm_right_instability_vs_modal_repeat.txt")
)

save_cor_test(
  instab_all$modal_pm,
  instab_all$instability_abs_mean,
  file.path(OUT_DIR_RES, "spearman_abs_instability_vs_modal_repeat.txt")
)

save_cor_test(
  instab_all$modal_pm,
  instab_all$instability_right_mean,
  file.path(OUT_DIR_RES, "spearman_right_instability_vs_modal_repeat.txt")
)

# ----------------------------
# 5) Plots
# ----------------------------
p_abs <- ggplot(
  instab_all,
  aes(x = modal_pm, y = instability_abs_mean, color = Sex)
) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  theme_classic() +
  labs(
    x = "Modal premutation repeat length",
    y = "Instability index: mean |Δ from modal|",
    color = "Sex",
    title = "Bidirectional instability vs premutation repeat length"
  )

ggsave(
  file.path(OUT_DIR_FIG, "instability_abs_vs_modal_repeat.pdf"),
  p_abs,
  width = 4.5,
  height = 3.5,
  useDingbats = FALSE
)

p_right <- ggplot(
  instab_all,
  aes(x = modal_pm, y = instability_right_mean, color = Sex)
) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  theme_classic() +
  labs(
    x = "Modal premutation repeat length",
    y = "Right-tail instability: mean(max(Δ, 0))",
    color = "Sex",
    title = "Right-tail instability vs premutation repeat length"
  )

ggsave(
  file.path(OUT_DIR_FIG, "instability_right_vs_modal_repeat.pdf"),
  p_right,
  width = 4.5,
  height = 3.5,
  useDingbats = FALSE
)

message("Done. Instability metrics saved to: ", OUT_DIR_RES)
