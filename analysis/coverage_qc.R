#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Rsamtools)
  library(GenomicRanges)
  library(IRanges)
  library(dplyr)
  library(purrr)
  library(stringr)
  library(tibble)
  library(ggplot2)
  library(ggbeeswarm)
  library(readr)
})

# ----------------------------
# 0) Config
# ----------------------------
BAM_DIR <- "data/spanning_bams"
BED_PATH <- "data/PureTarget_repeat_expansion_panel.bed"

OUT_DIR_FIG <- "figures"
OUT_DIR_RES <- "results"

dir.create(OUT_DIR_FIG, showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_DIR_RES, showWarnings = FALSE, recursive = TRUE)

# ----------------------------
# 1) Helper functions
# ----------------------------
sample_name <- function(path) {
  str_replace(basename(path), "\\.bam$", "")
}

read_bed_loci <- function(bed_path) {
  x <- readLines(bed_path)
  x <- x[nzchar(x)]
  x <- x[!grepl("^#", x)]
  
  spl <- strsplit(x, "\\s+")
  spl <- spl[lengths(spl) >= 3]
  
  bed_df <- tibble(
    seqnames = vapply(spl, `[[`, character(1), 1),
    start0   = as.integer(vapply(spl, `[[`, character(1), 2)),
    end0     = as.integer(vapply(spl, `[[`, character(1), 3)),
    info     = ifelse(
      lengths(spl) >= 4,
      vapply(spl, function(z) paste(z[4:length(z)], collapse = " "), character(1)),
      NA_character_
    )
  ) %>%
    mutate(
      start = start0 + 1L,
      end = end0
    )
  
  bed_df %>%
    mutate(
      locus = str_match(info, "ID=([^;]+)")[, 2],
      locus = ifelse(is.na(locus), paste0(seqnames, ":", start, "-", end), locus)
    ) %>%
    transmute(locus, seqnames, start, end)
}

assign_group <- function(sample) {
  bc <- str_extract(sample, "bc\\d+")
  
  fxpm_bcs <- c(
    "bc2033", "bc2034", "bc2037", "bc2038", "bc2039", "bc2040", "bc2041", "bc2042",
    "bc2006", "bc2007", "bc2008", "bc2009", "bc2010", "bc2011", "bc2012", "bc2013",
    "bc2014", "bc2018", "bc2019", "bc2020", "bc2023", "bc2024", "bc2043"
  )
  
  fxs_bcs <- c(
    "bc2015", "bc2022", "bc2035", "bc2036", "bc2044",
    "bc2002", "bc2003", "bc2004", "bc2005"
  )
  
  con_bcs <- c(
    "bc2001", "bc2016", "bc2017", "bc2021",
    "bc2045", "bc2046", "bc2047", "bc2048"
  )
  
  case_when(
    bc %in% fxs_bcs ~ "FXS",
    bc %in% fxpm_bcs ~ "FXPM",
    bc %in% con_bcs ~ "CON",
    TRUE ~ "UNKNOWN"
  )
}

cov_one <- function(bam, locus_row, total_reads, min_mapq = 0) {
  sn <- sample_name(bam)
  
  gr <- GRanges(
    seqnames = locus_row$seqnames,
    ranges = IRanges(start = locus_row$start, end = locus_row$end)
  )
  
  pu <- Rsamtools::pileup(
    file = bam,
    scanBamParam = ScanBamParam(which = gr),
    pileupParam = PileupParam(
      distinguish_nucleotides = FALSE,
      distinguish_strands = FALSE,
      min_mapq = min_mapq,
      include_deletions = TRUE,
      include_insertions = TRUE
    )
  )
  
  pu_depth <- as_tibble(pu) %>%
    group_by(seqnames, pos) %>%
    summarise(depth = sum(count), .groups = "drop")
  
  full_pos <- tibble(
    seqnames = as.character(locus_row$seqnames),
    pos = locus_row$start:locus_row$end
  )
  
  full_pos %>%
    left_join(pu_depth, by = c("seqnames", "pos")) %>%
    mutate(
      depth = ifelse(is.na(depth), 0L, as.integer(depth)),
      sample = sn,
      locus = locus_row$locus,
      rel_pos = pos - locus_row$start,
      depth_norm = depth / total_reads
    ) %>%
    select(sample, locus, seqnames, pos, rel_pos, depth, depth_norm)
}

# ----------------------------
# 2) Load BAMs and BED
# ----------------------------
bams <- list.files(BAM_DIR, pattern = "\\.bam$", full.names = TRUE)
stopifnot(length(bams) > 0)

bai1 <- paste0(bams, ".bai")
bai2 <- str_replace(bams, "\\.bam$", ".bai")
has_index <- file.exists(bai1) | file.exists(bai2)

if (any(!has_index)) {
  stop("Missing BAM index for:\n", paste(basename(bams[!has_index]), collapse = "\n"))
}

message("All BAMs indexed.")

loci <- read_bed_loci(BED_PATH)

meta <- tibble(sample = sample_name(bams)) %>%
  mutate(
    bc = str_extract(sample, "bc\\d+"),
    group = assign_group(sample)
  )

if (any(meta$group == "UNKNOWN")) {
  warning("Some samples were assigned UNKNOWN group.")
}

# ----------------------------
# 3) Calculate normalized coverage
# ----------------------------
get_total_reads <- function(bam) countBam(bam)$records

total_spanning <- setNames(
  vapply(bams, get_total_reads, numeric(1)),
  sample_name(bams)
)

cov_df <- map_dfr(bams, function(bam) {
  sn <- sample_name(bam)
  tr <- total_spanning[[sn]]
  
  pmap_dfr(loci, function(locus, seqnames, start, end) {
    cov_one(
      bam = bam,
      locus_row = tibble(
        locus = locus,
        seqnames = seqnames,
        start = start,
        end = end
      ),
      total_reads = tr
    )
  })
})

# ----------------------------
# 4) Summarize coverage
# ----------------------------
locus_sample <- cov_df %>%
  group_by(locus, sample) %>%
  summarise(
    mean_norm_cov = mean(depth_norm, na.rm = TRUE),
    median_norm_cov = median(depth_norm, na.rm = TRUE),
    frac_zero = mean(depth == 0, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(locus = factor(locus, levels = loci$locus))

locus_sample2 <- locus_sample %>%
  left_join(meta %>% select(sample, group), by = "sample") %>%
  mutate(
    group = factor(group, levels = c("CON", "FXPM", "FXS")),
    fmr1_flag = if_else(locus == "FMR1", "FMR1", "Other")
  )

locus_group_summary <- locus_sample2 %>%
  group_by(locus, group) %>%
  summarise(
    mean_cov = mean(mean_norm_cov, na.rm = TRUE),
    sd_cov = sd(mean_norm_cov, na.rm = TRUE),
    n = sum(!is.na(mean_norm_cov)),
    sem_cov = sd_cov / sqrt(n),
    .groups = "drop"
  )

# ----------------------------
# 5) Save tables
# ----------------------------
write_csv(cov_df, file.path(OUT_DIR_RES, "per_base_normalized_coverage.csv"))
write_csv(locus_sample2, file.path(OUT_DIR_RES, "locus_sample_coverage.csv"))
write_csv(locus_group_summary, file.path(OUT_DIR_RES, "locus_group_coverage_summary.csv"))

# ----------------------------
# 6) Generate QC plots
# ----------------------------
p_beeswarm <- ggplot(
  locus_sample2,
  aes(x = locus, y = mean_norm_cov, color = group)
) +
  geom_quasirandom(width = 0.25, alpha = 0.85, size = 1.8) +
  stat_summary(
    fun = mean,
    geom = "crossbar",
    width = 0.6,
    color = "black",
    linewidth = 0.6
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank()
  ) +
  labs(
    x = "Locus",
    y = "Mean normalized spanning coverage",
    color = "Group"
  )

ggsave(
  file.path(OUT_DIR_FIG, "TRGT_spanning_coverage_beeswarm_mean.pdf"),
  p_beeswarm,
  width = 8,
  height = 4,
  units = "in",
  dpi = 300
)

p_fmr1 <- ggplot(
  locus_sample2,
  aes(x = locus, y = mean_norm_cov, color = fmr1_flag)
) +
  geom_quasirandom(width = 0.25, alpha = 0.85, size = 1.8) +
  stat_summary(
    fun = mean,
    geom = "crossbar",
    width = 0.6,
    color = "black",
    linewidth = 0.6
  ) +
  facet_wrap(~ group, scales = "fixed") +
  scale_color_manual(values = c("FMR1" = "#D55E00", "Other" = "grey70")) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank()
  ) +
  labs(
    x = "Locus",
    y = "Mean normalized spanning coverage",
    color = NULL
  )

ggsave(
  file.path(OUT_DIR_FIG, "TRGT_spanning_coverage_FMR1_highlight.pdf"),
  p_fmr1,
  width = 10,
  height = 4,
  units = "in",
  dpi = 300
)

p_bar <- ggplot(
  locus_group_summary,
  aes(x = group, y = mean_cov, fill = group)
) +
  geom_col(width = 0.7, color = "black") +
  geom_errorbar(
    aes(ymin = mean_cov - sem_cov, ymax = mean_cov + sem_cov),
    width = 0.2,
    linewidth = 0.6
  ) +
  facet_wrap(~ locus) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.major.x = element_blank(),
    legend.position = "none"
  ) +
  labs(
    x = "Condition",
    y = "Mean normalized spanning coverage"
  )

ggsave(
  file.path(OUT_DIR_FIG, "TRGT_spanning_coverage_bar_by_condition_facet_locus.pdf"),
  p_bar,
  width = 12,
  height = 7,
  units = "in",
  dpi = 300
)

message("Done. Coverage QC outputs saved to: ", OUT_DIR_RES, " and ", OUT_DIR_FIG)
