### Combined Script: Boolean Filtering + Normalization (01_ + 02_)

# ─────────────────────────────────────────────────────
# Load Libraries
# ─────────────────────────────────────────────────────
library(DESeq2)
library(tidyverse)
library(eulerr)
library(GenomicRanges)
library(VennDiagram)
library(pheatmap)
library(akima)
library(ggrepel)
library(UpSetR)
library(RColorBrewer)

# ─────────────────────────────────────────────────────
# Utility Functions
# ─────────────────────────────────────────────────────
w <- function(a){ print(head(a)); print(nrow(a)) }
create_dir <- function(dir_path) { if (!dir.exists(dir_path)) dir.create(dir_path, recursive = TRUE) }
base_theme <- theme_classic() + theme(
  text = element_text(color = "black"),
  plot.title = element_text(hjust = 0.5, size = rel(1.5)),
  axis.text.x = element_text(size = rel(1.25)),
  axis.text.y = element_text(size = rel(1.25), angle = 90, hjust = 0.5),
  axis.title.y = element_text(size = rel(1.5))
)

# ─────────────────────────────────────────────────────
# Set Paths
# ─────────────────────────────────────────────────────
main_wd   <- "/Users/brandiatteberry/Desktop/Bioinformatics/ATACseq_Nextflow/"; create_dir(main_wd)
macs_wd   <- file.path(main_wd, "macs");  create_dir(macs_wd)
deseq_wd  <- file.path(main_wd, "R");     create_dir(deseq_wd)
scale_wd  <- file.path(main_wd, "scale"); create_dir(scale_wd)
output_wd <- file.path(deseq_wd, toupper(format(Sys.Date(), "%d%b%y"))); create_dir(output_wd)
homer_wd  <- file.path(deseq_wd, "Homer"); create_dir(homer_wd)

current_date  <- toupper(format(Sys.Date(), "%d%b%y"))
boolean_file  <- "/Users/brandiatteberry/Desktop/Bioinformatics/RBC_Lysis_Bioinformatics_training/consensus_peaks.mRp.clN.boolean.txt"
sample_map    <- "/Users/brandiatteberry/Desktop/Bioinformatics/RBC_Lysis_Bioinformatics_training/deseq_input_all.csv"
feature_counts_file <- "/Users/brandiatteberry/Desktop/Bioinformatics/RBC_Lysis_Bioinformatics_training/consensus_peaks.mRp.clN.featureCounts.txt"

# ─────────────────────────────────────────────────────
# Boolean Filtering + UpSetR Visualization (01_)
# ─────────────────────────────────────────────────────
sample_renaming <- read.csv(sample_map, check.names = FALSE)
consensus_samples_a <- sample_renaming %>%
  filter(Comp_A == "Untreated") %>%
  distinct(analysis_name_simple, .keep_all = TRUE)
consensus_samples_b <- consensus_samples_a[, c("sample", "analysis_name_simple")]
consensus_samples   <- consensus_samples_b$analysis_name_simple
consensus_peak_overlap_samples <- 3

main_df <- read.table(boolean_file, header = TRUE, sep = "\t", as.is = TRUE)

bool_data <- main_df %>% select(ends_with(".bool"))
# Strip ANY stem like ".mRp.clN.bool" or ".mLb.clN.bool" → keep sample name only
colnames(bool_data) <- gsub("\\.[^.]+\\.clN\\.bool$", "", colnames(bool_data))
# Map to analysis_name_simple where possible
bool_data <- bool_data %>%
  rename_with(~ consensus_samples_b$analysis_name_simple[match(.x, consensus_samples_b$sample)])
bool_data <- if (length(consensus_samples) == 0) bool_data else bool_data %>% select(any_of(consensus_samples))

filtered_data <- bool_data %>%
  mutate(across(everything(), as.integer)) %>%
  filter(rowSums(.) >= 1)

# UpSetR plot
png(file.path(output_wd, paste0("upsetR_plot_", current_date, ".png")), width = 1800, height = 1000, res = 120)
upset(data = filtered_data,
      sets = names(filtered_data),
      order.by = "freq",
      main.bar.color = "grey",
      matrix.color = "#7E2945",
      text.scale = 1.25,
      keep.order = TRUE,
      mb.ratio = c(0.75, 0.25))
dev.off()

# Overlap bar chart
upset_overlaps <- tibble()
for (i in 2:ncol(bool_data)) {
  n <- bool_data %>%
    mutate(across(everything(), as.integer)) %>%
    filter(rowSums(.) >= i) %>% nrow()
  upset_overlaps <- bind_rows(upset_overlaps, tibble(Overlaps = i, Peaks = n))
}
upset_overlap_plot <- ggplot(upset_overlaps, aes(x = Overlaps, y = Peaks)) +
  geom_bar(stat = "identity", fill = "grey") +
  geom_text(aes(label = Peaks), vjust = -0.5, color = "black", size = 5) +
  labs(title = "Peaks vs. Overlaps", x = "Number of Overlaps", y = "Number of Peaks") +
  base_theme
ggsave(file.path(output_wd, paste0("upset_overlap_per_sample_", current_date, ".png")),
       plot = upset_overlap_plot, width = 6, height = 8, dpi = 500)

# BED output (all intervals)
consensus_peaks <- data.frame(main_df[, 1:4], Height = 1, Strand = "+")
write.table(consensus_peaks,
            file.path(output_wd, paste0("consensus_peaks_", consensus_peak_overlap_samples, "samples_", current_date, ".bed")),
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

# ─────────────────────────────────────────────────────
# featureCounts Processing for Normalization Input
# ─────────────────────────────────────────────────────
fc_raw <- read.table(feature_counts_file, header = TRUE, comment.char = "#", check.names = FALSE)

file_first_six   <- file.path(scale_wd, "consensus_peaks_first_six_columns.txt")
file_bam_columns <- file.path(scale_wd, "consensus_peaks_sorted_bam_columns.txt")

write.table(fc_raw[, 1:6],            file_first_six,   sep = "\t", row.names = FALSE, quote = FALSE)
write.table(fc_raw[, 7:ncol(fc_raw)], file_bam_columns, sep = "\t", row.names = FALSE, quote = FALSE)
message("✅ featureCounts extraction complete: Files written to `scale` directory")

# ─────────────────────────────────────────────────────
# Generate Scale Factor Files (robust & copy-safe)
# ─────────────────────────────────────────────────────
counts <- read.table(file_bam_columns, header = TRUE, sep = "\t", check.names = FALSE)

# Clean BAM column names to "<sample>" by stripping ".[stem].clN.sorted.bam"
bam_samples_clean <- sub("\\.[^.]+\\.clN\\.sorted\\.bam$", "", colnames(counts))

# Sanity: counts has one column per sample
stopifnot(length(bam_samples_clean) == ncol(counts))

sample_totals <- colSums(counts)
max_total <- max(sample_totals)

for (j in seq_along(bam_samples_clean)) {
  sample_name <- bam_samples_clean[j]
  sf <- round(sample_totals[[j]] / max_total, 5)
  fname <- file.path(scale_wd, sprintf("%s.auto.clN.scale_factor.txt", sample_name))
  # use writeLines to avoid odd console prints / factors
  writeLines(format(sf, scientific = FALSE), con = fname)
}

message("✅ Scale factor files generated to: ", scale_wd)

# ─────────────────────────────────────────────────────
# Normalization & Overlap Filtering (02_)
# ─────────────────────────────────────────────────────
bool_timepoints <- c("_T10_CI","_T30_CI",
                     "_T10_DMSO","_T30_DMSO","_T60_DMSO","_T90_DMSO","_T120_DMSO",
                     "_T30_PMA","_T60_PMA","_T90_PMA","_T120_PMA")
num_overlaps <- 2

# Guard: (re)create split files if missing (e.g., new session)
if (!file.exists(file_first_six) || !file.exists(file_bam_columns)) {
  fc_raw <- read.table(feature_counts_file, header = TRUE, comment.char = "#", check.names = FALSE)
  write.table(fc_raw[, 1:6],            file_first_six,   sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(fc_raw[, 7:ncol(fc_raw)], file_bam_columns, sep = "\t", row.names = FALSE, quote = FALSE)
  message("✅ Regenerated split featureCounts files in `scale/`")
}

# Read & harmonize
data_a <- read.table(file_first_six, header = TRUE, sep = "\t", check.names = FALSE) %>%
  rename(chr = Chr, start = Start, end = End)
data_b <- read.table(file_bam_columns, header = TRUE, sep = "\t", check.names = FALSE)

# Clean BAM column names to "<sample>"
colnames(data_b) <- gsub("\\.[^.]+\\.clN\\.sorted\\.bam$", "", colnames(data_b))
rownames(data_b) <- data_a$Geneid

# ── Prefer our new .auto scale-factor files; else fallback to legacy files ──
scale_files_all <- list.files(scale_wd, pattern = "\\.scale_factor\\.txt$", full.names = TRUE)
base_all <- basename(scale_files_all)

auto_idx <- grepl("\\.auto\\.clN\\.scale_factor\\.txt$", base_all)
if (any(auto_idx)) {
  scale_files <- scale_files_all[auto_idx]
  base_sf <- base_all[auto_idx]
  # <sample>.auto.clN.scale_factor.txt  → extract <sample>
  sf_samples <- sub("\\.auto\\.clN\\.scale_factor\\.txt$", "", base_sf)
} else {
  scale_files <- scale_files_all
  base_sf <- base_all
  # Legacy: may look like <sample>.<stem>.clN.scale_factor.txt OR
  #         <sample>.<stem>.clN.sorted.bam.clN.scale_factor.txt (seen in your dir)
  # Start by removing the trailing ".clN.scale_factor.txt"
  temp <- sub("\\.clN\\.scale_factor\\.txt$", "", base_sf)
  # If it still ends with ".<stem>.clN.sorted.bam", strip that too
  temp <- sub("\\.[^.]+\\.clN\\.sorted\\.bam$", "", temp)
  # If it still ends with ".<stem>", strip that
  temp <- sub("\\.[^.]+$", "", temp)
  sf_samples <- temp
}

# Read numeric scale values
sf_values <- vapply(scale_files, function(f) {
  x <- scan(f, what = numeric(), quiet = TRUE)
  if (length(x) == 0) NA_real_ else x[1]
}, numeric(1))

scaling_factors_df <- tibble(name = sf_samples, scaling_factor = sf_values) %>%
  filter(!is.na(scaling_factor))

# Build a one-to-one mapping by prefix: choose the first file whose name starts with the cleaned sample
overlap <- intersect(scaling_factors_df$name, colnames(data_b))
if (length(overlap) == 0) {
  # Try prefix match if exact names don't intersect
  matched <- tibble(sample = character(), scaling_factor = numeric())
  for (s in colnames(data_b)) {
    idx <- which(startsWith(scaling_factors_df$name, paste0(s)))
    if (length(idx) == 0) next
    # Prefer the shortest name match (most specific)
    idx <- idx[order(nchar(scaling_factors_df$name[idx]))][1]
    matched <- bind_rows(matched, tibble(sample = s, scaling_factor = scaling_factors_df$scaling_factor[idx]))
  }
  if (nrow(matched) == 0) {
    stop("No overlapping sample names between BAM columns and scale-factor files. ",
         "Cleaned BAM names (e.g., ", paste(head(colnames(data_b)), collapse=", "), 
         ") must prefix-match scale-factor files (e.g., ", paste(head(scaling_factors_df$name), collapse=", "), ").")
  }
  overlap <- matched$sample
  # Rebuild scaling_factors_df to only matched samples
  scaling_factors_df <- matched %>% rename(name = sample)
} else {
  # filter to exact overlap only
  scaling_factors_df <- scaling_factors_df %>% filter(name %in% overlap)
}

# Apply scaling
new_data_b <- purrr::map_dfc(scaling_factors_df$name, function(sample) {
  sc <- scaling_factors_df %>% filter(name == sample) %>% pull(scaling_factor)
  round(100 * data_b[[sample]] * sc)
})
colnames(new_data_b) <- scaling_factors_df$name

# Carry interval IDs explicitly and add coordinates
data_b_export <- tibble(Interval = data_a$Geneid) %>% bind_cols(new_data_b)
data_b_df <- left_join(data_b_export, data_a, by = c("Interval" = "Geneid"))

# Overlap filtering on boolean columns (keep .bool suffix so ends_with('.bool') works)
bool_treated_data <- main_df %>% select(1:6, ends_with(".bool")) %>%
  mutate(interval_id = paste(chr, start, end, sep = ":"))

bool_timepoints <- c("_T10_CI","_T30_CI",
                     "_T10_DMSO","_T30_DMSO","_T60_DMSO","_T90_DMSO","_T120_DMSO",
                     "_T30_PMA","_T60_PMA","_T90_PMA","_T120_PMA")

bool_treated_overlaps <- purrr::map_dfr(bool_timepoints, ~{
  df <- bool_treated_data %>%
    select(1:6, contains(.x)) %>%
    mutate(Overlaps = rowSums(select(., ends_with(".bool")) == TRUE))
  df %>% filter(Overlaps >= num_overlaps) %>% select(interval_id)
})

filtered_treated_intervals <- bool_treated_data %>%
  filter(interval_id %in% unique(bool_treated_overlaps$interval_id))

colnames(filtered_treated_intervals)[1:4] <- c("chr", "start", "end", "Interval")
filtered_treated_intervals <- filtered_treated_intervals %>% select(-Interval)  # avoid duplicate in join

# Join scaled counts with filtered intervals by coordinates
normalized_counts <- inner_join(data_b_df, filtered_treated_intervals, by = c("chr","start","end"))

# Keep only coords, Interval (from data_b_df), and target timepoints
normalized_counts <- normalized_counts %>%
  select(chr, start, end, Interval, matches(paste(bool_timepoints, collapse = "|")))

# ─────────────────────────────────────────────────────
# Write outputs
# ─────────────────────────────────────────────────────
write.table(data_b_df,
            file = file.path(output_wd, paste0("consensus_peaks_sorted_scale_normalized_counts_", current_date, ".txt")),
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(normalized_counts,
            file = file.path(output_wd, paste0("consensus_peaks_sorted_with_", num_overlaps, "_overlaps_", current_date, ".txt")),
            sep = "\t", quote = FALSE, row.names = FALSE)

message("✅ Full script complete: Outputs written to ", output_wd)

