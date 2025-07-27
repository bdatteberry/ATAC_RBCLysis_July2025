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
theme <- theme_classic() + theme(
  text = element_text(color = "black"),
  plot.title = element_text(hjust = 0.5, size=rel(1.5)),
  axis.text.x = element_text(size=rel(1.25)), 
  axis.text.y = element_text(size=rel(1.25), angle=90, hjust=0.5),
  axis.title.y = element_text(size=rel(1.5))
)

# ─────────────────────────────────────────────────────
# Set Paths
# ─────────────────────────────────────────────────────
main_wd <- "/Users/brandiatteberry/Desktop/ATACseq_Nextflow/"; create_dir(main_wd)
macs_wd <- file.path(main_wd, "macs"); create_dir(macs_wd)
deseq_wd <- file.path(main_wd, "R"); create_dir(deseq_wd)
scale_wd <- file.path(main_wd, "scale"); create_dir(scale_wd)
output_wd <- file.path(deseq_wd, toupper(format(Sys.Date(), "%d%b%y"))); create_dir(output_wd)
homer_wd <- file.path(deseq_wd, "Homer"); create_dir(homer_wd)

current_date <- toupper(format(Sys.Date(), "%d%b%y"))
boolean_file <- "/Users/brandiatteberry/Desktop/consensus_peaks.mRp.clN.boolean.txt"
sample_map <- "/Users/brandiatteberry/Desktop/deseq_input_all.csv"

# ─────────────────────────────────────────────────────
# Boolean Filtering + UpSetR Visualization (01_)
# ─────────────────────────────────────────────────────
sample_renaming <- read.csv(sample_map, check.names=FALSE)
consensus_samples_a <- sample_renaming %>% filter(Comp_A == "Untreated") %>% distinct(analysis_name_simple, .keep_all = TRUE)
consensus_samples_b <- consensus_samples_a[, c("sample", "analysis_name_simple")]
consensus_samples <- consensus_samples_b$analysis_name_simple
consensus_peak_overlap_samples <- 3

main_df <- read.table(boolean_file, header=TRUE, sep="\t", as.is=TRUE)
bool_data <- main_df %>% select(ends_with(".bool"))
colnames(bool_data) <- gsub("\\.mRp\\.clN\\.bool$", "", colnames(bool_data))
bool_data <- bool_data %>% rename_with(~ consensus_samples_b$analysis_name_simple[match(.x, consensus_samples_b$sample)])
bool_data <- if(length(consensus_samples) == 0) {bool_data} else {bool_data %>% select(any_of(consensus_samples))}
filtered_data <- bool_data %>% mutate(across(everything(), as.integer)) %>% filter(rowSums(.) >= 1)

# UpSetR plot
png(file.path(output_wd, paste0("upsetR_plot_", current_date, ".png")), width=1800, height=1000, res=120)
upset(data = filtered_data, sets = names(filtered_data), order.by = "freq", main.bar.color = "grey", matrix.color = "#7E2945", text.scale = 1.25, keep.order = TRUE, mb.ratio = c(0.75, 0.25))
dev.off()

# Overlap bar chart
upset_overlaps <- NULL
for(i in 2:ncol(bool_data)){
  peak_number_overlaps <- bool_data %>% mutate(across(everything(), as.integer)) %>% filter(rowSums(.) >= i)
  rows_overlaps <- cbind("Overlaps"=i, "Peaks"=nrow(peak_number_overlaps))
  upset_overlaps <- rbind(upset_overlaps, rows_overlaps)
}
upset_overlap_plot <- ggplot(upset_overlaps, aes(x = Overlaps, y = Peaks)) +
  geom_bar(stat = "identity", fill = "grey") +
  geom_text(aes(label = Peaks), vjust = -0.5, color = "black", size=5) +
  labs(title = "Peaks vs. Overlaps", x = "Number of Overlaps", y = "Number of Peaks") +
  theme
ggsave(file.path(output_wd, paste0("upset_overlap_per_sample_", current_date, ".png")), plot = upset_overlap_plot, width = 6, height = 8, dpi = 500)

# BED output
consensus_peaks <- data.frame(main_df[,1:4], Height=1, Strand="+")
write.table(consensus_peaks, file.path(output_wd, paste0("consensus_peaks_", consensus_peak_overlap_samples, "samples_", current_date, ".bed")), col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")

# ─────────────────────────────────────────────────────
# featureCounts Processing for Normalization Input
# ─────────────────────────────────────────────────────
feature_counts_file <- "/Users/brandiatteberry/Desktop/consensus_peaks.mRp.clN.featureCounts.txt"
fc_raw <- read.table(feature_counts_file, header = TRUE, comment.char = "#", check.names = FALSE)
file_first_six <- file.path(scale_wd, "consensus_peaks_first_six_columns.txt")
file_bam_columns <- file.path(scale_wd, "consensus_peaks_sorted_bam_columns.txt")
write.table(fc_raw[, 1:6], file_first_six, sep = "\t", row.names = FALSE, quote = FALSE)
write.table(fc_raw[, 7:ncol(fc_raw)], file_bam_columns, sep = "\t", row.names = FALSE, quote = FALSE)
message("✅ featureCounts extraction complete: Files written to `scale` directory")

# ─────────────────────────────────────────────────────
# Generate Scale Factor Files (NEW SECTION)
# ─────────────────────────────────────────────────────
counts <- read.table(file_bam_columns, header = TRUE, sep = "\t", check.names = FALSE)
sample_totals <- colSums(counts)
scale_factors <- round(sample_totals / max(sample_totals), 5)
for (sample in names(scale_factors)) {
  fname <- file.path(scale_wd, paste0(sample, ".mRp.clN.scale_factor.txt"))
  write(scale_factors[[sample]], file = fname)
}
message("✅ Scale factor files generated to: ", scale_wd)

# ─────────────────────────────────────────────────────
# Normalization & Overlap Filtering (02_)
# ─────────────────────────────────────────────────────
bool_timepoints <- c("_T10_CI", "_T30_CI", "_T10_DMSO", "_T30_DMSO", "_T60_DMSO", "_T90_DMSO", "_T120_DMSO", "_T30_PMA", "_T60_PMA", "_T90_PMA", "_T120_PMA")
num_overlaps <- 2

file_a <- file.path(scale_wd, "consensus_peaks_first_six_columns.txt")
file_b <- file.path(scale_wd, "consensus_peaks_sorted_bam_columns.txt")
if (!file.exists(file_a) | !file.exists(file_b)) stop("Missing one or both required files for normalization")

data_a <- read.table(file_a, header = TRUE)
data_b <- read.table(file_b, header = TRUE, check.names = FALSE)
colnames(data_b) <- gsub("\\.mLb\\.clN\\.sorted\\.bam", "", colnames(data_b))
rownames(data_b) <- data_a$Geneid

scale_files <- list.files(scale_wd, pattern = "\\.scale_factor\\.txt$", full.names = TRUE)
scaling_factors_df <- bind_rows(lapply(scale_files, function(file) {
  name <- gsub("\\.mRp\\.clN\\.scale_factor\\.txt$", "", basename(file))
  scale <- read.table(file, header = FALSE, col.names = "scaling_factor")
  scale$name <- name
  scale
})) %>% select(name, scaling_factor)

new_data_b <- map_dfc(scaling_factors_df$name, function(sample) {
  if (sample %in% colnames(data_b)) {
    scale <- scaling_factors_df %>% filter(name == sample) %>% pull(scaling_factor)
    round(100 * data_b[[sample]] * scale)
  } else {
    NULL
  }
})

colnames(new_data_b) <- scaling_factors_df$name[scaling_factors_df$name %in% colnames(data_b)]
data_b_export <- cbind(Interval = rownames(new_data_b), new_data_b)
data_b_df <- left_join(data_b_export, data_a, by = c("Interval" = "Geneid"))

# Overlap filtering
bool_treated_data <- main_df %>% select(1:6, ends_with(".bool"))
bool_treated_data <- bool_treated_data %>% mutate(interval_id = paste(chr, start, end, sep = ":"))
colnames(bool_treated_data) <- gsub("\\.mRp\\.clN\\.bool", "", colnames(bool_treated_data))
bool_treated_overlaps <- map_dfr(bool_timepoints, ~ {
  df <- bool_treated_data %>% select(1:6, contains(.x)) %>% mutate(Overlaps = rowSums(select(., ends_with(".bool")) == TRUE))
  df %>% filter(Overlaps >= num_overlaps) %>% select(interval_id)
})

filtered_treated_intervals <- bool_treated_data %>% filter(interval_id %in% unique(bool_treated_overlaps$interval_id))
colnames(filtered_treated_intervals)[1:4] <- c("chr", "start", "end", "Interval")

normalized_counts <- inner_join(data_b_df, filtered_treated_intervals, by = c("chr", "start", "end"))
normalized_counts <- normalized_counts[, c("chr", "start", "end", "Interval", matches(paste(bool_timepoints, collapse = "|")))]

write.table(data_b_df,
            file = file.path(output_wd, paste0("consensus_peaks_sorted_scale_normalized_counts_", current_date, ".txt")),
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(normalized_counts,
            file = file.path(output_wd, paste0("consensus_peaks_sorted_with_", num_overlaps, "_overlaps_", current_date, ".txt")),
            sep = "\t", quote = FALSE, row.names = FALSE)

message("✅ Full script complete: Outputs written to ", output_wd)
