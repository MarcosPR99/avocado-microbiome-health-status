############################################################
# 16S FINAL SCRIPT
# SECTION 1 — BASIC DATASET METRICS
#
# Objetivos:
# 1) Resumir el número de muestras por mes, compartimento y estado sanitario
# 2) Resumir la profundidad de secuenciación
# 3) Contar ASVs detectados por grupo
# 4) Contar géneros detectados por grupo
# 5) Exportar tablas limpias para Methods / Results / Supplementary
############################################################

suppressPackageStartupMessages({
  library(phyloseq)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
})

############################################################
# 0) WORKING DIRECTORY, INPUTS AND OUTPUT FOLDER
############################################################

setwd("path/to/QIIME2_files_to_R_16S")

# Input files
physeq_path <- "Phyloseq_16S.rds"
metadata_path <- "metadata_samples.tsv"

# Output folder for this section
outdir_basic <- "00_basic_metrics_16S"
if (!dir.exists(outdir_basic)) dir.create(outdir_basic, recursive = TRUE)

############################################################
# 1) LOAD PHYLOSEQ AND BASIC CHECKS
############################################################

ps <- readRDS(physeq_path)
ps <- prune_samples(sample_sums(ps) > 0, ps)

sample_df <- data.frame(sample_data(ps))
sample_df$SampleID <- rownames(sample_df)
sample_df$Reads <- sample_sums(ps)

cat("\n====================================================\n")
cat("BASIC CHECKS\n")
cat("====================================================\n")
cat("Working directory:\n")
print(getwd())

cat("\nFiles in working directory:\n")
print(list.files())

cat("\nNumber of samples in phyloseq:\n")
print(nsamples(ps))

cat("\nNumber of ASVs in phyloseq:\n")
print(ntaxa(ps))

cat("\nSample-data columns:\n")
print(colnames(sample_df))

required_cols <- c("Month", "Compartment", "Health_status")
missing_cols <- setdiff(required_cols, colnames(sample_df))

if (length(missing_cols) > 0) {
  stop("Faltan columnas requeridas en sample_data: ", paste(missing_cols, collapse = ", "))
}

############################################################
# 2) CLEAN FACTOR LEVELS AND ORDER
############################################################

month_levels <- c("July", "October")

compartment_levels <- c(
  "Bulk soil",
  "Rhizosphere",
  "Endosphere",
  "Phyllosphere"
)

health_levels <- c(
  "Bulk soil",
  "Healthy Rhizosphere", "Diseased Rhizosphere",
  "Healthy Endosphere",  "Diseased Endosphere",
  "Healthy Phyllosphere","Diseased Phyllosphere"
)

sample_df$Month <- factor(as.character(sample_df$Month), levels = month_levels)
sample_df$Compartment <- factor(as.character(sample_df$Compartment), levels = compartment_levels)
sample_df$Health_status <- factor(as.character(sample_df$Health_status), levels = health_levels)

cat("\n====================================================\n")
cat("LEVEL COUNTS\n")
cat("====================================================\n")
cat("Month:\n")
print(table(sample_df$Month, useNA = "ifany"))

cat("\nCompartment:\n")
print(table(sample_df$Compartment, useNA = "ifany"))

cat("\nHealth_status:\n")
print(table(sample_df$Health_status, useNA = "ifany"))

############################################################
# 3) SAMPLE COUNTS
############################################################

sample_count_by_group <- sample_df %>%
  group_by(Month, Compartment, Health_status) %>%
  summarise(
    n_samples = n(),
    .groups = "drop"
  ) %>%
  arrange(Month, Compartment, Health_status)

sample_count_by_compartment <- sample_df %>%
  group_by(Month, Compartment) %>%
  summarise(
    n_samples = n(),
    .groups = "drop"
  ) %>%
  arrange(Month, Compartment)

sample_count_by_month <- sample_df %>%
  group_by(Month) %>%
  summarise(
    n_samples = n(),
    .groups = "drop"
  ) %>%
  arrange(Month)

write_csv(sample_count_by_group,
          file.path(outdir_basic, "01_sample_count_by_month_compartment_health.csv"))
write_csv(sample_count_by_compartment,
          file.path(outdir_basic, "02_sample_count_by_month_compartment.csv"))
write_csv(sample_count_by_month,
          file.path(outdir_basic, "03_sample_count_by_month.csv"))

############################################################
# 4) SEQUENCING DEPTH SUMMARY
############################################################

reads_by_group <- sample_df %>%
  group_by(Month, Compartment, Health_status) %>%
  summarise(
    n_samples    = n(),
    total_reads  = sum(Reads, na.rm = TRUE),
    mean_reads   = mean(Reads, na.rm = TRUE),
    median_reads = median(Reads, na.rm = TRUE),
    min_reads    = min(Reads, na.rm = TRUE),
    max_reads    = max(Reads, na.rm = TRUE),
    sd_reads     = sd(Reads, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(Month, Compartment, Health_status)

reads_by_compartment <- sample_df %>%
  group_by(Month, Compartment) %>%
  summarise(
    n_samples    = n(),
    total_reads  = sum(Reads, na.rm = TRUE),
    mean_reads   = mean(Reads, na.rm = TRUE),
    median_reads = median(Reads, na.rm = TRUE),
    min_reads    = min(Reads, na.rm = TRUE),
    max_reads    = max(Reads, na.rm = TRUE),
    sd_reads     = sd(Reads, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(Month, Compartment)

reads_by_month <- sample_df %>%
  group_by(Month) %>%
  summarise(
    n_samples    = n(),
    total_reads  = sum(Reads, na.rm = TRUE),
    mean_reads   = mean(Reads, na.rm = TRUE),
    median_reads = median(Reads, na.rm = TRUE),
    min_reads    = min(Reads, na.rm = TRUE),
    max_reads    = max(Reads, na.rm = TRUE),
    sd_reads     = sd(Reads, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(Month)

write_csv(reads_by_group,
          file.path(outdir_basic, "04_reads_summary_by_month_compartment_health.csv"))
write_csv(reads_by_compartment,
          file.path(outdir_basic, "05_reads_summary_by_month_compartment.csv"))
write_csv(reads_by_month,
          file.path(outdir_basic, "06_reads_summary_by_month.csv"))

############################################################
# 5) PREPARE OTU TABLE AND TAXONOMY TABLE
############################################################

otu_mat <- as(otu_table(ps), "matrix")
if (!taxa_are_rows(ps)) otu_mat <- t(otu_mat)

tax_df <- as.data.frame(tax_table(ps))
tax_df$ASV_ID <- rownames(tax_df)

cat("\n====================================================\n")
cat("TAXONOMY TABLE CHECK\n")
cat("====================================================\n")
print(colnames(tax_df))

if (!"Genus" %in% colnames(tax_df)) {
  stop("No encuentro la columna 'Genus' en tax_table(ps).")
}

############################################################
# 6) HELPER FUNCTIONS
############################################################

count_detected_asvs <- function(sample_ids, otu_matrix) {
  submat <- otu_matrix[, sample_ids, drop = FALSE]
  sum(rowSums(submat > 0) > 0)
}

count_detected_genera <- function(sample_ids, otu_matrix, taxonomy_df) {
  submat <- otu_matrix[, sample_ids, drop = FALSE]
  present_asvs <- rownames(submat)[rowSums(submat > 0) > 0]
  
  taxonomy_df %>%
    filter(ASV_ID %in% present_asvs) %>%
    filter(!is.na(Genus), Genus != "") %>%
    summarise(n_genera = n_distinct(Genus)) %>%
    pull(n_genera)
}

############################################################
# 7) GROUPING TABLE
############################################################

grouping_df <- sample_df %>%
  select(SampleID, Month, Compartment, Health_status)

############################################################
# 8) HELPER TO SUMMARIZE DETECTED ASVS AND GENERA
############################################################

summarise_detected_features <- function(df_subset, otu_matrix, taxonomy_df) {
  sample_ids <- df_subset$SampleID
  
  tibble(
    n_samples = length(sample_ids),
    n_ASVs_detected = count_detected_asvs(sample_ids, otu_matrix),
    n_genera_detected = count_detected_genera(sample_ids, otu_matrix, taxonomy_df)
  )
}

############################################################
# 9) ASV AND GENUS COUNTS BY MONTH + COMPARTMENT + HEALTH
############################################################

asv_genus_by_group <- grouping_df %>%
  group_by(Month, Compartment, Health_status) %>%
  group_modify(~ summarise_detected_features(.x, otu_mat, tax_df)) %>%
  ungroup() %>%
  arrange(Month, Compartment, Health_status)

write_csv(
  asv_genus_by_group,
  file.path(outdir_basic, "07_asv_genus_counts_by_month_compartment_health.csv")
)

############################################################
# 10) ASV AND GENUS COUNTS BY MONTH + COMPARTMENT
############################################################

asv_genus_by_compartment <- grouping_df %>%
  group_by(Month, Compartment) %>%
  group_modify(~ summarise_detected_features(.x, otu_mat, tax_df)) %>%
  ungroup() %>%
  arrange(Month, Compartment)

write_csv(
  asv_genus_by_compartment,
  file.path(outdir_basic, "08_asv_genus_counts_by_month_compartment.csv")
)

############################################################
# 11) ASV AND GENUS COUNTS BY MONTH
############################################################

asv_genus_by_month <- grouping_df %>%
  group_by(Month) %>%
  group_modify(~ summarise_detected_features(.x, otu_mat, tax_df)) %>%
  ungroup() %>%
  arrange(Month)

write_csv(
  asv_genus_by_month,
  file.path(outdir_basic, "09_asv_genus_counts_by_month.csv")
)

############################################################
# 12) MERGED FINAL SUMMARY TABLES
############################################################

summary_group_final <- reads_by_group %>%
  left_join(asv_genus_by_group,
            by = c("Month", "Compartment", "Health_status", "n_samples")) %>%
  arrange(Month, Compartment, Health_status)

summary_compartment_final <- reads_by_compartment %>%
  left_join(asv_genus_by_compartment,
            by = c("Month", "Compartment", "n_samples")) %>%
  arrange(Month, Compartment)

summary_month_final <- reads_by_month %>%
  left_join(asv_genus_by_month,
            by = c("Month", "n_samples")) %>%
  arrange(Month)

write_csv(summary_group_final,
          file.path(outdir_basic, "11_final_summary_by_month_compartment_health.csv"))
write_csv(summary_compartment_final,
          file.path(outdir_basic, "12_final_summary_by_month_compartment.csv"))
write_csv(summary_month_final,
          file.path(outdir_basic, "13_final_summary_by_month.csv"))

############################################################
# 13) OPTIONAL METHODS/RESULTS-READY TABLE
############################################################

methods_table <- summary_group_final %>%
  mutate(
    reads_mean_range = paste0(
      round(mean_reads, 1), " (",
      min_reads, "–", max_reads, ")"
    )
  ) %>%
  select(
    Month,
    Compartment,
    Health_status,
    n_samples,
    total_reads,
    mean_reads,
    median_reads,
    min_reads,
    max_reads,
    n_ASVs_detected,
    n_genera_detected,
    reads_mean_range
  )

write_csv(methods_table,
          file.path(outdir_basic, "14_methods_results_ready_table.csv"))

############################################################
# 14) PRINT RESULTS TO CONSOLE
############################################################

global_totals <- tibble(
  n_samples_total = nsamples(ps),
  total_reads = sum(sample_sums(ps)),
  total_ASVs = ntaxa(ps),
  total_genera = tax_df %>%
    filter(!is.na(Genus), Genus != "") %>%
    summarise(n_genera = n_distinct(Genus)) %>%
    pull(n_genera)
)

write_csv(
  global_totals,
  file.path(outdir_basic, "10_global_totals_16S.csv")
)

cat("\n====================================================\n")
cat("GLOBAL TOTALS\n")
cat("====================================================\n")
print(global_totals)

cat("\n====================================================\n")
cat("SUMMARY BY MONTH\n")
cat("====================================================\n")
print(summary_month_final, n = Inf)

cat("\n====================================================\n")
cat("SUMMARY BY MONTH AND COMPARTMENT\n")
cat("====================================================\n")
print(summary_compartment_final, n = Inf)

cat("\n====================================================\n")
cat("SUMMARY BY MONTH, COMPARTMENT, AND HEALTH STATUS\n")
cat("====================================================\n")
print(summary_group_final, n = Inf)

cat("\n====================================================\n")
cat("MANUSCRIPT-STYLE SENTENCES BY MONTH AND COMPARTMENT\n")
cat("====================================================\n")

for (i in seq_len(nrow(summary_compartment_final))) {
  row_i <- summary_compartment_final[i, ]
  
  cat(
    as.character(row_i$Month), " — ",
    as.character(row_i$Compartment), ": ",
    row_i$n_samples, " samples, ",
    row_i$total_reads, " reads, ",
    row_i$n_ASVs_detected, " detected ASVs, ",
    row_i$n_genera_detected, " detected genera.\n",
    sep = ""
  )
}

cat("\n====================================================\n")
cat("FILES WRITTEN TO:\n")
cat("====================================================\n")
print(normalizePath(outdir_basic))





############################################################
# SECTION 1B — BASIC METRICS PER SAMPLE
#
# Objetivos:
# 1) Calcular reads por muestra
# 2) Calcular ASVs detectados por muestra
# 3) Calcular géneros detectados por muestra
# 4) Exportar tablas por muestra y resúmenes por grupo
############################################################

############################################################
# 15) HELPER FUNCTIONS FOR PER-SAMPLE COUNTS
############################################################

count_asvs_one_sample <- function(sample_id, otu_matrix) {
  sum(otu_matrix[, sample_id] > 0, na.rm = TRUE)
}

count_genera_one_sample <- function(sample_id, otu_matrix, taxonomy_df) {
  present_asvs <- rownames(otu_matrix)[otu_matrix[, sample_id] > 0]
  
  taxonomy_df %>%
    filter(ASV_ID %in% present_asvs) %>%
    filter(!is.na(Genus), Genus != "") %>%
    summarise(n_genera = n_distinct(Genus)) %>%
    pull(n_genera)
}

############################################################
# 16) PER-SAMPLE BASIC METRICS
############################################################

per_sample_basic_metrics <- sample_df %>%
  mutate(
    n_ASVs_detected = vapply(
      SampleID,
      count_asvs_one_sample,
      integer(1),
      otu_matrix = otu_mat
    ),
    n_genera_detected = vapply(
      SampleID,
      count_genera_one_sample,
      integer(1),
      otu_matrix = otu_mat,
      taxonomy_df = tax_df
    )
  ) %>%
  select(
    SampleID,
    Month,
    Compartment,
    Health_status,
    Reads,
    n_ASVs_detected,
    n_genera_detected
  ) %>%
  arrange(Month, Compartment, Health_status, SampleID)

write_csv(
  per_sample_basic_metrics,
  file.path(outdir_basic, "15_per_sample_basic_metrics.csv")
)

############################################################
# 17) SUMMARY OF PER-SAMPLE METRICS BY MONTH + COMPARTMENT + HEALTH
############################################################

per_sample_summary_by_group <- per_sample_basic_metrics %>%
  group_by(Month, Compartment, Health_status) %>%
  summarise(
    n_samples = n(),
    
    reads_mean   = mean(Reads, na.rm = TRUE),
    reads_median = median(Reads, na.rm = TRUE),
    reads_min    = min(Reads, na.rm = TRUE),
    reads_max    = max(Reads, na.rm = TRUE),
    reads_sd     = sd(Reads, na.rm = TRUE),
    
    ASVs_mean    = mean(n_ASVs_detected, na.rm = TRUE),
    ASVs_median  = median(n_ASVs_detected, na.rm = TRUE),
    ASVs_min     = min(n_ASVs_detected, na.rm = TRUE),
    ASVs_max     = max(n_ASVs_detected, na.rm = TRUE),
    ASVs_sd      = sd(n_ASVs_detected, na.rm = TRUE),
    
    genera_mean   = mean(n_genera_detected, na.rm = TRUE),
    genera_median = median(n_genera_detected, na.rm = TRUE),
    genera_min    = min(n_genera_detected, na.rm = TRUE),
    genera_max    = max(n_genera_detected, na.rm = TRUE),
    genera_sd     = sd(n_genera_detected, na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  arrange(Month, Compartment, Health_status)

write_csv(
  per_sample_summary_by_group,
  file.path(outdir_basic, "16_per_sample_summary_by_month_compartment_health.csv")
)

############################################################
# 18) SUMMARY OF PER-SAMPLE METRICS BY MONTH + COMPARTMENT
############################################################

per_sample_summary_by_compartment <- per_sample_basic_metrics %>%
  group_by(Month, Compartment) %>%
  summarise(
    n_samples = n(),
    
    reads_mean   = mean(Reads, na.rm = TRUE),
    reads_median = median(Reads, na.rm = TRUE),
    reads_min    = min(Reads, na.rm = TRUE),
    reads_max    = max(Reads, na.rm = TRUE),
    reads_sd     = sd(Reads, na.rm = TRUE),
    
    ASVs_mean    = mean(n_ASVs_detected, na.rm = TRUE),
    ASVs_median  = median(n_ASVs_detected, na.rm = TRUE),
    ASVs_min     = min(n_ASVs_detected, na.rm = TRUE),
    ASVs_max     = max(n_ASVs_detected, na.rm = TRUE),
    ASVs_sd      = sd(n_ASVs_detected, na.rm = TRUE),
    
    genera_mean   = mean(n_genera_detected, na.rm = TRUE),
    genera_median = median(n_genera_detected, na.rm = TRUE),
    genera_min    = min(n_genera_detected, na.rm = TRUE),
    genera_max    = max(n_genera_detected, na.rm = TRUE),
    genera_sd     = sd(n_genera_detected, na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  arrange(Month, Compartment)

write_csv(
  per_sample_summary_by_compartment,
  file.path(outdir_basic, "17_per_sample_summary_by_month_compartment.csv")
)

############################################################
# 19) OPTIONAL OUTLIER CHECK TABLE
############################################################

per_sample_outlier_flags <- per_sample_basic_metrics %>%
  group_by(Month, Compartment, Health_status) %>%
  mutate(
    reads_z = ifelse(sd(Reads, na.rm = TRUE) > 0,
                     (Reads - mean(Reads, na.rm = TRUE)) / sd(Reads, na.rm = TRUE),
                     0),
    ASVs_z = ifelse(sd(n_ASVs_detected, na.rm = TRUE) > 0,
                    (n_ASVs_detected - mean(n_ASVs_detected, na.rm = TRUE)) / sd(n_ASVs_detected, na.rm = TRUE),
                    0),
    genera_z = ifelse(sd(n_genera_detected, na.rm = TRUE) > 0,
                      (n_genera_detected - mean(n_genera_detected, na.rm = TRUE)) / sd(n_genera_detected, na.rm = TRUE),
                      0),
    potential_outlier = ifelse(
      abs(reads_z) >= 3 | abs(ASVs_z) >= 3 | abs(genera_z) >= 3,
      "YES",
      "NO"
    )
  ) %>%
  ungroup() %>%
  arrange(desc(potential_outlier), Month, Compartment, Health_status, SampleID)

write_csv(
  per_sample_outlier_flags,
  file.path(outdir_basic, "18_per_sample_outlier_flags.csv")
)

############################################################
# 20) PRINT TO CONSOLE
############################################################

cat("\n====================================================\n")
cat("PER-SAMPLE BASIC METRICS (first 20 rows)\n")
cat("====================================================\n")
print(dplyr::slice_head(per_sample_basic_metrics, n = 20), width = Inf)

cat("\n====================================================\n")
cat("PER-SAMPLE SUMMARY BY MONTH + COMPARTMENT + HEALTH STATUS\n")
cat("====================================================\n")
print(per_sample_summary_by_group, n = Inf)

cat("\n====================================================\n")
cat("PER-SAMPLE SUMMARY BY MONTH + COMPARTMENT\n")
cat("====================================================\n")
print(per_sample_summary_by_compartment, n = Inf)

cat("\n====================================================\n")
cat("POTENTIAL OUTLIERS\n")
cat("====================================================\n")
print(
  per_sample_outlier_flags %>%
    filter(potential_outlier == "YES"),
  n = Inf
)




###############################################################################################################################################################
###############################################################################################################################################################
###############################################################################################################################################################
###############################################################################################################################################################
###############################################################################################################################################################




############################################################
# SECTION 2 — ALPHA DIVERSITY (RAREFIED DATA)
#
# Objetivos:
# 1) Extraer métricas de alpha diversity desde el archivo
#    alpha-rarefaction_withPhy_depth_18288.qzv
# 2) Construir una tabla por muestra para July + October
# 3) Resumir las métricas por mes, compartimento y estado sanitario
# 4) Realizar comparaciones estadísticas de Healthy vs Diseased
#    dentro de cada compartimento y mes
# 5) Exportar tablas y figuras listas para el manuscrito
############################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(multcompView)
  library(patchwork)
  library(readr)
})

############################################################
# 0) INPUTS AND OUTPUT FOLDERS
############################################################

qzv_file <- "alpha-rarefaction_withPhy_depth_18288.qzv"
alpha_depth_target <- 18288

outdir_alpha <- "01_alpha_diversity_16S"
if (!dir.exists(outdir_alpha)) dir.create(outdir_alpha, recursive = TRUE)

############################################################
# 1) BASIC CHECKS
############################################################

cat("\n====================================================\n")
cat("ALPHA DIVERSITY — BASIC CHECKS\n")
cat("====================================================\n")
cat("Working directory:\n")
print(getwd())

cat("\nQZV file exists?\n")
print(file.exists(qzv_file))

if (!file.exists(qzv_file)) {
  stop("No encuentro el archivo: ", qzv_file)
}

zip_listing <- tryCatch(
  unzip(qzv_file, list = TRUE),
  error = function(e) e
)

cat("\nZIP LISTING TEST\n")
print(zip_listing)

if (inherits(zip_listing, "error")) {
  stop("R no puede leer el .qzv como archivo zip. Revisa el archivo.")
}

############################################################
# 2) HELPERS
############################################################

norm_names <- function(n) {
  n <- tolower(n)
  n <- gsub("^#+", "", n)
  n <- gsub("[[:space:]]+", "_", n)
  n <- gsub("-", "_", n)
  n
}

safe_num <- function(x) {
  if (is.numeric(x)) return(x)
  x <- as.character(x)
  x <- gsub("\\s", "", x)
  
  has_both <- grepl(",", x) & grepl("\\.", x)
  x[has_both] <- gsub("\\.", "", x[has_both])
  x[has_both] <- sub(",", ".", x[has_both], fixed = TRUE)
  
  only_comma <- grepl(",", x) & !grepl("\\.", x)
  x[only_comma] <- sub(",", ".", x[only_comma], fixed = TRUE)
  
  many_dots <- vapply(regmatches(x, gregexpr("\\.", x)), length, integer(1)) > 1
  x[many_dots] <- gsub("\\.", "", x[many_dots])
  
  suppressWarnings(as.numeric(x))
}

digits_by_metric <- function(metric) dplyr::case_when(
  grepl("Shannon", metric, ignore.case = TRUE) ~ 2L,
  grepl("Faith", metric, ignore.case = TRUE) ~ 1L,
  grepl("Chao1|Observed", metric, ignore.case = TRUE) ~ 0L,
  TRUE ~ 2L
)

sig_mark <- function(p) ifelse(
  is.na(p), NA_character_,
  ifelse(p < 0.001, "***",
         ifelse(p < 0.01, "**",
                ifelse(p < 0.05, "*", "ns")))
)

label_stacked <- function(x) gsub(" ", "\n", x)

############################################################
# 3) READ ONE METRIC TABLE FROM QZV
############################################################

read_metric_csv <- function(qzv, filename_regex) {
  tmpdir <- file.path(tempdir(), "alpha_raref_tmp_16S")
  
  if (dir.exists(tmpdir)) unlink(tmpdir, recursive = TRUE, force = TRUE)
  dir.create(tmpdir, recursive = TRUE, showWarnings = FALSE)
  
  unzip(qzv, exdir = tmpdir)
  
  files <- list.files(
    tmpdir,
    pattern = filename_regex,
    recursive = TRUE,
    full.names = TRUE
  )
  
  cat("\nSearching metric file with pattern:", filename_regex, "\n")
  print(files)
  
  if (length(files) < 1) {
    stop("No se encontró ningún archivo dentro del qzv para el patrón: ", filename_regex)
  }
  
  prefer <- grepl("[/\\\\]data[/\\\\][^/\\\\]*$", files)
  csv_file <- if (any(prefer)) files[which(prefer)[1]] else files[1]
  
  cat("Using file:\n")
  print(csv_file)
  
  out <- readr::read_csv(csv_file, show_col_types = FALSE)
  names(out) <- norm_names(names(out))
  out
}

############################################################
# 4) DEFINE METRICS TO EXTRACT
############################################################

metrics <- list(
  "Observed features"   = "observed_features\\.csv$",
  "Chao1"               = "chao1\\.csv$",
  "Faith's PD"          = "faith_pd\\.csv$",
  "Shannon's Diversity" = "shannon\\.csv$"
)

############################################################
# 5) EXTRACT ONE METRIC FOR JULY + OCTOBER
############################################################

df_for_metric_all_months <- function(metric_regex, metric_name, target_depth) {
  wide <- read_metric_csv(qzv_file, metric_regex)
  
  cat("\nColumns detected for", metric_name, ":\n")
  print(names(wide))
  
  need_cols <- c("sample_id", "name", "health_status", "compartment", "month")
  missing_cols <- setdiff(need_cols, names(wide))
  if (length(missing_cols) > 0) {
    stop(
      "Faltan columnas requeridas en la tabla de alpha diversity para ",
      metric_name, ": ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  depth_iter_cols <- grep("^depth_[0-9]+_iter_[0-9]+$", names(wide), value = TRUE)
  
  if (length(depth_iter_cols) == 0) {
    stop("No se encontraron columnas depth_*_iter_* para ", metric_name)
  }
  
  long <- wide %>%
    pivot_longer(
      cols = all_of(depth_iter_cols),
      names_to = c("depth", "iter"),
      names_pattern = "^depth_([0-9]+)_iter_([0-9]+)$",
      values_to = "value"
    ) %>%
    mutate(
      depth = as.numeric(depth),
      iter  = as.integer(iter)
    ) %>%
    filter(tolower(month) %in% c("july", "october"))
  
  available_depths <- sort(unique(long$depth))
  depth_use <- available_depths[which.min(abs(available_depths - target_depth))]
  
  cat("\nAvailable depths for", metric_name, ":\n")
  print(available_depths)
  cat("Target depth:", target_depth, "\n")
  cat("Depth used:", depth_use, "\n")
  
  alpha <- long %>%
    filter(depth == depth_use) %>%
    group_by(sample_id, name, health_status, compartment, month) %>%
    summarise(alpha_value = mean(value, na.rm = TRUE), .groups = "drop") %>%
    mutate(
      Month = factor(month, levels = c("July", "October")),
      Compartment = factor(compartment,
                           levels = c("Bulk soil", "Rhizosphere", "Endosphere", "Phyllosphere")),
      Health_status = factor(
        health_status,
        levels = c(
          "Bulk soil",
          "Healthy Rhizosphere", "Diseased Rhizosphere",
          "Healthy Endosphere",  "Diseased Endosphere",
          "Healthy Phyllosphere","Diseased Phyllosphere"
        )
      ),
      Metric = metric_name,
      alpha_value = safe_num(alpha_value)
    ) %>%
    select(
      sample_id,
      Metric,
      Month,
      Compartment,
      Health_status,
      alpha_value
    )
  
  alpha
}

############################################################
# 6) BUILD FULL LONG TABLE
############################################################

alpha_values_list <- list()

for (mtr in names(metrics)) {
  cat("\n========================================\n")
  cat("PROCESSING:", mtr, "\n")
  cat("========================================\n")
  
  alpha_values_list[[mtr]] <- df_for_metric_all_months(
    metric_regex = metrics[[mtr]],
    metric_name = mtr,
    target_depth = alpha_depth_target
  )
}

alpha_values_long <- bind_rows(alpha_values_list) %>%
  arrange(Metric, Month, Compartment, Health_status, sample_id)

write_csv(
  alpha_values_long,
  file.path(outdir_alpha, "01_alpha_values_long_July_October.csv")
)

############################################################
# 7) SUMMARY TABLES
############################################################

alpha_summary_by_group <- alpha_values_long %>%
  group_by(Metric, Month, Compartment, Health_status) %>%
  summarise(
    n       = n(),
    mean    = mean(alpha_value, na.rm = TRUE),
    sd      = sd(alpha_value, na.rm = TRUE),
    se      = sd / sqrt(n),
    median  = median(alpha_value, na.rm = TRUE),
    q1      = quantile(alpha_value, 0.25, na.rm = TRUE),
    q3      = quantile(alpha_value, 0.75, na.rm = TRUE),
    min     = min(alpha_value, na.rm = TRUE),
    max     = max(alpha_value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(Metric, Month, Compartment, Health_status)

write_csv(
  alpha_summary_by_group,
  file.path(outdir_alpha, "02_alpha_summary_by_group_July_October.csv")
)

alpha_summary_by_compartment <- alpha_values_long %>%
  group_by(Metric, Month, Compartment) %>%
  summarise(
    n       = n(),
    mean    = mean(alpha_value, na.rm = TRUE),
    sd      = sd(alpha_value, na.rm = TRUE),
    se      = sd / sqrt(n),
    median  = median(alpha_value, na.rm = TRUE),
    q1      = quantile(alpha_value, 0.25, na.rm = TRUE),
    q3      = quantile(alpha_value, 0.75, na.rm = TRUE),
    min     = min(alpha_value, na.rm = TRUE),
    max     = max(alpha_value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(Metric, Month, Compartment)

write_csv(
  alpha_summary_by_compartment,
  file.path(outdir_alpha, "03_alpha_summary_by_compartment_July_October.csv")
)

############################################################
# 8) GLOBAL ANOVA BY MONTH AND METRIC
############################################################

run_global_alpha_anova <- function(df_metric_month) {
  fit <- aov(alpha_value ~ Health_status, data = df_metric_month)
  anova_tbl <- summary(fit)[[1]]
  
  tibble(
    df1 = anova_tbl$Df[1],
    df2 = anova_tbl$Df[2],
    F_value = anova_tbl$`F value`[1],
    p_value = anova_tbl$`Pr(>F)`[1]
  )
}

alpha_global_anova <- alpha_values_long %>%
  group_by(Metric, Month) %>%
  group_modify(~ run_global_alpha_anova(.x)) %>%
  ungroup()

write_csv(
  alpha_global_anova,
  file.path(outdir_alpha, "04_alpha_global_anova_by_month.csv")
)

############################################################
# 9) PAIRWISE TUKEY BY MONTH AND METRIC
############################################################

run_global_tukey <- function(df_metric_month) {
  fit <- aov(alpha_value ~ Health_status, data = df_metric_month)
  tk  <- TukeyHSD(fit)$Health_status
  
  tibble::as_tibble(tk, rownames = "contrast") %>%
    transmute(
      contrast,
      estimate = diff,
      conf.low = lwr,
      conf.high = upr,
      p_adj = `p adj`
    )
}

alpha_global_tukey <- alpha_values_long %>%
  group_by(Metric, Month) %>%
  group_modify(~ run_global_tukey(.x)) %>%
  ungroup() %>%
  mutate(
    sig = sig_mark(p_adj)
  )

write_csv(
  alpha_global_tukey,
  file.path(outdir_alpha, "05_alpha_global_tukey_by_month.csv")
)

############################################################
# 10) HEALTHY VS DISEASED WITHIN EACH COMPARTMENT AND MONTH
############################################################

extract_two_group_test <- function(df_sub) {
  fit <- aov(alpha_value ~ Health_status, data = df_sub)
  tk  <- TukeyHSD(fit)$Health_status
  
  if (nrow(tk) != 1) {
    return(tibble(
      contrast = NA_character_,
      estimate = NA_real_,
      conf.low = NA_real_,
      conf.high = NA_real_,
      p_adj = NA_real_
    ))
  }
  
  tibble(
    contrast = rownames(tk)[1],
    estimate = tk[1, "diff"],
    conf.low = tk[1, "lwr"],
    conf.high = tk[1, "upr"],
    p_adj = tk[1, "p adj"]
  )
}

alpha_within_compartment_tests <- alpha_values_long %>%
  filter(Compartment != "Bulk soil") %>%
  group_by(Metric, Month, Compartment) %>%
  group_modify(~ extract_two_group_test(.x)) %>%
  ungroup() %>%
  mutate(
    sig = sig_mark(p_adj)
  ) %>%
  arrange(Metric, Month, Compartment)

write_csv(
  alpha_within_compartment_tests,
  file.path(outdir_alpha, "06_alpha_within_compartment_healthy_vs_diseased.csv")
)

############################################################
# 11) OPTIONAL WILCOXON VALIDATION WITHIN COMPARTMENT
############################################################

run_wilcox_two_groups <- function(df_sub) {
  df_sub <- df_sub %>%
    filter(!is.na(alpha_value), !is.na(Health_status)) %>%
    mutate(Health_status = droplevels(as.factor(Health_status)))
  
  present_groups <- unique(as.character(df_sub$Health_status))
  
  if (length(present_groups) != 2) {
    return(tibble(
      group1 = NA_character_,
      group2 = NA_character_,
      wilcox_p = NA_real_
    ))
  }
  
  x1 <- df_sub$alpha_value[df_sub$Health_status == present_groups[1]]
  x2 <- df_sub$alpha_value[df_sub$Health_status == present_groups[2]]
  
  pval <- tryCatch(
    wilcox.test(x1, x2, exact = FALSE)$p.value,
    error = function(e) NA_real_
  )
  
  tibble(
    group1 = present_groups[1],
    group2 = present_groups[2],
    wilcox_p = pval
  )
}

alpha_within_compartment_wilcox <- alpha_values_long %>%
  filter(Compartment != "Bulk soil") %>%
  group_by(Metric, Month, Compartment) %>%
  group_modify(~ run_wilcox_two_groups(.x)) %>%
  ungroup() %>%
  arrange(Metric, Month, Compartment)

write_csv(
  alpha_within_compartment_wilcox,
  file.path(outdir_alpha, "07_alpha_within_compartment_wilcox_validation.csv")
)


############################################################
# 12) PLOT FUNCTION
############################################################

hs_colors <- c(
  "Bulk soil"              = "#896054",
  "Healthy Rhizosphere"    = "#4B9EA7",
  "Diseased Rhizosphere"   = "#801723",
  "Healthy Endosphere"     = "#6B9863",
  "Diseased Endosphere"    = "#E7C313",
  "Healthy Phyllosphere"   = "#4A5263",
  "Diseased Phyllosphere"  = "#6B4EA7"
)

plot_alpha_metric_month <- function(metric_name, month_name) {
  dfm <- alpha_values_long %>%
    filter(Metric == metric_name, Month == month_name)
  
  sumdf <- dfm %>%
    group_by(Health_status) %>%
    summarise(
      n    = n(),
      mean = mean(alpha_value, na.rm = TRUE),
      sd   = sd(alpha_value, na.rm = TRUE),
      se   = sd / sqrt(n),
      ymax_group = max(alpha_value, na.rm = TRUE),
      .groups = "drop"
    )
  
  fit <- aov(alpha_value ~ Health_status, data = dfm)
  tk  <- TukeyHSD(fit)
  cld <- multcompView::multcompLetters4(fit, tk)
  
  letters_df <- data.frame(
    Health_status = names(cld$Health_status$Letters),
    letters = cld$Health_status$Letters
  ) %>%
    left_join(sumdf %>% select(Health_status, ymax_group), by = "Health_status") %>%
    mutate(y_pos = 1.08 * ymax_group)
  
  y_top <- max(letters_df$y_pos, na.rm = TRUE)
  
  ggplot(dfm, aes(x = Health_status, y = alpha_value, fill = Health_status)) +
    geom_violin(trim = FALSE, linewidth = 0.6, alpha = 0.9) +
    geom_jitter(
      width = 0.18, size = 1.8, shape = 21, stroke = 0.4,
      color = "black", alpha = 0.85
    ) +
    geom_errorbar(
      data = sumdf,
      aes(x = Health_status, y = mean, ymin = mean - se, ymax = mean + se),
      width = 0.15, linewidth = 0.7, inherit.aes = FALSE
    ) +
    geom_point(
      data = sumdf,
      aes(x = Health_status, y = mean),
      size = 2.8, shape = 23, color = "black", fill = "white",
      stroke = 0.6, inherit.aes = FALSE
    ) +
    geom_text(
      data = letters_df,
      aes(x = Health_status, y = y_pos, label = letters),
      vjust = -0.2, size = 5, fontface = "bold",
      inherit.aes = FALSE
    ) +
    scale_fill_manual(values = hs_colors, breaks = names(hs_colors), drop = FALSE) +
    scale_x_discrete(labels = label_stacked) +
    geom_vline(xintercept = c(1.5, 3.5, 5.5), linetype = "dashed", color = "grey30") +
    coord_cartesian(ylim = c(NA, y_top * 1.08), clip = "off") +
    labs(
      x = NULL,
      y = metric_name,
      subtitle = paste0(month_name, " — mean ± SE; Tukey letters")
    ) +
    theme_bw() +
    theme(
      legend.position = "none",
      axis.title.y = element_text(size = 15, face = "bold"),
      axis.text.x = element_text(size = 11, lineheight = 0.9),
      panel.grid = element_blank(),
      plot.margin = margin(10, 10, 20, 10, "pt")
    )
}

############################################################
# 13) BUILD PANELS FOR JULY AND OCTOBER
############################################################

make_alpha_panel_for_month <- function(month_name) {
  (plot_alpha_metric_month("Observed features", month_name) +
     plot_alpha_metric_month("Chao1", month_name)) /
    (plot_alpha_metric_month("Shannon's Diversity", month_name) +
       plot_alpha_metric_month("Faith's PD", month_name)) +
    patchwork::plot_annotation(
      title = paste0("16S alpha diversity — ", month_name),
      tag_levels = "a",
      tag_prefix = "(",
      tag_suffix = ")"
    )
}

alpha_panel_july <- make_alpha_panel_for_month("July")
alpha_panel_october <- make_alpha_panel_for_month("October")

print(alpha_panel_july)
print(alpha_panel_october)

ggsave(
  file.path(outdir_alpha, "08_alpha_panel_July.png"),
  plot = alpha_panel_july,
  width = 14, height = 10, dpi = 600
)
ggsave(
  file.path(outdir_alpha, "08_alpha_panel_July.tiff"),
  plot = alpha_panel_july,
  width = 14, height = 10, dpi = 600, compression = "lzw"
)

ggsave(
  file.path(outdir_alpha, "09_alpha_panel_October.png"),
  plot = alpha_panel_october,
  width = 14, height = 10, dpi = 600
)
ggsave(
  file.path(outdir_alpha, "09_alpha_panel_October.tiff"),
  plot = alpha_panel_october,
  width = 14, height = 10, dpi = 600, compression = "lzw"
)

############################################################
# 14) TABLE READY FOR WRITING
############################################################

alpha_table_for_writing <- alpha_summary_by_group %>%
  mutate(
    mean_se = paste0(round(mean, 2), " ± ", round(se, 2))
  ) %>%
  select(
    Metric,
    Month,
    Compartment,
    Health_status,
    n,
    mean,
    se,
    mean_se,
    median,
    q1,
    q3,
    min,
    max
  )

write_csv(
  alpha_table_for_writing,
  file.path(outdir_alpha, "10_alpha_table_for_writing.csv")
)

############################################################
# 15) PRINT KEY OUTPUTS TO CONSOLE
############################################################

cat("\n====================================================\n")
cat("ALPHA VALUES LONG\n")
cat("====================================================\n")
print(dplyr::slice_head(alpha_values_long, n = 20), width = Inf)

cat("\n====================================================\n")
cat("ALPHA SUMMARY BY GROUP\n")
cat("====================================================\n")
print(alpha_summary_by_group, n = Inf, width = Inf)

cat("\n====================================================\n")
cat("WITHIN-COMPARTMENT HEALTHY VS DISEASED TESTS\n")
cat("====================================================\n")
print(alpha_within_compartment_tests, n = Inf, width = Inf)

cat("\n====================================================\n")
cat("WITHIN-COMPARTMENT WILCOX VALIDATION\n")
cat("====================================================\n")
print(alpha_within_compartment_wilcox, n = Inf, width = Inf)

cat("\n====================================================\n")
cat("FILES WRITTEN TO:\n")
cat("====================================================\n")
print(normalizePath(outdir_alpha))



############################################################
# SECTION 2B — KEY ALPHA CONTRASTS TABLE
#
# Objetivos:
# 1) Unir los contrasts principales de alpha diversity
# 2) Integrar Tukey + Wilcoxon en una sola tabla
# 3) Crear una versión compacta para Results / Supplementary
############################################################

############################################################
# 16) MERGE TUKEY + WILCOX FOR KEY CONTRASTS
############################################################

alpha_key_contrasts <- alpha_within_compartment_tests %>%
  left_join(
    alpha_within_compartment_wilcox %>%
      select(Metric, Month, Compartment, wilcox_p),
    by = c("Metric", "Month", "Compartment")
  ) %>%
  mutate(
    direction = case_when(
      is.na(p_adj) ~ NA_character_,
      p_adj >= 0.05 ~ "No significant difference",
      estimate > 0 ~ "Diseased > Healthy",
      estimate < 0 ~ "Healthy > Diseased",
      TRUE ~ "No significant difference"
    ),
    tukey_sig = case_when(
      is.na(p_adj) ~ NA_character_,
      p_adj < 0.001 ~ "***",
      p_adj < 0.01  ~ "**",
      p_adj < 0.05  ~ "*",
      TRUE ~ "ns"
    ),
    wilcox_sig = case_when(
      is.na(wilcox_p) ~ NA_character_,
      wilcox_p < 0.001 ~ "***",
      wilcox_p < 0.01  ~ "**",
      wilcox_p < 0.05  ~ "*",
      TRUE ~ "ns"
    )
  ) %>%
  arrange(Month, Compartment, Metric)

write_csv(
  alpha_key_contrasts,
  file.path(outdir_alpha, "11_alpha_key_contrasts_tukey_wilcox.csv")
)

############################################################
# 17) PAPER-READY COMPACT TABLE
############################################################

alpha_key_contrasts_compact <- alpha_key_contrasts %>%
  mutate(
    estimate_round = round(estimate, digits_by_metric(Metric)),
    conf.low_round = round(conf.low, digits_by_metric(Metric)),
    conf.high_round = round(conf.high, digits_by_metric(Metric)),
    p_adj_round = signif(p_adj, 3),
    wilcox_p_round = signif(wilcox_p, 3),
    estimate_ci = paste0(
      estimate_round,
      " [",
      conf.low_round,
      ", ",
      conf.high_round,
      "]"
    )
  ) %>%
  select(
    Metric,
    Month,
    Compartment,
    contrast,
    direction,
    estimate_ci,
    p_adj_round,
    tukey_sig,
    wilcox_p_round,
    wilcox_sig
  ) %>%
  rename(
    Tukey_estimate_CI = estimate_ci,
    Tukey_p = p_adj_round,
    Wilcox_p = wilcox_p_round
  ) %>%
  arrange(Month, Compartment, Metric)

write_csv(
  alpha_key_contrasts_compact,
  file.path(outdir_alpha, "12_alpha_key_contrasts_compact_for_paper.csv")
)

############################################################
# 18) EVEN SHORTER WRITING TABLE
############################################################

alpha_key_contrasts_writing <- alpha_key_contrasts %>%
  mutate(
    result_summary = case_when(
      direction == "Diseased > Healthy" & !is.na(p_adj) & p_adj < 0.05 ~
        "Diseased significantly higher than healthy",
      direction == "Healthy > Diseased" & !is.na(p_adj) & p_adj < 0.05 ~
        "Healthy significantly higher than diseased",
      TRUE ~
        "No significant difference"
    )
  ) %>%
  select(
    Metric,
    Month,
    Compartment,
    result_summary,
    estimate,
    conf.low,
    conf.high,
    p_adj,
    wilcox_p
  ) %>%
  arrange(Month, Compartment, Metric)

write_csv(
  alpha_key_contrasts_writing,
  file.path(outdir_alpha, "13_alpha_key_contrasts_writing_table.csv")
)

############################################################
# 19) OPTIONAL WIDE TABLE FOR QUICK READING
############################################################

alpha_key_contrasts_wide <- alpha_key_contrasts_compact %>%
  unite("Comparison", Month, Compartment, Metric, sep = " | ") %>%
  pivot_wider(
    names_from = Comparison,
    values_from = c(
      direction,
      Tukey_estimate_CI,
      Tukey_p,
      tukey_sig,
      Wilcox_p,
      wilcox_sig
    )
  )

write_csv(
  alpha_key_contrasts_wide,
  file.path(outdir_alpha, "14_alpha_key_contrasts_wide.csv")
)

############################################################
# 20) PRINT TO CONSOLE
############################################################

cat("\n====================================================\n")
cat("ALPHA KEY CONTRASTS — FULL TABLE\n")
cat("====================================================\n")
print(alpha_key_contrasts, n = Inf, width = Inf)

cat("\n====================================================\n")
cat("ALPHA KEY CONTRASTS — COMPACT PAPER TABLE\n")
cat("====================================================\n")
print(alpha_key_contrasts_compact, n = Inf, width = Inf)

cat("\n====================================================\n")
cat("ALPHA KEY CONTRASTS — WRITING TABLE\n")
cat("====================================================\n")
print(alpha_key_contrasts_writing, n = Inf, width = Inf)







setwd("path/to/QIIME2_files_to_R_16S")

library(readr)
library(dplyr)

alpha16_master <- read_csv("01_alpha_diversity_16S/10_alpha_table_for_writing.csv", show_col_types = FALSE) %>%
  left_join(
    read_csv("01_alpha_diversity_16S/13_alpha_key_contrasts_writing_table.csv", show_col_types = FALSE),
    by = c("Metric", "Month", "Compartment")
  )

write_csv(alpha16_master, "01_alpha_diversity_16S/Table_S2_alpha_diversity_16S_master.csv")

print(alpha16_master, n = Inf, width = Inf)





###############################################################################################################################################################
###############################################################################################################################################################
###############################################################################################################################################################
###############################################################################################################################################################
###############################################################################################################################################################



############################################################
# SECTION 3 — BETA DIVERSITY (16S)
#
# Objetivos:
# 1) Calcular beta diversity global por mes
# 2) Calcular beta diversity dentro de cada compartimento por mes
# 3) Usar Bray-Curtis, Jaccard y UniFrac (si hay árbol)
# 4) Generar PCoAs y tablas estadísticas listas para paper
############################################################

suppressPackageStartupMessages({
  library(phyloseq)
  library(speedyseq)
  library(tidyverse)
  library(vegan)
  library(ape)
  library(patchwork)
  library(readr)
  library(ggplot2)
  library(reshape2)
})

############################################################
# 0) INPUTS AND OUTPUT FOLDERS
############################################################

physeq_path <- "Phyloseq_16S.rds"

outdir_beta <- "02_beta_diversity_16S"
outdir_beta_global <- file.path(outdir_beta, "global")
outdir_beta_by_comp <- file.path(outdir_beta, "by_compartment")
outdir_beta_tables <- file.path(outdir_beta, "tables")
outdir_beta_figs <- file.path(outdir_beta, "figures")

dir.create(outdir_beta, recursive = TRUE, showWarnings = FALSE)
dir.create(outdir_beta_global, recursive = TRUE, showWarnings = FALSE)
dir.create(outdir_beta_by_comp, recursive = TRUE, showWarnings = FALSE)
dir.create(outdir_beta_tables, recursive = TRUE, showWarnings = FALSE)
dir.create(outdir_beta_figs, recursive = TRUE, showWarnings = FALSE)

############################################################
# 1) LOAD PHYLOSEQ AND BASIC CHECKS
############################################################

ps0 <- readRDS(physeq_path)
ps0 <- prune_samples(sample_sums(ps0) > 0, ps0)

meta0 <- as(sample_data(ps0), "data.frame")
meta0$SampleID <- rownames(meta0)

required_cols <- c("Month", "Compartment", "Health_status")
missing_cols <- setdiff(required_cols, colnames(meta0))
if (length(missing_cols) > 0) {
  stop("Faltan columnas requeridas en sample_data: ", paste(missing_cols, collapse = ", "))
}

meta0$Month <- factor(as.character(meta0$Month), levels = c("July", "October"))
meta0$Compartment <- factor(
  as.character(meta0$Compartment),
  levels = c("Bulk soil", "Rhizosphere", "Endosphere", "Phyllosphere")
)
meta0$Health_status <- factor(
  as.character(meta0$Health_status),
  levels = c(
    "Bulk soil",
    "Healthy Rhizosphere", "Diseased Rhizosphere",
    "Healthy Endosphere",  "Diseased Endosphere",
    "Healthy Phyllosphere","Diseased Phyllosphere"
  )
)

sample_data(ps0) <- sample_data(meta0)

cat("\n====================================================\n")
cat("BETA DIVERSITY — BASIC CHECKS\n")
cat("====================================================\n")
cat("Working directory:\n")
print(getwd())

cat("\nNumber of samples in phyloseq:\n")
print(nsamples(ps0))

cat("\nNumber of taxa in phyloseq:\n")
print(ntaxa(ps0))

cat("\nSample-data columns:\n")
print(colnames(as(sample_data(ps0), "data.frame")))

############################################################
# 2) COLORS
############################################################

group_colors <- c(
  "Bulk soil"              = "#896054",
  "Healthy Rhizosphere"    = "#4B9EA7",
  "Diseased Rhizosphere"   = "#801723",
  "Healthy Endosphere"     = "#6B9863",
  "Diseased Endosphere"    = "#E7C313",
  "Healthy Phyllosphere"   = "#4A5263",
  "Diseased Phyllosphere"  = "#6B4EA7"
)

binary_colors <- c(
  "Healthy" = "#4B9EA7",
  "Diseased" = "#801723"
)

############################################################
# 3) HELPERS
############################################################

pcoa_axis_labels <- function(pcoa_obj) {
  vals <- pcoa_obj$values$Relative_eig
  c(
    paste0("PCoA1 (", round(100 * vals[1], 1), "%)"),
    paste0("PCoA2 (", round(100 * vals[2], 1), "%)")
  )
}

make_health_binary <- function(x) {
  dplyr::case_when(
    grepl("^Healthy", x) ~ "Healthy",
    grepl("^Diseased", x) ~ "Diseased",
    TRUE ~ as.character(x)
  )
}

prepare_metadata_for_join <- function(md) {
  if (!"SampleID" %in% colnames(md)) {
    md$SampleID <- rownames(md)
  }
  md
}

extract_distance_matrix <- function(ps, method = c("bray", "jaccard")) {
  method <- match.arg(method)
  
  otu_mat <- as(otu_table(ps), "matrix")
  if (taxa_are_rows(ps)) otu_mat <- t(otu_mat)
  
  if (method == "bray") {
    return(vegan::vegdist(otu_mat, method = "bray"))
  }
  
  if (method == "jaccard") {
    return(vegan::vegdist((otu_mat > 0) * 1, method = "jaccard", binary = TRUE))
  }
}

mean_intra_inter <- function(D, metadata, group_var) {
  M <- as.matrix(D)
  diag(M) <- NA_real_
  
  grp <- metadata[[group_var]]
  names(grp) <- rownames(metadata)
  
  long <- reshape2::melt(M, varnames = c("i", "j"), value.name = "dist", na.rm = TRUE)
  long$g1 <- grp[long$i]
  long$g2 <- grp[long$j]
  
  long <- long %>%
    dplyr::filter(!is.na(g1), !is.na(g2)) %>%
    dplyr::mutate(Comparison = ifelse(g1 == g2, "Intra", "Inter"))
  
  long %>%
    dplyr::group_by(Comparison) %>%
    dplyr::summarise(
      mean = mean(dist, na.rm = TRUE),
      sd = sd(dist, na.rm = TRUE),
      .groups = "drop"
    )
}

plot_pcoa_global <- function(pcoa_df, title_text, xlab_text, ylab_text) {
  ggplot(pcoa_df, aes(x = Axis.1, y = Axis.2, color = Health_status, shape = Compartment)) +
    stat_ellipse(
      aes(fill = Health_status),
      geom = "polygon",
      type = "norm",
      level = 0.80,
      alpha = 0.15,
      color = NA
    ) +
    stat_ellipse(type = "norm", level = 0.80, linewidth = 0.8, show.legend = FALSE) +
    geom_point(size = 3, alpha = 0.9) +
    scale_color_manual(values = group_colors, drop = FALSE) +
    scale_fill_manual(values = group_colors, drop = FALSE) +
    guides(fill = "none") +
    theme_bw() +
    theme(
      legend.position = "right",
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold")
    ) +
    labs(
      title = title_text,
      x = xlab_text,
      y = ylab_text,
      color = "Group",
      shape = "Compartment"
    )
}

plot_pcoa_binary <- function(pcoa_df, title_text, xlab_text, ylab_text) {
  ggplot(pcoa_df, aes(x = Axis.1, y = Axis.2, color = Health_binary)) +
    stat_ellipse(
      aes(fill = Health_binary),
      geom = "polygon",
      type = "norm",
      level = 0.80,
      alpha = 0.15,
      color = NA
    ) +
    stat_ellipse(type = "norm", level = 0.80, linewidth = 0.8, show.legend = FALSE) +
    geom_point(size = 3, alpha = 0.9) +
    scale_color_manual(values = binary_colors, drop = FALSE) +
    scale_fill_manual(values = binary_colors, drop = FALSE) +
    guides(fill = "none") +
    theme_bw() +
    theme(
      legend.position = "right",
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold")
    ) +
    labs(
      title = title_text,
      x = xlab_text,
      y = ylab_text,
      color = "Health status"
    )
}

############################################################
# 4) GLOBAL ANALYSIS BY MONTH
############################################################

run_global_beta_one_month <- function(ps_in, month_name) {
  md0 <- as(sample_data(ps_in), "data.frame")
  keep_ids <- rownames(md0)[as.character(md0$Month) == month_name]
  
  ps_sub <- prune_samples(keep_ids, ps_in)
  ps_sub <- prune_samples(sample_sums(ps_sub) > 0, ps_sub)
  
  md <- as(sample_data(ps_sub), "data.frame")
  md <- prepare_metadata_for_join(md)
  
  D_bray <- extract_distance_matrix(ps_sub, "bray")
  D_jaccard <- extract_distance_matrix(ps_sub, "jaccard")
  
  pcoa_bray <- ape::pcoa(D_bray)
  pcoa_jaccard <- ape::pcoa(D_jaccard)
  
  pcoa_bray_df <- as.data.frame(pcoa_bray$vectors[, 1:2, drop = FALSE])
  pcoa_bray_df$SampleID <- rownames(pcoa_bray_df)
  pcoa_bray_df <- dplyr::left_join(pcoa_bray_df, md, by = "SampleID")
  
  pcoa_jaccard_df <- as.data.frame(pcoa_jaccard$vectors[, 1:2, drop = FALSE])
  pcoa_jaccard_df$SampleID <- rownames(pcoa_jaccard_df)
  pcoa_jaccard_df <- dplyr::left_join(pcoa_jaccard_df, md, by = "SampleID")
  
  bray_labels <- pcoa_axis_labels(pcoa_bray)
  jac_labels <- pcoa_axis_labels(pcoa_jaccard)
  
  p_bray <- plot_pcoa_global(
    pcoa_bray_df,
    paste0("Global PCoA — Bray-Curtis (", month_name, ")"),
    bray_labels[1],
    bray_labels[2]
  )
  
  p_jaccard <- plot_pcoa_global(
    pcoa_jaccard_df,
    paste0("Global PCoA — Jaccard (", month_name, ")"),
    jac_labels[1],
    jac_labels[2]
  )
  
  ggsave(
    file.path(outdir_beta_global, paste0("global_PCoA_Bray_", month_name, ".png")),
    plot = p_bray, width = 8, height = 6, dpi = 600
  )
  ggsave(
    file.path(outdir_beta_global, paste0("global_PCoA_Jaccard_", month_name, ".png")),
    plot = p_jaccard, width = 8, height = 6, dpi = 600
  )
  
  per_bray <- vegan::adonis2(D_bray ~ Health_status, data = md, permutations = 999)
  ano_bray <- vegan::anosim(D_bray, md$Health_status, permutations = 999)
  bd_bray <- vegan::betadisper(D_bray, md$Health_status)
  bdpt_bray <- vegan::permutest(bd_bray, permutations = 999)
  
  per_jaccard <- vegan::adonis2(D_jaccard ~ Health_status, data = md, permutations = 999)
  ano_jaccard <- vegan::anosim(D_jaccard, md$Health_status, permutations = 999)
  bd_jaccard <- vegan::betadisper(D_jaccard, md$Health_status)
  bdpt_jaccard <- vegan::permutest(bd_jaccard, permutations = 999)
  
  summary_tbl <- tibble(
    Month = month_name,
    Metric = c("Bray-Curtis", "Jaccard"),
    n_samples = nsamples(ps_sub),
    n_taxa = ntaxa(ps_sub),
    PCoA_axis1_var = c(
      100 * pcoa_bray$values$Relative_eig[1],
      100 * pcoa_jaccard$values$Relative_eig[1]
    ),
    PCoA_axis2_var = c(
      100 * pcoa_bray$values$Relative_eig[2],
      100 * pcoa_jaccard$values$Relative_eig[2]
    ),
    PERMANOVA_R2 = c(as.numeric(per_bray$R2[1]), as.numeric(per_jaccard$R2[1])),
    PERMANOVA_F  = c(as.numeric(per_bray$F[1]),  as.numeric(per_jaccard$F[1])),
    PERMANOVA_p  = c(as.numeric(per_bray$`Pr(>F)`[1]), as.numeric(per_jaccard$`Pr(>F)`[1])),
    ANOSIM_R     = c(as.numeric(ano_bray$statistic), as.numeric(ano_jaccard$statistic)),
    ANOSIM_p     = c(as.numeric(ano_bray$signif),    as.numeric(ano_jaccard$signif)),
    betadisper_F = c(as.numeric(bdpt_bray$tab$F[1]), as.numeric(bdpt_jaccard$tab$F[1])),
    betadisper_p = c(as.numeric(bdpt_bray$tab$`Pr(>F)`[1]), as.numeric(bdpt_jaccard$tab$`Pr(>F)`[1]))
  )
  
  write_csv(
    summary_tbl,
    file.path(outdir_beta_global, paste0("summary_global_", month_name, ".csv"))
  )
  
  list(
    summary = summary_tbl,
    p_bray = p_bray,
    p_jaccard = p_jaccard
  )
}

global_july <- run_global_beta_one_month(ps0, "July")
global_october <- run_global_beta_one_month(ps0, "October")

global_beta_summary <- dplyr::bind_rows(global_july$summary, global_october$summary)

write_csv(
  global_beta_summary,
  file.path(outdir_beta_tables, "01_beta_global_summary_July_October.csv")
)

############################################################
# 5) WITHIN-COMPARTMENT ANALYSIS BY MONTH
############################################################

run_compartment_beta <- function(ps_in, month_name, compartment_name) {
  md0 <- as(sample_data(ps_in), "data.frame")
  keep_ids <- rownames(md0)[
    as.character(md0$Month) == month_name &
      as.character(md0$Compartment) == compartment_name
  ]
  
  ps_sub <- prune_samples(keep_ids, ps_in)
  ps_sub <- prune_samples(sample_sums(ps_sub) > 0, ps_sub)
  
  md <- as(sample_data(ps_sub), "data.frame")
  md <- prepare_metadata_for_join(md)
  md$Health_binary <- factor(make_health_binary(md$Health_status), levels = c("Healthy", "Diseased"))
  sample_data(ps_sub) <- sample_data(md)
  
  D_bray <- extract_distance_matrix(ps_sub, "bray")
  D_jaccard <- extract_distance_matrix(ps_sub, "jaccard")
  
  pcoa_bray <- ape::pcoa(D_bray)
  pcoa_jaccard <- ape::pcoa(D_jaccard)
  
  pcoa_bray_df <- as.data.frame(pcoa_bray$vectors[, 1:2, drop = FALSE])
  pcoa_bray_df$SampleID <- rownames(pcoa_bray_df)
  pcoa_bray_df <- dplyr::left_join(pcoa_bray_df, md, by = "SampleID")
  
  pcoa_jaccard_df <- as.data.frame(pcoa_jaccard$vectors[, 1:2, drop = FALSE])
  pcoa_jaccard_df$SampleID <- rownames(pcoa_jaccard_df)
  pcoa_jaccard_df <- dplyr::left_join(pcoa_jaccard_df, md, by = "SampleID")
  
  bray_labels <- pcoa_axis_labels(pcoa_bray)
  jac_labels <- pcoa_axis_labels(pcoa_jaccard)
  
  p_bray <- plot_pcoa_binary(
    pcoa_bray_df,
    paste0(compartment_name, " — PCoA Bray-Curtis (", month_name, ")"),
    bray_labels[1],
    bray_labels[2]
  )
  
  p_jaccard <- plot_pcoa_binary(
    pcoa_jaccard_df,
    paste0(compartment_name, " — PCoA Jaccard (", month_name, ")"),
    jac_labels[1],
    jac_labels[2]
  )
  
  fig_prefix <- paste0(month_name, "_", gsub(" ", "_", compartment_name))
  
  ggsave(
    file.path(outdir_beta_by_comp, paste0(fig_prefix, "_PCoA_Bray.png")),
    plot = p_bray, width = 7, height = 6, dpi = 600
  )
  ggsave(
    file.path(outdir_beta_by_comp, paste0(fig_prefix, "_PCoA_Jaccard.png")),
    plot = p_jaccard, width = 7, height = 6, dpi = 600
  )
  
  per_bray <- vegan::adonis2(D_bray ~ Health_binary, data = md, permutations = 999)
  ano_bray <- vegan::anosim(D_bray, md$Health_binary, permutations = 999)
  bd_bray <- vegan::betadisper(D_bray, md$Health_binary)
  bdpt_bray <- vegan::permutest(bd_bray, permutations = 999)
  
  per_jaccard <- vegan::adonis2(D_jaccard ~ Health_binary, data = md, permutations = 999)
  ano_jaccard <- vegan::anosim(D_jaccard, md$Health_binary, permutations = 999)
  bd_jaccard <- vegan::betadisper(D_jaccard, md$Health_binary)
  bdpt_jaccard <- vegan::permutest(bd_jaccard, permutations = 999)
  
  intra_inter_bray <- mean_intra_inter(D_bray, md, "Health_binary") %>%
    tidyr::pivot_wider(names_from = Comparison, values_from = c(mean, sd))
  
  intra_inter_jaccard <- mean_intra_inter(D_jaccard, md, "Health_binary") %>%
    tidyr::pivot_wider(names_from = Comparison, values_from = c(mean, sd))
  
  summary_tbl <- tibble(
    Month = month_name,
    Compartment = compartment_name,
    n_samples = nsamples(ps_sub),
    n_Healthy = sum(md$Health_binary == "Healthy"),
    n_Diseased = sum(md$Health_binary == "Diseased"),
    
    Bray_PCoA_axis1_var = 100 * pcoa_bray$values$Relative_eig[1],
    Bray_PCoA_axis2_var = 100 * pcoa_bray$values$Relative_eig[2],
    Bray_PERMANOVA_R2 = as.numeric(per_bray$R2[1]),
    Bray_PERMANOVA_F  = as.numeric(per_bray$F[1]),
    Bray_PERMANOVA_p  = as.numeric(per_bray$`Pr(>F)`[1]),
    Bray_ANOSIM_R     = as.numeric(ano_bray$statistic),
    Bray_ANOSIM_p     = as.numeric(ano_bray$signif),
    Bray_betadisper_F = as.numeric(bdpt_bray$tab$F[1]),
    Bray_betadisper_p = as.numeric(bdpt_bray$tab$`Pr(>F)`[1]),
    Bray_mean_Intra = intra_inter_bray$mean_Intra,
    Bray_sd_Intra   = intra_inter_bray$sd_Intra,
    Bray_mean_Inter = intra_inter_bray$mean_Inter,
    Bray_sd_Inter   = intra_inter_bray$sd_Inter,
    
    Jaccard_PCoA_axis1_var = 100 * pcoa_jaccard$values$Relative_eig[1],
    Jaccard_PCoA_axis2_var = 100 * pcoa_jaccard$values$Relative_eig[2],
    Jaccard_PERMANOVA_R2 = as.numeric(per_jaccard$R2[1]),
    Jaccard_PERMANOVA_F  = as.numeric(per_jaccard$F[1]),
    Jaccard_PERMANOVA_p  = as.numeric(per_jaccard$`Pr(>F)`[1]),
    Jaccard_ANOSIM_R     = as.numeric(ano_jaccard$statistic),
    Jaccard_ANOSIM_p     = as.numeric(ano_jaccard$signif),
    Jaccard_betadisper_F = as.numeric(bdpt_jaccard$tab$F[1]),
    Jaccard_betadisper_p = as.numeric(bdpt_jaccard$tab$`Pr(>F)`[1]),
    Jaccard_mean_Intra = intra_inter_jaccard$mean_Intra,
    Jaccard_sd_Intra   = intra_inter_jaccard$sd_Intra,
    Jaccard_mean_Inter = intra_inter_jaccard$mean_Inter,
    Jaccard_sd_Inter   = intra_inter_jaccard$sd_Inter
  )
  
  write_csv(
    summary_tbl,
    file.path(outdir_beta_by_comp, paste0(fig_prefix, "_summary.csv"))
  )
  
  summary_tbl
}

compartments_to_test <- c("Rhizosphere", "Endosphere", "Phyllosphere")
months_to_test <- c("July", "October")

beta_by_compartment_summary <- purrr::map_dfr(
  months_to_test,
  function(mo) {
    purrr::map_dfr(
      compartments_to_test,
      function(co) run_compartment_beta(ps0, mo, co)
    )
  }
)

write_csv(
  beta_by_compartment_summary,
  file.path(outdir_beta_tables, "02_beta_by_compartment_summary_July_October.csv")
)

############################################################
# 6) OPTIONAL UNIFRAC ANALYSIS
############################################################

has_tree <- !is.null(phyloseq::phy_tree(ps0, errorIfNULL = FALSE))

if (isTRUE(has_tree)) {
  run_unifrac_compartment <- function(ps_in, month_name, compartment_name) {
    md0 <- as(sample_data(ps_in), "data.frame")
    keep_ids <- rownames(md0)[
      as.character(md0$Month) == month_name &
        as.character(md0$Compartment) == compartment_name
    ]
    
    ps_sub <- prune_samples(keep_ids, ps_in)
    ps_sub <- prune_samples(sample_sums(ps_sub) > 0, ps_sub)
    ps_sub <- prune_taxa(taxa_sums(ps_sub) > 0, ps_sub)
    
    tr <- phyloseq::phy_tree(ps_sub)
    common <- intersect(taxa_names(ps_sub), tr$tip.label)
    ps_sub <- prune_taxa(common, ps_sub)
    
    md <- as(sample_data(ps_sub), "data.frame")
    md <- prepare_metadata_for_join(md)
    md$Health_binary <- factor(make_health_binary(md$Health_status), levels = c("Healthy", "Diseased"))
    sample_data(ps_sub) <- sample_data(md)
    
    D_wu  <- phyloseq::distance(ps_sub, method = "unifrac", weighted = TRUE, fast = TRUE)
    D_uwu <- phyloseq::distance(ps_sub, method = "unifrac", weighted = FALSE, fast = TRUE)
    
    per_wu  <- vegan::adonis2(D_wu  ~ Health_binary, data = md, permutations = 999)
    per_uwu <- vegan::adonis2(D_uwu ~ Health_binary, data = md, permutations = 999)
    
    ano_wu  <- vegan::anosim(D_wu,  md$Health_binary, permutations = 999)
    ano_uwu <- vegan::anosim(D_uwu, md$Health_binary, permutations = 999)
    
    bd_wu   <- vegan::betadisper(D_wu,  md$Health_binary)
    bd_uwu  <- vegan::betadisper(D_uwu, md$Health_binary)
    
    bdpt_wu  <- vegan::permutest(bd_wu, permutations = 999)
    bdpt_uwu <- vegan::permutest(bd_uwu, permutations = 999)
    
    pcoa_wu <- ape::pcoa(D_wu)
    pcoa_uwu <- ape::pcoa(D_uwu)
    
    tibble(
      Month = month_name,
      Compartment = compartment_name,
      Metric = c("Weighted UniFrac", "Unweighted UniFrac"),
      PCoA_axis1_var = c(
        100 * pcoa_wu$values$Relative_eig[1],
        100 * pcoa_uwu$values$Relative_eig[1]
      ),
      PCoA_axis2_var = c(
        100 * pcoa_wu$values$Relative_eig[2],
        100 * pcoa_uwu$values$Relative_eig[2]
      ),
      PERMANOVA_R2 = c(as.numeric(per_wu$R2[1]), as.numeric(per_uwu$R2[1])),
      PERMANOVA_F  = c(as.numeric(per_wu$F[1]),  as.numeric(per_uwu$F[1])),
      PERMANOVA_p  = c(as.numeric(per_wu$`Pr(>F)`[1]), as.numeric(per_uwu$`Pr(>F)`[1])),
      ANOSIM_R     = c(as.numeric(ano_wu$statistic), as.numeric(ano_uwu$statistic)),
      ANOSIM_p     = c(as.numeric(ano_wu$signif),    as.numeric(ano_uwu$signif)),
      betadisper_F = c(as.numeric(bdpt_wu$tab$F[1]), as.numeric(bdpt_uwu$tab$F[1])),
      betadisper_p = c(as.numeric(bdpt_wu$tab$`Pr(>F)`[1]), as.numeric(bdpt_uwu$tab$`Pr(>F)`[1]))
    )
  }
  
  beta_unifrac_summary <- purrr::map_dfr(
    months_to_test,
    function(mo) {
      purrr::map_dfr(
        compartments_to_test,
        function(co) run_unifrac_compartment(ps0, mo, co)
      )
    }
  )
  
  write_csv(
    beta_unifrac_summary,
    file.path(outdir_beta_tables, "03_beta_unifrac_by_compartment_July_October.csv")
  )
} else {
  beta_unifrac_summary <- tibble()
}

############################################################
# 7) COMPACT PAPER TABLES
############################################################

beta_key_table <- beta_by_compartment_summary %>%
  transmute(
    Month,
    Compartment,
    
    Bray_result = case_when(
      Bray_PERMANOVA_p < 0.05 & Bray_betadisper_p >= 0.05 ~ "Significant centroid shift",
      Bray_PERMANOVA_p < 0.05 & Bray_betadisper_p < 0.05  ~ "Significant, but dispersion differs",
      TRUE ~ "Not significant"
    ),
    Bray_R2 = round(Bray_PERMANOVA_R2, 3),
    Bray_F = round(Bray_PERMANOVA_F, 3),
    Bray_p = signif(Bray_PERMANOVA_p, 3),
    Bray_ANOSIM_R = round(Bray_ANOSIM_R, 3),
    Bray_betadisper_p = signif(Bray_betadisper_p, 3),
    
    Jaccard_result = case_when(
      Jaccard_PERMANOVA_p < 0.05 & Jaccard_betadisper_p >= 0.05 ~ "Significant centroid shift",
      Jaccard_PERMANOVA_p < 0.05 & Jaccard_betadisper_p < 0.05  ~ "Significant, but dispersion differs",
      TRUE ~ "Not significant"
    ),
    Jaccard_R2 = round(Jaccard_PERMANOVA_R2, 3),
    Jaccard_F = round(Jaccard_PERMANOVA_F, 3),
    Jaccard_p = signif(Jaccard_PERMANOVA_p, 3),
    Jaccard_ANOSIM_R = round(Jaccard_ANOSIM_R, 3),
    Jaccard_betadisper_p = signif(Jaccard_betadisper_p, 3)
  ) %>%
  arrange(Month, Compartment)

write_csv(
  beta_key_table,
  file.path(outdir_beta_tables, "04_beta_key_table_for_paper.csv")
)

############################################################
# 8) PRINT TO CONSOLE
############################################################

cat("\n====================================================\n")
cat("GLOBAL BETA SUMMARY\n")
cat("====================================================\n")
print(global_beta_summary, n = Inf, width = Inf)

cat("\n====================================================\n")
cat("BETA SUMMARY BY COMPARTMENT\n")
cat("====================================================\n")
print(beta_by_compartment_summary, n = Inf, width = Inf)

if (nrow(beta_unifrac_summary) > 0) {
  cat("\n====================================================\n")
  cat("UNIFRAC SUMMARY BY COMPARTMENT\n")
  cat("====================================================\n")
  print(beta_unifrac_summary, n = Inf, width = Inf)
}

cat("\n====================================================\n")
cat("BETA KEY TABLE FOR PAPER\n")
cat("====================================================\n")
print(beta_key_table, n = Inf, width = Inf)

cat("\n====================================================\n")
cat("FILES WRITTEN TO:\n")
cat("====================================================\n")
print(normalizePath(outdir_beta))


############################################################
# 9) SHOW EXPORTED PNG FIGURES ON SCREEN
############################################################

suppressPackageStartupMessages({
  library(png)
  library(grid)
})

show_pngs <- function(folder) {
  png_files <- list.files(folder, pattern = "\\.png$", full.names = TRUE)
  
  if (length(png_files) == 0) {
    cat("\nNo PNG files found in:", folder, "\n")
    return(invisible(NULL))
  }
  
  cat("\nShowing PNG files from:\n", folder, "\n")
  print(basename(png_files))
  
  for (f in png_files) {
    img <- png::readPNG(f)
    grid::grid.newpage()
    grid::grid.raster(img)
    grid::grid.text(
      basename(f),
      x = grid::unit(0.02, "npc"),
      y = grid::unit(0.98, "npc"),
      just = c("left", "top")
    )
  }
}

cat("\n====================================================\n")
cat("SHOWING GLOBAL BETA FIGURES\n")
cat("====================================================\n")
show_pngs(outdir_beta_global)

cat("\n====================================================\n")
cat("SHOWING BY-COMPARTMENT BETA FIGURES\n")
cat("====================================================\n")
show_pngs(outdir_beta_by_comp)

cat("\n====================================================\n")
cat("SHOWING FIGURES FOLDER\n")
cat("====================================================\n")
show_pngs(outdir_beta_figs)




setwd("path/to/QIIME2_files_to_R_16S")

library(readr)
library(dplyr)

beta16_global <- read_csv("02_beta_diversity_16S/tables/01_beta_global_summary_July_October.csv", show_col_types = FALSE)
beta16_comp   <- read_csv("02_beta_diversity_16S/tables/02_beta_by_compartment_summary_July_October.csv", show_col_types = FALSE)
beta16_unifrac <- read_csv("02_beta_diversity_16S/tables/03_beta_unifrac_by_compartment_July_October.csv", show_col_types = FALSE)
beta16_key <- read_csv("02_beta_diversity_16S/tables/04_beta_key_table_for_paper.csv", show_col_types = FALSE)

write_csv(beta16_global, "02_beta_diversity_16S/tables/Table_S4_beta_16S_global_master.csv")
write_csv(beta16_comp, "02_beta_diversity_16S/tables/Table_S5_beta_16S_compartment_master.csv")
write_csv(beta16_unifrac, "02_beta_diversity_16S/tables/Table_S6_beta_16S_unifrac_master.csv")

print(beta16_global, n = Inf, width = Inf)
print(beta16_key, n = Inf, width = Inf)
print(beta16_unifrac, n = Inf, width = Inf)










###############################################################################################################################################################
###############################################################################################################################################################
###############################################################################################################################################################
###############################################################################################################################################################
###############################################################################################################################################################


############################################################
# SECTION 4 — TAXONOMIC COMPOSITION (16S)
#
# Objetivos:
# 1) Calcular abundancia relativa taxonómica por muestra
# 2) Resumir composición por mes, compartimento y estado sanitario
# 3) Generar tablas de top phyla y top genera
# 4) Exportar figuras descriptivas listas para paper
############################################################

suppressPackageStartupMessages({
  library(phyloseq)
  library(speedyseq)
  library(tidyverse)
  library(ggplot2)
  library(patchwork)
  library(readr)
  library(stringr)
  library(forcats)
})

############################################################
# 0) INPUTS AND OUTPUT FOLDERS
############################################################

physeq_path <- "Phyloseq_16S.rds"

outdir_tax <- "03_taxonomy_16S"
outdir_tax_tables <- file.path(outdir_tax, "tables")
outdir_tax_figs <- file.path(outdir_tax, "figures")

dir.create(outdir_tax, recursive = TRUE, showWarnings = FALSE)
dir.create(outdir_tax_tables, recursive = TRUE, showWarnings = FALSE)
dir.create(outdir_tax_figs, recursive = TRUE, showWarnings = FALSE)

############################################################
# 1) LOAD PHYLOSEQ AND BASIC CHECKS
############################################################

ps_tax0 <- readRDS(physeq_path)
ps_tax0 <- prune_samples(sample_sums(ps_tax0) > 0, ps_tax0)

md0 <- as(sample_data(ps_tax0), "data.frame")
md0$SampleID <- rownames(md0)

required_cols <- c("Month", "Compartment", "Health_status")
missing_cols <- setdiff(required_cols, colnames(md0))
if (length(missing_cols) > 0) {
  stop("Faltan columnas requeridas en sample_data: ", paste(missing_cols, collapse = ", "))
}

tax0 <- as.data.frame(tax_table(ps_tax0))
required_tax_cols <- c("Phylum", "Genus")
missing_tax_cols <- setdiff(required_tax_cols, colnames(tax0))
if (length(missing_tax_cols) > 0) {
  stop("Faltan columnas requeridas en tax_table: ", paste(missing_tax_cols, collapse = ", "))
}

md0$Month <- factor(as.character(md0$Month), levels = c("July", "October"))
md0$Compartment <- factor(
  as.character(md0$Compartment),
  levels = c("Bulk soil", "Rhizosphere", "Endosphere", "Phyllosphere")
)
md0$Health_status <- factor(
  as.character(md0$Health_status),
  levels = c(
    "Bulk soil",
    "Healthy Rhizosphere", "Diseased Rhizosphere",
    "Healthy Endosphere",  "Diseased Endosphere",
    "Healthy Phyllosphere","Diseased Phyllosphere"
  )
)
sample_data(ps_tax0) <- sample_data(md0)

cat("\n====================================================\n")
cat("TAXONOMY — BASIC CHECKS\n")
cat("====================================================\n")
cat("Number of samples:\n")
print(nsamples(ps_tax0))

cat("\nNumber of taxa:\n")
print(ntaxa(ps_tax0))

cat("\nTaxonomy columns:\n")
print(colnames(tax0))

############################################################
# 2) HELPERS
############################################################

clean_taxon_name <- function(x, rank_name = "Taxon") {
  x <- as.character(x)
  x[is.na(x) | x == ""] <- paste0("Unclassified_", rank_name)
  x
}

make_relative_phyloseq <- function(ps) {
  transform_sample_counts(ps, function(x) {
    if (sum(x) == 0) x else x / sum(x)
  })
}

summarise_rel_abundance <- function(df_long) {
  df_long %>%
    group_by(Month, Compartment, Health_status, Taxon) %>%
    summarise(
      mean_rel_abundance = mean(Abundance, na.rm = TRUE),
      sd_rel_abundance = sd(Abundance, na.rm = TRUE),
      se_rel_abundance = sd_rel_abundance / sqrt(dplyr::n()),
      prevalence = sum(Abundance > 0, na.rm = TRUE),
      prevalence_pct = 100 * prevalence / dplyr::n(),
      .groups = "drop"
    )
}

get_top_taxa <- function(df_summary, n_top = 15) {
  df_summary %>%
    group_by(Taxon) %>%
    summarise(global_mean = mean(mean_rel_abundance, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(global_mean)) %>%
    slice_head(n = n_top) %>%
    pull(Taxon)
}

make_tax_barplot_grouped <- function(df_summary, top_taxa, title_text, ylab_text = "Mean relative abundance (%)") {
  df_plot <- df_summary %>%
    mutate(
      Taxon2 = ifelse(Taxon %in% top_taxa, Taxon, "Other")
    ) %>%
    group_by(Month, Health_status, Taxon2) %>%
    summarise(
      mean_rel_abundance = sum(mean_rel_abundance, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      Taxon2 = forcats::fct_reorder(Taxon2, mean_rel_abundance, .fun = mean, .desc = TRUE)
    )
  
  ggplot(df_plot, aes(x = Health_status, y = 100 * mean_rel_abundance, fill = Taxon2)) +
    geom_col(color = "black", linewidth = 0.2) +
    facet_wrap(~ Month, ncol = 1, scales = "free_x") +
    labs(
      title = title_text,
      x = NULL,
      y = ylab_text,
      fill = NULL
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold"),
      legend.position = "right"
    )
}

make_tax_barplot_compartment <- function(df_summary, top_taxa, compartment_name, title_text, ylab_text = "Mean relative abundance (%)") {
  df_plot <- df_summary %>%
    filter(Compartment == compartment_name) %>%
    mutate(
      Taxon2 = ifelse(Taxon %in% top_taxa, Taxon, "Other")
    ) %>%
    group_by(Month, Health_status, Taxon2) %>%
    summarise(
      mean_rel_abundance = sum(mean_rel_abundance, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      Taxon2 = forcats::fct_reorder(Taxon2, mean_rel_abundance, .fun = mean, .desc = TRUE)
    )
  
  ggplot(df_plot, aes(x = Health_status, y = 100 * mean_rel_abundance, fill = Taxon2)) +
    geom_col(color = "black", linewidth = 0.2) +
    facet_wrap(~ Month, ncol = 1, scales = "free_x") +
    labs(
      title = title_text,
      x = NULL,
      y = ylab_text,
      fill = NULL
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold"),
      legend.position = "right"
    )
}

############################################################
# 3) BUILD RELATIVE-ABUNDANCE PHYLOSEQ OBJECTS
############################################################

ps_rel <- make_relative_phyloseq(ps_tax0)

############################################################
# 4) PHYLUM-LEVEL TABLES
############################################################

ps_phylum <- tax_glom(ps_rel, taxrank = "Phylum", NArm = FALSE)

phylum_long <- psmelt(ps_phylum) %>%
  mutate(
    Taxon = clean_taxon_name(Phylum, "Phylum")
  ) %>%
  transmute(
    SampleID = Sample,
    Month,
    Compartment,
    Health_status,
    Taxon,
    Abundance
  ) %>%
  arrange(Month, Compartment, Health_status, SampleID, desc(Abundance))

write_csv(
  phylum_long,
  file.path(outdir_tax_tables, "01_phylum_relative_abundance_by_sample.csv")
)

phylum_summary <- summarise_rel_abundance(phylum_long) %>%
  arrange(Month, Compartment, Health_status, desc(mean_rel_abundance))

write_csv(
  phylum_summary,
  file.path(outdir_tax_tables, "02_phylum_relative_abundance_summary.csv")
)

top15_phyla <- get_top_taxa(phylum_summary, n_top = 15)

top15_phylum_summary <- phylum_summary %>%
  filter(Taxon %in% top15_phyla) %>%
  arrange(Month, Compartment, Health_status, desc(mean_rel_abundance))

write_csv(
  top15_phylum_summary,
  file.path(outdir_tax_tables, "03_top15_phyla_summary.csv")
)

############################################################
# 5) GENUS-LEVEL TABLES
############################################################

ps_genus <- tax_glom(ps_rel, taxrank = "Genus", NArm = FALSE)

genus_long <- psmelt(ps_genus) %>%
  mutate(
    Taxon = clean_taxon_name(Genus, "Genus")
  ) %>%
  transmute(
    SampleID = Sample,
    Month,
    Compartment,
    Health_status,
    Taxon,
    Abundance
  ) %>%
  arrange(Month, Compartment, Health_status, SampleID, desc(Abundance))

write_csv(
  genus_long,
  file.path(outdir_tax_tables, "04_genus_relative_abundance_by_sample.csv")
)

genus_summary <- summarise_rel_abundance(genus_long) %>%
  arrange(Month, Compartment, Health_status, desc(mean_rel_abundance))

write_csv(
  genus_summary,
  file.path(outdir_tax_tables, "05_genus_relative_abundance_summary.csv")
)

top20_genera <- get_top_taxa(genus_summary, n_top = 20)

top20_genus_summary <- genus_summary %>%
  filter(Taxon %in% top20_genera) %>%
  arrange(Month, Compartment, Health_status, desc(mean_rel_abundance))

write_csv(
  top20_genus_summary,
  file.path(outdir_tax_tables, "06_top20_genera_summary.csv")
)

############################################################
# 6) TAXONOMY TABLES FOR WRITING
############################################################

dominant_phylum_by_group <- phylum_summary %>%
  group_by(Month, Compartment, Health_status) %>%
  slice_max(order_by = mean_rel_abundance, n = 5, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(mean_rel_abundance_pct = 100 * mean_rel_abundance) %>%
  arrange(Month, Compartment, Health_status, desc(mean_rel_abundance_pct))

write_csv(
  dominant_phylum_by_group,
  file.path(outdir_tax_tables, "07_dominant_phyla_by_group_top5.csv")
)

dominant_genus_by_group <- genus_summary %>%
  group_by(Month, Compartment, Health_status) %>%
  slice_max(order_by = mean_rel_abundance, n = 10, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(mean_rel_abundance_pct = 100 * mean_rel_abundance) %>%
  arrange(Month, Compartment, Health_status, desc(mean_rel_abundance_pct))

write_csv(
  dominant_genus_by_group,
  file.path(outdir_tax_tables, "08_dominant_genera_by_group_top10.csv")
)

############################################################
# 7) PAPER-READY COMPACT TABLES
############################################################

phylum_compact_table <- top15_phylum_summary %>%
  mutate(
    mean_rel_abundance_pct = round(100 * mean_rel_abundance, 2),
    prevalence_pct = round(prevalence_pct, 1)
  ) %>%
  select(
    Month,
    Compartment,
    Health_status,
    Taxon,
    mean_rel_abundance_pct,
    prevalence_pct
  ) %>%
  arrange(Month, Compartment, Health_status, desc(mean_rel_abundance_pct))

write_csv(
  phylum_compact_table,
  file.path(outdir_tax_tables, "09_phylum_compact_table_for_paper.csv")
)

genus_compact_table <- top20_genus_summary %>%
  mutate(
    mean_rel_abundance_pct = round(100 * mean_rel_abundance, 2),
    prevalence_pct = round(prevalence_pct, 1)
  ) %>%
  select(
    Month,
    Compartment,
    Health_status,
    Taxon,
    mean_rel_abundance_pct,
    prevalence_pct
  ) %>%
  arrange(Month, Compartment, Health_status, desc(mean_rel_abundance_pct))

write_csv(
  genus_compact_table,
  file.path(outdir_tax_tables, "10_genus_compact_table_for_paper.csv")
)

############################################################
# 8) FIGURES — GLOBAL TAXONOMIC COMPOSITION
############################################################

p_phylum_global <- make_tax_barplot_grouped(
  phylum_summary,
  top_taxa = top15_phyla,
  title_text = "16S taxonomic composition at phylum level",
  ylab_text = "Mean relative abundance (%)"
)

p_genus_global <- make_tax_barplot_grouped(
  genus_summary,
  top_taxa = top20_genera,
  title_text = "16S taxonomic composition at genus level",
  ylab_text = "Mean relative abundance (%)"
)

ggsave(
  file.path(outdir_tax_figs, "01_phylum_global_barplot.png"),
  plot = p_phylum_global,
  width = 12, height = 8, dpi = 600
)

ggsave(
  file.path(outdir_tax_figs, "02_genus_global_barplot.png"),
  plot = p_genus_global,
  width = 13, height = 9, dpi = 600
)

############################################################
# 9) FIGURES — BY COMPARTMENT
############################################################

compartments_to_plot <- c("Rhizosphere", "Endosphere", "Phyllosphere")

for (comp_name in compartments_to_plot) {
  p_phylum_comp <- make_tax_barplot_compartment(
    phylum_summary,
    top_taxa = top15_phyla,
    compartment_name = comp_name,
    title_text = paste0("16S phylum-level composition — ", comp_name),
    ylab_text = "Mean relative abundance (%)"
  )
  
  p_genus_comp <- make_tax_barplot_compartment(
    genus_summary,
    top_taxa = top20_genera,
    compartment_name = comp_name,
    title_text = paste0("16S genus-level composition — ", comp_name),
    ylab_text = "Mean relative abundance (%)"
  )
  
  comp_file <- gsub(" ", "_", comp_name)
  
  ggsave(
    file.path(outdir_tax_figs, paste0("03_", comp_file, "_phylum_barplot.png")),
    plot = p_phylum_comp,
    width = 12, height = 8, dpi = 600
  )
  
  ggsave(
    file.path(outdir_tax_figs, paste0("04_", comp_file, "_genus_barplot.png")),
    plot = p_genus_comp,
    width = 13, height = 9, dpi = 600
  )
}

############################################################
# 10) OPTIONAL: TOP-TAXA HEATMAP TABLES
############################################################

phylum_heatmap_table <- phylum_summary %>%
  filter(Taxon %in% top15_phyla) %>%
  mutate(mean_rel_abundance_pct = 100 * mean_rel_abundance) %>%
  select(Month, Compartment, Health_status, Taxon, mean_rel_abundance_pct)

write_csv(
  phylum_heatmap_table,
  file.path(outdir_tax_tables, "11_phylum_heatmap_table.csv")
)

genus_heatmap_table <- genus_summary %>%
  filter(Taxon %in% top20_genera) %>%
  mutate(mean_rel_abundance_pct = 100 * mean_rel_abundance) %>%
  select(Month, Compartment, Health_status, Taxon, mean_rel_abundance_pct)

write_csv(
  genus_heatmap_table,
  file.path(outdir_tax_tables, "12_genus_heatmap_table.csv")
)

############################################################
# 11) PRINT TO CONSOLE
############################################################

cat("\n====================================================\n")
cat("TOP 15 PHYLA\n")
cat("====================================================\n")
print(top15_phyla)

cat("\n====================================================\n")
cat("TOP 20 GENERA\n")
cat("====================================================\n")
print(top20_genera)

cat("\n====================================================\n")
cat("DOMINANT PHYLA BY GROUP\n")
cat("====================================================\n")
print(dominant_phylum_by_group, n = Inf, width = Inf)

cat("\n====================================================\n")
cat("DOMINANT GENERA BY GROUP\n")
cat("====================================================\n")
print(dominant_genus_by_group, n = 50, width = Inf)

cat("\n====================================================\n")
cat("FILES WRITTEN TO:\n")
cat("====================================================\n")
print(normalizePath(outdir_tax))




setwd("path/to/QIIME2_files_to_R_16S")

library(readr)
library(dplyr)

phylum16 <- read_csv("03_taxonomy_16S/tables/09_phylum_compact_table_for_paper.csv", show_col_types = FALSE)
genus16  <- read_csv("03_taxonomy_16S/tables/10_genus_compact_table_for_paper.csv", show_col_types = FALSE)

cat("\n=== 16S TOP PHYLA BY GROUP ===\n")
phylum16 %>%
  group_by(Month, Compartment, Health_status) %>%
  slice_max(mean_rel_abundance_pct, n = 5, with_ties = FALSE) %>%
  arrange(Month, Compartment, Health_status, desc(mean_rel_abundance_pct)) %>%
  print(n = Inf, width = Inf)

cat("\n=== 16S TOP GENERA BY GROUP ===\n")
genus16 %>%
  group_by(Month, Compartment, Health_status) %>%
  slice_max(mean_rel_abundance_pct, n = 8, with_ties = FALSE) %>%
  arrange(Month, Compartment, Health_status, desc(mean_rel_abundance_pct)) %>%
  print(n = Inf, width = Inf)











###############################################################################################################################################################
###############################################################################################################################################################
###############################################################################################################################################################
###############################################################################################################################################################
###############################################################################################################################################################


############################################################
# SECTION 5 — DIFFERENTIAL ABUNDANCE (16S) WITH DESEQ2
#
# Comparaciones:
# - Healthy vs Diseased dentro de cada compartimento y mes
# - July / October
# - Rhizosphere / Endosphere / Phyllosphere
#
# Niveles:
# - ASV
# - Genus
#
# Notas:
# - A nivel ASV se conserva el hash real del ASV
# - A nivel Genus se agregan conteos manualmente con rowsum()
#   para evitar problemas de duplicados
############################################################

suppressPackageStartupMessages({
  library(phyloseq)
  library(DESeq2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(stringr)
  library(purrr)
})

############################################################
# 0) INPUTS AND OUTPUT FOLDERS
############################################################

physeq_path <- "Phyloseq_16S.rds"

outdir_deseq <- "04_differential_abundance_16S_DESeq2"
outdir_deseq_asv <- file.path(outdir_deseq, "ASV")
outdir_deseq_genus <- file.path(outdir_deseq, "Genus")
outdir_deseq_tables <- file.path(outdir_deseq, "tables")

dir.create(outdir_deseq, recursive = TRUE, showWarnings = FALSE)
dir.create(outdir_deseq_asv, recursive = TRUE, showWarnings = FALSE)
dir.create(outdir_deseq_genus, recursive = TRUE, showWarnings = FALSE)
dir.create(outdir_deseq_tables, recursive = TRUE, showWarnings = FALSE)

############################################################
# 1) SETTINGS
############################################################

comparisons_tbl <- tibble::tribble(
  ~Month,     ~Compartment,    ~Group_Healthy,             ~Group_Diseased,
  "July",     "Rhizosphere",   "Healthy Rhizosphere",      "Diseased Rhizosphere",
  "October",  "Rhizosphere",   "Healthy Rhizosphere",      "Diseased Rhizosphere",
  "July",     "Endosphere",    "Healthy Endosphere",       "Diseased Endosphere",
  "October",  "Endosphere",    "Healthy Endosphere",       "Diseased Endosphere",
  "July",     "Phyllosphere",  "Healthy Phyllosphere",     "Diseased Phyllosphere",
  "October",  "Phyllosphere",  "Healthy Phyllosphere",     "Diseased Phyllosphere"
)

DESEQ_ALPHA <- 0.05
DESEQ_FIT_TYPE <- "parametric"

# Prefiltros
MIN_TOTAL_COUNT_ASV <- 10
MIN_TOTAL_COUNT_GENUS <- 10

# Género: conservar solo ASVs con género asignado
KEEP_ONLY_ASSIGNED_GENUS <- TRUE

############################################################
# 2) LOAD PHYLOSEQ
############################################################

ps0 <- readRDS(physeq_path)
ps0 <- prune_samples(sample_sums(ps0) > 0, ps0)

md0 <- as(sample_data(ps0), "data.frame")
md0$SampleID <- rownames(md0)

required_cols <- c("Month", "Compartment", "Health_status")
missing_cols <- setdiff(required_cols, colnames(md0))
if (length(missing_cols) > 0) {
  stop("Missing required columns in sample_data: ", paste(missing_cols, collapse = ", "))
}

############################################################
# 3) HELPERS
############################################################

clean_genus_label <- function(x) {
  x <- as.character(x)
  x[is.na(x) | x == ""] <- NA_character_
  x
}

collapse_unique <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & x != ""]
  x <- unique(x)
  if (length(x) == 0) {
    return(NA_character_)
  }
  paste(x, collapse = " | ")
}

prepare_subset_for_comparison <- function(ps_in, month_name, compartment_name, healthy_label, diseased_label) {
  md <- as(sample_data(ps_in), "data.frame")
  
  keep_ids <- rownames(md)[
    as.character(md$Month) == month_name &
      as.character(md$Compartment) == compartment_name &
      as.character(md$Health_status) %in% c(healthy_label, diseased_label)
  ]
  
  ps_sub <- prune_samples(keep_ids, ps_in)
  ps_sub <- prune_samples(sample_sums(ps_sub) > 0, ps_sub)
  
  md_sub <- as(sample_data(ps_sub), "data.frame")
  md_sub$Condition <- ifelse(
    as.character(md_sub$Health_status) == healthy_label,
    "Healthy",
    "Diseased"
  )
  md_sub$Condition <- factor(md_sub$Condition, levels = c("Healthy", "Diseased"))
  sample_data(ps_sub) <- sample_data(md_sub)
  
  ps_sub
}

run_deseq_asv_one <- function(ps_sub, alpha = 0.05, min_total_count = 10, fit_type = "parametric") {
  count_mat <- as(otu_table(ps_sub), "matrix")
  if (!taxa_are_rows(ps_sub)) {
    count_mat <- t(count_mat)
  }
  
  tax_df <- as.data.frame(tax_table(ps_sub)) %>%
    rownames_to_column("ASV_hash")
  
  keep_asvs <- rowSums(count_mat) >= min_total_count
  count_mat_filt <- count_mat[keep_asvs, , drop = FALSE]
  tax_df_filt <- tax_df %>% filter(ASV_hash %in% rownames(count_mat_filt))
  
  coldata <- as(sample_data(ps_sub), "data.frame")
  coldata <- coldata[colnames(count_mat_filt), , drop = FALSE]
  
  dds <- DESeqDataSetFromMatrix(
    countData = round(count_mat_filt),
    colData = coldata,
    design = ~ Condition
  )
  
  dds <- dds[rowSums(counts(dds)) > 0, ]
  dds <- DESeq(dds, fitType = fit_type)
  
  res <- results(
    dds,
    contrast = c("Condition", "Diseased", "Healthy"),
    alpha = alpha
  )
  
  res_df <- as.data.frame(res) %>%
    rownames_to_column("ASV_hash") %>%
    as_tibble() %>%
    left_join(tax_df_filt, by = "ASV_hash") %>%
    mutate(
      Direction = case_when(
        !is.na(padj) & padj < alpha & log2FoldChange > 0 ~ "Diseased > Healthy",
        !is.na(padj) & padj < alpha & log2FoldChange < 0 ~ "Healthy > Diseased",
        TRUE ~ "No significant difference"
      )
    ) %>%
    arrange(padj, desc(abs(log2FoldChange)))
  
  sig_res_df <- res_df %>%
    filter(!is.na(padj), padj < alpha) %>%
    arrange(padj, desc(abs(log2FoldChange)))
  
  list(
    full = res_df,
    sig = sig_res_df,
    n_tested = nrow(res_df)
  )
}

run_deseq_genus_one <- function(ps_sub,
                                alpha = 0.05,
                                min_total_count = 10,
                                fit_type = "parametric",
                                keep_only_assigned_genus = TRUE) {
  count_mat <- as(otu_table(ps_sub), "matrix")
  if (!taxa_are_rows(ps_sub)) {
    count_mat <- t(count_mat)
  }
  
  tax_df <- as.data.frame(tax_table(ps_sub)) %>%
    rownames_to_column("ASV_hash")
  
  if (!"Genus" %in% colnames(tax_df)) {
    stop("Taxonomy table has no Genus column.")
  }
  
  tax_df$Genus <- clean_genus_label(tax_df$Genus)
  
  # Alinear tax_df con count_mat
  tax_df <- tax_df[match(rownames(count_mat), tax_df$ASV_hash), , drop = FALSE]
  
  if (!identical(rownames(count_mat), tax_df$ASV_hash)) {
    stop("ASV hash alignment failed between count matrix and taxonomy table.")
  }
  
  if (keep_only_assigned_genus) {
    keep_rows <- !is.na(tax_df$Genus)
  } else {
    tax_df$Genus[is.na(tax_df$Genus)] <- "Unclassified_Genus"
    keep_rows <- rep(TRUE, nrow(tax_df))
  }
  
  count_mat_genus_input <- count_mat[keep_rows, , drop = FALSE]
  tax_df_genus_input <- tax_df[keep_rows, , drop = FALSE]
  
  # Agregación segura por género
  genus_count_mat <- rowsum(
    x = count_mat_genus_input,
    group = tax_df_genus_input$Genus,
    reorder = FALSE
  )
  
  # Tabla taxonómica por género
  genus_lookup <- tax_df_genus_input %>%
    group_by(Genus) %>%
    summarise(
      Kingdom = collapse_unique(Kingdom),
      Phylum  = collapse_unique(Phylum),
      Class   = collapse_unique(Class),
      Order   = collapse_unique(Order),
      Family  = collapse_unique(Family),
      Species = collapse_unique(Species),
      .groups = "drop"
    )
  
  keep_genera <- rowSums(genus_count_mat) >= min_total_count
  genus_count_mat_filt <- genus_count_mat[keep_genera, , drop = FALSE]
  genus_lookup_filt <- genus_lookup %>%
    filter(Genus %in% rownames(genus_count_mat_filt))
  
  coldata <- as(sample_data(ps_sub), "data.frame")
  coldata <- coldata[colnames(genus_count_mat_filt), , drop = FALSE]
  
  dds <- DESeqDataSetFromMatrix(
    countData = round(genus_count_mat_filt),
    colData = coldata,
    design = ~ Condition
  )
  
  dds <- dds[rowSums(counts(dds)) > 0, ]
  dds <- DESeq(dds, fitType = fit_type)
  
  res <- results(
    dds,
    contrast = c("Condition", "Diseased", "Healthy"),
    alpha = alpha
  )
  
  res_df <- as.data.frame(res) %>%
    rownames_to_column("Genus") %>%
    as_tibble() %>%
    left_join(genus_lookup_filt, by = "Genus") %>%
    mutate(
      Direction = case_when(
        !is.na(padj) & padj < alpha & log2FoldChange > 0 ~ "Diseased > Healthy",
        !is.na(padj) & padj < alpha & log2FoldChange < 0 ~ "Healthy > Diseased",
        TRUE ~ "No significant difference"
      )
    ) %>%
    arrange(padj, desc(abs(log2FoldChange)))
  
  sig_res_df <- res_df %>%
    filter(!is.na(padj), padj < alpha) %>%
    arrange(padj, desc(abs(log2FoldChange)))
  
  list(
    full = res_df,
    sig = sig_res_df,
    n_tested = nrow(res_df)
  )
}

############################################################
# 4) RUN ALL COMPARISONS — ASV LEVEL
############################################################

all_asv_full <- list()
all_asv_sig <- list()

for (i in seq_len(nrow(comparisons_tbl))) {
  row_i <- comparisons_tbl[i, ]
  
  month_i <- row_i$Month
  comp_i <- row_i$Compartment
  healthy_i <- row_i$Group_Healthy
  diseased_i <- row_i$Group_Diseased
  
  comp_tag <- paste(month_i, comp_i, sep = "_") %>% gsub(" ", "_", .)
  
  cat("\n====================================================\n")
  cat("RUNNING DESEQ2 — ASV LEVEL:", comp_tag, "\n")
  cat("====================================================\n")
  
  ps_sub <- prepare_subset_for_comparison(
    ps_in = ps0,
    month_name = month_i,
    compartment_name = comp_i,
    healthy_label = healthy_i,
    diseased_label = diseased_i
  )
  
  res_list <- run_deseq_asv_one(
    ps_sub = ps_sub,
    alpha = DESEQ_ALPHA,
    min_total_count = MIN_TOTAL_COUNT_ASV,
    fit_type = DESEQ_FIT_TYPE
  )
  
  res_full <- res_list$full %>%
    mutate(
      Taxonomic_level = "ASV",
      Month = month_i,
      Compartment = comp_i,
      Comparison = paste0(healthy_i, " vs ", diseased_i),
      n_tested = res_list$n_tested
    ) %>%
    relocate(Taxonomic_level, Month, Compartment, Comparison, n_tested)
  
  res_sig <- res_list$sig %>%
    mutate(
      Taxonomic_level = "ASV",
      Month = month_i,
      Compartment = comp_i,
      Comparison = paste0(healthy_i, " vs ", diseased_i),
      n_tested = res_list$n_tested
    ) %>%
    relocate(Taxonomic_level, Month, Compartment, Comparison, n_tested)
  
  all_asv_full[[comp_tag]] <- res_full
  all_asv_sig[[comp_tag]] <- res_sig
  
  write_csv(
    res_full,
    file.path(outdir_deseq_asv, paste0("DESeq2_ASV_", comp_tag, "_full.csv"))
  )
  
  write_csv(
    res_sig,
    file.path(outdir_deseq_asv, paste0("DESeq2_ASV_", comp_tag, "_significant.csv"))
  )
}

asv_results_all <- bind_rows(all_asv_full)
asv_results_sig_all <- bind_rows(all_asv_sig)

write_csv(
  asv_results_all,
  file.path(outdir_deseq_tables, "01_DESeq2_ASV_all_comparisons_full.csv")
)

write_csv(
  asv_results_sig_all,
  file.path(outdir_deseq_tables, "02_DESeq2_ASV_all_comparisons_significant.csv")
)

############################################################
# 5) RUN ALL COMPARISONS — GENUS LEVEL
############################################################

all_genus_full <- list()
all_genus_sig <- list()

for (i in seq_len(nrow(comparisons_tbl))) {
  row_i <- comparisons_tbl[i, ]
  
  month_i <- row_i$Month
  comp_i <- row_i$Compartment
  healthy_i <- row_i$Group_Healthy
  diseased_i <- row_i$Group_Diseased
  
  comp_tag <- paste(month_i, comp_i, sep = "_") %>% gsub(" ", "_", .)
  
  cat("\n====================================================\n")
  cat("RUNNING DESEQ2 — GENUS LEVEL:", comp_tag, "\n")
  cat("====================================================\n")
  
  ps_sub <- prepare_subset_for_comparison(
    ps_in = ps0,
    month_name = month_i,
    compartment_name = comp_i,
    healthy_label = healthy_i,
    diseased_label = diseased_i
  )
  
  res_list <- run_deseq_genus_one(
    ps_sub = ps_sub,
    alpha = DESEQ_ALPHA,
    min_total_count = MIN_TOTAL_COUNT_GENUS,
    fit_type = DESEQ_FIT_TYPE,
    keep_only_assigned_genus = KEEP_ONLY_ASSIGNED_GENUS
  )
  
  res_full <- res_list$full %>%
    mutate(
      Taxonomic_level = "Genus",
      Month = month_i,
      Compartment = comp_i,
      Comparison = paste0(healthy_i, " vs ", diseased_i),
      n_tested = res_list$n_tested
    ) %>%
    relocate(Taxonomic_level, Month, Compartment, Comparison, n_tested)
  
  res_sig <- res_list$sig %>%
    mutate(
      Taxonomic_level = "Genus",
      Month = month_i,
      Compartment = comp_i,
      Comparison = paste0(healthy_i, " vs ", diseased_i),
      n_tested = res_list$n_tested
    ) %>%
    relocate(Taxonomic_level, Month, Compartment, Comparison, n_tested)
  
  all_genus_full[[comp_tag]] <- res_full
  all_genus_sig[[comp_tag]] <- res_sig
  
  write_csv(
    res_full,
    file.path(outdir_deseq_genus, paste0("DESeq2_Genus_", comp_tag, "_full.csv"))
  )
  
  write_csv(
    res_sig,
    file.path(outdir_deseq_genus, paste0("DESeq2_Genus_", comp_tag, "_significant.csv"))
  )
}

genus_results_all <- bind_rows(all_genus_full)
genus_results_sig_all <- bind_rows(all_genus_sig)

write_csv(
  genus_results_all,
  file.path(outdir_deseq_tables, "03_DESeq2_Genus_all_comparisons_full.csv")
)

write_csv(
  genus_results_sig_all,
  file.path(outdir_deseq_tables, "04_DESeq2_Genus_all_comparisons_significant.csv")
)

############################################################
# 6) COMPACT TABLES
############################################################

asv_results_compact <- asv_results_sig_all %>%
  mutate(
    log2FoldChange = round(log2FoldChange, 3),
    lfcSE = round(lfcSE, 3),
    stat = round(stat, 3),
    pvalue = signif(pvalue, 3),
    padj = signif(padj, 3)
  ) %>%
  select(
    Taxonomic_level,
    Month,
    Compartment,
    Comparison,
    n_tested,
    ASV_hash,
    Kingdom,
    Phylum,
    Class,
    Order,
    Family,
    Genus,
    Species,
    baseMean,
    log2FoldChange,
    lfcSE,
    stat,
    pvalue,
    padj,
    Direction
  ) %>%
  arrange(Month, Compartment, padj, desc(abs(log2FoldChange)))

write_csv(
  asv_results_compact,
  file.path(outdir_deseq_tables, "05_DESeq2_ASV_compact_for_paper.csv")
)

genus_results_compact <- genus_results_sig_all %>%
  mutate(
    log2FoldChange = round(log2FoldChange, 3),
    lfcSE = round(lfcSE, 3),
    stat = round(stat, 3),
    pvalue = signif(pvalue, 3),
    padj = signif(padj, 3)
  ) %>%
  select(
    Taxonomic_level,
    Month,
    Compartment,
    Comparison,
    n_tested,
    Genus,
    Kingdom,
    Phylum,
    Class,
    Order,
    Family,
    Species,
    baseMean,
    log2FoldChange,
    lfcSE,
    stat,
    pvalue,
    padj,
    Direction
  ) %>%
  arrange(Month, Compartment, padj, desc(abs(log2FoldChange)))

write_csv(
  genus_results_compact,
  file.path(outdir_deseq_tables, "06_DESeq2_Genus_compact_for_paper.csv")
)

############################################################
# 7) SUMMARY COUNTS
############################################################

asv_sig_counts <- asv_results_sig_all %>%
  group_by(Month, Compartment, Direction) %>%
  summarise(n_significant = n(), .groups = "drop") %>%
  arrange(Month, Compartment, Direction)

genus_sig_counts <- genus_results_sig_all %>%
  group_by(Month, Compartment, Direction) %>%
  summarise(n_significant = n(), .groups = "drop") %>%
  arrange(Month, Compartment, Direction)

write_csv(
  asv_sig_counts,
  file.path(outdir_deseq_tables, "07_DESeq2_ASV_significant_counts.csv")
)

write_csv(
  genus_sig_counts,
  file.path(outdir_deseq_tables, "08_DESeq2_Genus_significant_counts.csv")
)

############################################################
# 8) OPTIONAL WRITING TABLES
############################################################

asv_writing_table <- asv_results_sig_all %>%
  mutate(
    result_summary = case_when(
      Direction == "Diseased > Healthy" ~ "Diseased significantly enriched",
      Direction == "Healthy > Diseased" ~ "Healthy significantly enriched",
      TRUE ~ "No significant difference"
    )
  ) %>%
  select(
    Month,
    Compartment,
    ASV_hash,
    Genus,
    Phylum,
    log2FoldChange,
    padj,
    result_summary
  ) %>%
  arrange(Month, Compartment, padj, desc(abs(log2FoldChange)))

genus_writing_table <- genus_results_sig_all %>%
  mutate(
    result_summary = case_when(
      Direction == "Diseased > Healthy" ~ "Diseased significantly enriched",
      Direction == "Healthy > Diseased" ~ "Healthy significantly enriched",
      TRUE ~ "No significant difference"
    )
  ) %>%
  select(
    Month,
    Compartment,
    Genus,
    Phylum,
    log2FoldChange,
    padj,
    result_summary
  ) %>%
  arrange(Month, Compartment, padj, desc(abs(log2FoldChange)))

write_csv(
  asv_writing_table,
  file.path(outdir_deseq_tables, "09_DESeq2_ASV_writing_table.csv")
)

write_csv(
  genus_writing_table,
  file.path(outdir_deseq_tables, "10_DESeq2_Genus_writing_table.csv")
)

############################################################
# 9) PRINT TO CONSOLE
############################################################

cat("\n====================================================\n")
cat("ASV-LEVEL SIGNIFICANT COUNTS\n")
cat("====================================================\n")
print(asv_sig_counts, n = Inf, width = Inf)

cat("\n====================================================\n")
cat("GENUS-LEVEL SIGNIFICANT COUNTS\n")
cat("====================================================\n")
print(genus_sig_counts, n = Inf, width = Inf)

cat("\n====================================================\n")
cat("TOP ASV-LEVEL SIGNIFICANT RESULTS\n")
cat("====================================================\n")
print(dplyr::slice_head(asv_results_compact, n = 20), width = Inf)

cat("\n====================================================\n")
cat("TOP GENUS-LEVEL SIGNIFICANT RESULTS\n")
cat("====================================================\n")
print(dplyr::slice_head(genus_results_compact, n = 20), width = Inf)

cat("\n====================================================\n")
cat("FILES WRITTEN TO:\n")
cat("====================================================\n")
print(normalizePath(outdir_deseq))


############################################################
# SECTION 6 — FIGURES FOR DIFFERENTIAL ABUNDANCE (16S)
# DESEQ2 RESULTS — GENUS LEVEL
#
# Objetivos:
# 1) Graficar las 6 comparaciones Healthy vs Diseased
# 2) Mantener el estilo previo:
#    - Panel A: mean relative abundance (%)
#    - Panel B: log2FoldChange + error bars
# 3) Mostrar top 15 enriched in Healthy + top 15 enriched in Diseased
# 4) Excluir NA / uncultured / unclassified del gráfico principal
############################################################

suppressPackageStartupMessages({
  library(phyloseq)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(patchwork)
  library(ggnewscale)
  library(scales)
  library(stringr)
  library(tibble)
})

############################################################
# 0) INPUTS AND OUTPUT FOLDERS
############################################################

physeq_path <- "Phyloseq_16S.rds"

deseq_dir <- "04_differential_abundance_16S_DESeq2"
deseq_genus_dir <- file.path(deseq_dir, "Genus")

outdir_fig <- "07_DESeq2_figures_Genus"
dir.create(outdir_fig, recursive = TRUE, showWarnings = FALSE)

############################################################
# 1) SETTINGS
############################################################

comparisons_tbl <- tibble::tribble(
  ~Month,     ~Compartment,    ~Group_Healthy,             ~Group_Diseased,            ~short_label,
  "July",     "Rhizosphere",   "Healthy Rhizosphere",      "Diseased Rhizosphere",     "July_Rhizosphere",
  "October",  "Rhizosphere",   "Healthy Rhizosphere",      "Diseased Rhizosphere",     "October_Rhizosphere",
  "July",     "Endosphere",    "Healthy Endosphere",       "Diseased Endosphere",      "July_Endosphere",
  "October",  "Endosphere",    "Healthy Endosphere",       "Diseased Endosphere",      "October_Endosphere",
  "July",     "Phyllosphere",  "Healthy Phyllosphere",     "Diseased Phyllosphere",    "July_Phyllosphere",
  "October",  "Phyllosphere",  "Healthy Phyllosphere",     "Diseased Phyllosphere",    "October_Phyllosphere"
)

TOP_N_PER_DIRECTION <- 15

COL_HEALTHY  <- "#A6DBA0"
COL_DISEASED <- "#F4A6A6"
BAND_LIGHT   <- "#F2F2F2"
DOT_SIZE     <- 4.8
BASE_SIZE    <- 12

############################################################
# 2) HELPERS
############################################################

clean_genus_label <- function(x) {
  x <- as.character(x)
  x[is.na(x) | x == ""] <- NA_character_
  x
}

is_bad_genus <- function(x) {
  x <- as.character(x)
  is.na(x) |
    x == "" |
    str_detect(x, regex("^unclassified", ignore_case = TRUE)) |
    str_detect(x, regex("^uncultured$", ignore_case = TRUE)) |
    str_detect(x, regex("uncultured", ignore_case = TRUE))
}

p_label_fun <- function(x) {
  dplyr::case_when(
    is.na(x)  ~ "NA",
    x < 0.001 ~ "<0.001",
    x < 0.01  ~ "<0.01",
    x < 0.05  ~ "<0.05",
    TRUE      ~ paste0("=", formatC(x, format = "g", digits = 2))
  )
}

theme_big_labels <- function() {
  theme_minimal(base_size = BASE_SIZE) +
    theme(
      panel.background = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y  = element_text(size = 13, margin = margin(r = 6)),
      axis.title.x = element_text(size = 12),
      legend.position = "top",
      legend.title = element_blank(),
      plot.margin = unit(c(6, 6, 6, 24), "pt")
    )
}

############################################################
# 3) LOAD PHYLOSEQ
############################################################

ps0 <- readRDS(physeq_path)
ps0 <- prune_samples(sample_sums(ps0) > 0, ps0)

############################################################
# 4) FUNCTION TO BUILD ONE FIGURE
############################################################

make_deseq_genus_figure <- function(month_name,
                                    compartment_name,
                                    healthy_label,
                                    diseased_label,
                                    short_label,
                                    top_n_per_direction = 15) {
  
  # -------------------------------------------------------
  # 4.1 Read DESeq2 genus results for this comparison
  # -------------------------------------------------------
  res_file <- file.path(
    deseq_genus_dir,
    paste0("DESeq2_Genus_", short_label, "_significant.csv")
  )
  
  if (!file.exists(res_file)) {
    stop("Missing file: ", res_file)
  }
  
  res_df <- read_csv(res_file, show_col_types = FALSE)
  
  # -------------------------------------------------------
  # 4.2 Filter genera to plot
  # -------------------------------------------------------
  res_df <- res_df %>%
    mutate(
      Genus = clean_genus_label(Genus),
      Direction = as.character(Direction)
    ) %>%
    filter(!is_bad_genus(Genus))
  
  top_healthy <- res_df %>%
    filter(Direction == "Healthy > Diseased") %>%
    arrange(padj, desc(abs(log2FoldChange))) %>%
    slice_head(n = top_n_per_direction)
  
  top_diseased <- res_df %>%
    filter(Direction == "Diseased > Healthy") %>%
    arrange(padj, desc(abs(log2FoldChange))) %>%
    slice_head(n = top_n_per_direction)
  
  plot_df <- bind_rows(top_healthy, top_diseased) %>%
    distinct(Genus, .keep_all = TRUE)
  
  if (nrow(plot_df) == 0) {
    stop("No plottable genera remained for: ", short_label)
  }
  
  genus_order_desc <- plot_df %>%
    arrange(desc(log2FoldChange), padj) %>%
    pull(Genus)
  
  genus_order_plot <- rev(genus_order_desc)
  
  # -------------------------------------------------------
  # 4.3 Subset phyloseq for abundance panel
  # -------------------------------------------------------
  md <- as(sample_data(ps0), "data.frame")
  
  keep_ids <- rownames(md)[
    as.character(md$Month) == month_name &
      as.character(md$Compartment) == compartment_name &
      as.character(md$Health_status) %in% c(healthy_label, diseased_label)
  ]
  
  ps_sub <- prune_samples(keep_ids, ps0)
  ps_sub <- prune_samples(sample_sums(ps_sub) > 0, ps_sub)
  
  md_sub <- as(sample_data(ps_sub), "data.frame")
  md_sub$Condition <- ifelse(
    as.character(md_sub$Health_status) == healthy_label,
    "Healthy",
    "Diseased"
  )
  md_sub$Condition <- factor(md_sub$Condition, levels = c("Healthy", "Diseased"))
  sample_data(ps_sub) <- sample_data(md_sub)
  
  count_mat <- as(otu_table(ps_sub), "matrix")
  if (!taxa_are_rows(ps_sub)) {
    count_mat <- t(count_mat)
  }
  
  tax_df <- as.data.frame(tax_table(ps_sub)) %>%
    rownames_to_column("ASV_hash")
  
  tax_df$Genus <- clean_genus_label(tax_df$Genus)
  
  keep_rows <- !is_bad_genus(tax_df$Genus)
  count_mat <- count_mat[keep_rows, , drop = FALSE]
  tax_df <- tax_df[keep_rows, , drop = FALSE]
  
  tax_df <- tax_df[match(rownames(count_mat), tax_df$ASV_hash), , drop = FALSE]
  if (!identical(rownames(count_mat), tax_df$ASV_hash)) {
    stop("ASV hash alignment failed in abundance panel for: ", short_label)
  }
  
  genus_count_mat <- rowsum(
    x = count_mat,
    group = tax_df$Genus,
    reorder = FALSE
  )
  
  genus_rel_mat <- sweep(genus_count_mat, 2, colSums(genus_count_mat), "/")
  genus_rel_mat[is.na(genus_rel_mat)] <- 0
  
  abund_long <- as.data.frame(genus_rel_mat) %>%
    rownames_to_column("Genus") %>%
    pivot_longer(
      cols = -Genus,
      names_to = "SampleID",
      values_to = "RelAbund"
    ) %>%
    left_join(
      md_sub %>%
        rownames_to_column("SampleID") %>%
        select(SampleID, Condition),
      by = "SampleID"
    ) %>%
    filter(Genus %in% plot_df$Genus)
  
  bars_df <- abund_long %>%
    group_by(Genus, Condition) %>%
    summarise(
      n = n(),
      mean_val = mean(RelAbund, na.rm = TRUE),
      sd_val   = sd(RelAbund, na.rm = TRUE),
      se_val   = sd_val / sqrt(n),
      t_val    = qt(0.975, df = pmax(n - 1, 1)),
      lo_val   = mean_val - t_val * se_val,
      hi_val   = mean_val + t_val * se_val,
      .groups  = "drop"
    ) %>%
    mutate(
      mean_pct = 100 * mean_val,
      lo_pct   = 100 * pmax(lo_val, 0),
      hi_pct   = 100 * pmin(hi_val, 1),
      Condition = factor(Condition, levels = c("Healthy", "Diseased")),
      Genus = factor(Genus, levels = genus_order_plot)
    )
  
  # -------------------------------------------------------
  # 4.4 Data for effect-size panel
  # -------------------------------------------------------
  effect_df <- plot_df %>%
    mutate(
      Genus = factor(Genus, levels = genus_order_plot),
      p_label = p_label_fun(padj),
      Direction_fill = dplyr::case_when(
        Direction == "Healthy > Diseased" ~ "Healthy",
        Direction == "Diseased > Healthy" ~ "Diseased",
        TRUE ~ "NS"
      )
    )
  
  # -------------------------------------------------------
  # 4.5 Background bands
  # -------------------------------------------------------
  lvl <- levels(bars_df$Genus)
  bg_df <- data.frame(
    Genus = factor(lvl, levels = lvl),
    idx = seq_along(lvl)
  ) %>%
    mutate(
      is_alt = idx %% 2 == 0,
      xmin = idx - 0.5,
      xmax = idx + 0.5
    )
  
  # -------------------------------------------------------
  # 4.6 Axis ranges
  # -------------------------------------------------------
  max_pct <- max(bars_df$hi_pct, bars_df$mean_pct, na.rm = TRUE)
  y_limit_bar <- ceiling(max_pct * 1.20 * 10) / 10
  
  y_min_eff <- min(effect_df$log2FoldChange - 1.96 * effect_df$lfcSE, na.rm = TRUE)
  y_max_eff <- max(effect_df$log2FoldChange + 1.96 * effect_df$lfcSE, na.rm = TRUE)
  y_rng_eff <- y_max_eff - y_min_eff
  off_eff <- if (is.finite(y_rng_eff) && y_rng_eff > 0) 0.08 * y_rng_eff else 0.3
  
  # -------------------------------------------------------
  # 4.7 Panel A — Mean relative abundance
  # -------------------------------------------------------
  pA <- ggplot() +
    geom_rect(
      data = bg_df,
      aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = is_alt),
      inherit.aes = FALSE
    ) +
    scale_fill_manual(values = c(`TRUE` = BAND_LIGHT, `FALSE` = "white"), guide = "none") +
    ggnewscale::new_scale_fill() +
    geom_col(
      data = bars_df,
      aes(x = Genus, y = mean_pct, fill = Condition),
      position = position_dodge(width = 0.8),
      width = 0.75,
      colour = "black",
      linewidth = 0.55
    ) +
    geom_errorbar(
      data = bars_df,
      aes(x = Genus, ymin = lo_pct, ymax = hi_pct, group = Condition),
      position = position_dodge(width = 0.8),
      width = 0.22,
      colour = "black",
      linewidth = 0.45
    ) +
    scale_fill_manual(values = c("Healthy" = COL_HEALTHY, "Diseased" = COL_DISEASED)) +
    scale_y_continuous(
      limits = c(0, y_limit_bar),
      breaks = pretty_breaks(n = 4)
    ) +
    labs(
      x = NULL,
      y = "Mean proportion (%)",
      title = paste0(month_name, " · ", compartment_name, " — top significant genera"),
      subtitle = paste0("Top ", top_n_per_direction, " Healthy + top ", top_n_per_direction, " Diseased")
    ) +
    coord_flip(clip = "off") +
    theme_big_labels() +
    theme(legend.position = "top")
  
  # -------------------------------------------------------
  # 4.8 Panel B — log2FoldChange
  # -------------------------------------------------------
  pB <- ggplot() +
    geom_rect(
      data = bg_df,
      aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = is_alt),
      inherit.aes = FALSE
    ) +
    scale_fill_manual(values = c(`TRUE` = BAND_LIGHT, `FALSE` = "white"), guide = "none") +
    ggnewscale::new_scale_fill() +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_errorbar(
      data = effect_df,
      aes(
        x = Genus,
        ymin = log2FoldChange - 1.96 * lfcSE,
        ymax = log2FoldChange + 1.96 * lfcSE
      ),
      width = 0.2,
      colour = "black",
      linewidth = 0.45
    ) +
    geom_point(
      data = effect_df,
      aes(x = Genus, y = log2FoldChange, fill = Direction_fill),
      shape = 21,
      colour = "black",
      stroke = 0.55,
      size = DOT_SIZE
    ) +
    scale_fill_manual(values = c("Healthy" = COL_HEALTHY, "Diseased" = COL_DISEASED)) +
    geom_text(
      data = effect_df,
      aes(x = Genus, y = y_max_eff + off_eff, label = p_label),
      hjust = 0,
      vjust = 0.5,
      size = 3.1
    ) +
    scale_y_continuous(
      limits = c(y_min_eff, y_max_eff + 3 * off_eff),
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    labs(
      x = NULL,
      y = "log2FoldChange (Diseased vs Healthy)"
    ) +
    coord_flip(clip = "off") +
    theme_big_labels() +
    theme(
      legend.position = "top",
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      plot.title = element_blank(),
      plot.subtitle = element_blank()
    )
  
  # -------------------------------------------------------
  # 4.9 Combine
  # -------------------------------------------------------
  combined <- wrap_plots(
    pA, pB,
    nrow = 1,
    widths = c(2, 1),
    guides = "collect"
  )
  
  fig_height <- max(6.5, 0.24 * length(genus_order_plot) + 2.2)
  
  ggsave(
    filename = file.path(outdir_fig, paste0("FIG_DESeq2_Genus_", short_label, ".png")),
    plot = combined,
    width = 11,
    height = fig_height,
    dpi = 300
  )
  
  ggsave(
    filename = file.path(outdir_fig, paste0("FIG_DESeq2_Genus_", short_label, ".pdf")),
    plot = combined,
    width = 11,
    height = fig_height
  )
  
  cat("\n====================================================\n")
  cat("SELECTED GENERA FOR:", short_label, "\n")
  cat("====================================================\n")
  print(
    plot_df %>%
      select(Genus, log2FoldChange, padj, Direction) %>%
      arrange(desc(log2FoldChange), padj),
    n = Inf
  )
  
  return(combined)
}

############################################################
# 5) BUILD ALL 6 FIGURES
############################################################

plot_list <- list()

for (i in seq_len(nrow(comparisons_tbl))) {
  row_i <- comparisons_tbl[i, ]
  
  fig_i <- make_deseq_genus_figure(
    month_name = as.character(row_i$Month),
    compartment_name = as.character(row_i$Compartment),
    healthy_label = as.character(row_i$Group_Healthy),
    diseased_label = as.character(row_i$Group_Diseased),
    short_label = as.character(row_i$short_label),
    top_n_per_direction = TOP_N_PER_DIRECTION
  )
  
  plot_list[[as.character(row_i$short_label)]] <- fig_i
  print(fig_i)
}

############################################################
# 6) COMBINED PANEL OF THE 6 COMPARISONS
############################################################

panel_all <- wrap_plots(plot_list, ncol = 1)

ggsave(
  filename = file.path(outdir_fig, "FIG_DESeq2_Genus_PANEL_all_comparisons.png"),
  plot = panel_all,
  width = 11,
  height = 38,
  dpi = 300
)

ggsave(
  filename = file.path(outdir_fig, "FIG_DESeq2_Genus_PANEL_all_comparisons.pdf"),
  plot = panel_all,
  width = 11,
  height = 38
)

cat("\n====================================================\n")
cat("FILES WRITTEN TO:\n")
cat("====================================================\n")
print(normalizePath(outdir_fig))



############################################################
# SECTION 7 — FIGURES FOR DIFFERENTIAL ABUNDANCE (16S)
# DESEQ2 RESULTS — ASV LEVEL
#
# Objetivos:
# 1) Graficar las 6 comparaciones Healthy vs Diseased
# 2) Mantener el estilo previo:
#    - Panel A: mean relative abundance (%)
#    - Panel B: log2FoldChange + error bars
# 3) Mostrar top 15 ASVs enriched in Healthy + top 15 enriched in Diseased
# 4) Mostrar taxonomía en la etiqueta:
#    - si hay Genus informativo, mostrar solo Genus
#    - si no, usar el mejor fallback taxonómico disponible
#    - si hay etiquetas repetidas, añadir hash corto para distinguir
############################################################

suppressPackageStartupMessages({
  library(phyloseq)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(patchwork)
  library(ggnewscale)
  library(scales)
  library(stringr)
  library(tibble)
})

############################################################
# 0) INPUTS AND OUTPUT FOLDERS
############################################################

physeq_path <- "Phyloseq_16S.rds"

deseq_dir <- "04_differential_abundance_16S_DESeq2"
deseq_asv_dir <- file.path(deseq_dir, "ASV")

outdir_fig_asv <- "08_DESeq2_figures_ASV"
dir.create(outdir_fig_asv, recursive = TRUE, showWarnings = FALSE)

############################################################
# 1) SETTINGS
############################################################

comparisons_tbl <- tibble::tribble(
  ~Month,     ~Compartment,    ~Group_Healthy,             ~Group_Diseased,            ~short_label,
  "July",     "Rhizosphere",   "Healthy Rhizosphere",      "Diseased Rhizosphere",     "July_Rhizosphere",
  "October",  "Rhizosphere",   "Healthy Rhizosphere",      "Diseased Rhizosphere",     "October_Rhizosphere",
  "July",     "Endosphere",    "Healthy Endosphere",       "Diseased Endosphere",      "July_Endosphere",
  "October",  "Endosphere",    "Healthy Endosphere",       "Diseased Endosphere",      "October_Endosphere",
  "July",     "Phyllosphere",  "Healthy Phyllosphere",     "Diseased Phyllosphere",    "July_Phyllosphere",
  "October",  "Phyllosphere",  "Healthy Phyllosphere",     "Diseased Phyllosphere",    "October_Phyllosphere"
)

TOP_N_PER_DIRECTION <- 15

COL_HEALTHY  <- "#A6DBA0"
COL_DISEASED <- "#F4A6A6"
BAND_LIGHT   <- "#F2F2F2"
DOT_SIZE     <- 4.8
BASE_SIZE    <- 12

############################################################
# 2) HELPERS
############################################################

clean_tax_label <- function(x) {
  x <- as.character(x)
  x[is.na(x) | x == "" | x == "NA"] <- NA_character_
  x
}

is_bad_taxon <- function(x) {
  x <- as.character(x)
  is.na(x) |
    x == "" |
    str_detect(x, regex("^unclassified", ignore_case = TRUE)) |
    str_detect(x, regex("^uncultured$", ignore_case = TRUE)) |
    str_detect(x, regex("uncultured", ignore_case = TRUE))
}

p_label_fun <- function(x) {
  dplyr::case_when(
    is.na(x)  ~ "NA",
    x < 0.001 ~ "<0.001",
    x < 0.01  ~ "<0.01",
    x < 0.05  ~ "<0.05",
    TRUE      ~ paste0("=", formatC(x, format = "g", digits = 2))
  )
}

theme_big_labels <- function() {
  theme_minimal(base_size = BASE_SIZE) +
    theme(
      panel.background = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y  = element_text(size = 12, margin = margin(r = 6)),
      axis.title.x = element_text(size = 12),
      legend.position = "top",
      legend.title = element_blank(),
      plot.margin = unit(c(6, 6, 6, 24), "pt")
    )
}

# Etiqueta taxonómica para ASV:
# prioridad = Genus > Family > Order > Class > Phylum > hash corto
build_asv_tax_label <- function(df) {
  df %>%
    mutate(
      Genus  = clean_tax_label(Genus),
      Family = clean_tax_label(Family),
      Order  = clean_tax_label(Order),
      Class  = clean_tax_label(Class),
      Phylum = clean_tax_label(Phylum),
      Species = clean_tax_label(Species)
    ) %>%
    mutate(
      tax_label_raw = case_when(
        !is_bad_taxon(Genus)  ~ Genus,
        !is_bad_taxon(Family) ~ paste0("Family: ", Family),
        !is_bad_taxon(Order)  ~ paste0("Order: ", Order),
        !is_bad_taxon(Class)  ~ paste0("Class: ", Class),
        !is_bad_taxon(Phylum) ~ paste0("Phylum: ", Phylum),
        TRUE ~ paste0("ASV ", substr(ASV_hash, 1, 8))
      )
    ) %>%
    group_by(tax_label_raw) %>%
    mutate(
      n_same_label = n(),
      tax_label = if_else(
        n_same_label == 1,
        tax_label_raw,
        paste0(tax_label_raw, " (", substr(ASV_hash, 1, 8), ")")
      )
    ) %>%
    ungroup()
}

############################################################
# 3) LOAD PHYLOSEQ
############################################################

ps0 <- readRDS(physeq_path)
ps0 <- prune_samples(sample_sums(ps0) > 0, ps0)

############################################################
# 4) FUNCTION TO BUILD ONE FIGURE
############################################################

make_deseq_asv_figure <- function(month_name,
                                  compartment_name,
                                  healthy_label,
                                  diseased_label,
                                  short_label,
                                  top_n_per_direction = 15) {
  
  # -------------------------------------------------------
  # 4.1 Read DESeq2 ASV results for this comparison
  # -------------------------------------------------------
  res_file <- file.path(
    deseq_asv_dir,
    paste0("DESeq2_ASV_", short_label, "_significant.csv")
  )
  
  if (!file.exists(res_file)) {
    stop("Missing file: ", res_file)
  }
  
  res_df <- read_csv(res_file, show_col_types = FALSE)
  
  # -------------------------------------------------------
  # 4.2 Build plotting labels
  # -------------------------------------------------------
  res_df <- build_asv_tax_label(res_df)
  
  top_healthy <- res_df %>%
    filter(Direction == "Healthy > Diseased") %>%
    arrange(padj, desc(abs(log2FoldChange))) %>%
    slice_head(n = top_n_per_direction)
  
  top_diseased <- res_df %>%
    filter(Direction == "Diseased > Healthy") %>%
    arrange(padj, desc(abs(log2FoldChange))) %>%
    slice_head(n = top_n_per_direction)
  
  plot_df <- bind_rows(top_healthy, top_diseased) %>%
    distinct(ASV_hash, .keep_all = TRUE)
  
  if (nrow(plot_df) == 0) {
    stop("No plottable ASVs remained for: ", short_label)
  }
  
  asv_order_desc <- plot_df %>%
    arrange(desc(log2FoldChange), padj) %>%
    pull(tax_label)
  
  asv_order_plot <- rev(asv_order_desc)
  
  # -------------------------------------------------------
  # 4.3 Subset phyloseq for abundance panel
  # -------------------------------------------------------
  md <- as(sample_data(ps0), "data.frame")
  
  keep_ids <- rownames(md)[
    as.character(md$Month) == month_name &
      as.character(md$Compartment) == compartment_name &
      as.character(md$Health_status) %in% c(healthy_label, diseased_label)
  ]
  
  ps_sub <- prune_samples(keep_ids, ps0)
  ps_sub <- prune_samples(sample_sums(ps_sub) > 0, ps_sub)
  
  md_sub <- as(sample_data(ps_sub), "data.frame")
  md_sub$Condition <- ifelse(
    as.character(md_sub$Health_status) == healthy_label,
    "Healthy",
    "Diseased"
  )
  md_sub$Condition <- factor(md_sub$Condition, levels = c("Healthy", "Diseased"))
  sample_data(ps_sub) <- sample_data(md_sub)
  
  count_mat <- as(otu_table(ps_sub), "matrix")
  if (!taxa_are_rows(ps_sub)) {
    count_mat <- t(count_mat)
  }
  
  tax_df <- as.data.frame(tax_table(ps_sub)) %>%
    rownames_to_column("ASV_hash")
  
  tax_df <- build_asv_tax_label(tax_df)
  
  # Mantener solo los ASVs seleccionados
  count_mat <- count_mat[rownames(count_mat) %in% plot_df$ASV_hash, , drop = FALSE]
  tax_df <- tax_df %>% filter(ASV_hash %in% rownames(count_mat))
  
  tax_df <- tax_df[match(rownames(count_mat), tax_df$ASV_hash), , drop = FALSE]
  if (!identical(rownames(count_mat), tax_df$ASV_hash)) {
    stop("ASV hash alignment failed in abundance panel for: ", short_label)
  }
  
  # Relative abundance per sample
  rel_mat <- sweep(count_mat, 2, colSums(count_mat), "/")
  rel_mat[is.na(rel_mat)] <- 0
  
  abund_long <- as.data.frame(rel_mat) %>%
    rownames_to_column("ASV_hash") %>%
    pivot_longer(
      cols = -ASV_hash,
      names_to = "SampleID",
      values_to = "RelAbund"
    ) %>%
    left_join(
      tax_df %>% select(ASV_hash, tax_label),
      by = "ASV_hash"
    ) %>%
    left_join(
      md_sub %>%
        rownames_to_column("SampleID") %>%
        select(SampleID, Condition),
      by = "SampleID"
    )
  
  bars_df <- abund_long %>%
    group_by(ASV_hash, tax_label, Condition) %>%
    summarise(
      n = n(),
      mean_val = mean(RelAbund, na.rm = TRUE),
      sd_val   = sd(RelAbund, na.rm = TRUE),
      se_val   = sd_val / sqrt(n),
      t_val    = qt(0.975, df = pmax(n - 1, 1)),
      lo_val   = mean_val - t_val * se_val,
      hi_val   = mean_val + t_val * se_val,
      .groups  = "drop"
    ) %>%
    mutate(
      mean_pct = 100 * mean_val,
      lo_pct   = 100 * pmax(lo_val, 0),
      hi_pct   = 100 * pmin(hi_val, 1),
      Condition = factor(Condition, levels = c("Healthy", "Diseased")),
      tax_label = factor(tax_label, levels = asv_order_plot)
    )
  
  # -------------------------------------------------------
  # 4.4 Data for effect-size panel
  # -------------------------------------------------------
  effect_df <- plot_df %>%
    mutate(
      tax_label = factor(tax_label, levels = asv_order_plot),
      p_label = p_label_fun(padj),
      Direction_fill = dplyr::case_when(
        Direction == "Healthy > Diseased" ~ "Healthy",
        Direction == "Diseased > Healthy" ~ "Diseased",
        TRUE ~ "NS"
      )
    )
  
  # -------------------------------------------------------
  # 4.5 Background bands
  # -------------------------------------------------------
  lvl <- levels(bars_df$tax_label)
  bg_df <- data.frame(
    tax_label = factor(lvl, levels = lvl),
    idx = seq_along(lvl)
  ) %>%
    mutate(
      is_alt = idx %% 2 == 0,
      xmin = idx - 0.5,
      xmax = idx + 0.5
    )
  
  # -------------------------------------------------------
  # 4.6 Axis ranges
  # -------------------------------------------------------
  max_pct <- max(bars_df$hi_pct, bars_df$mean_pct, na.rm = TRUE)
  y_limit_bar <- ceiling(max_pct * 1.20 * 10) / 10
  
  y_min_eff <- min(effect_df$log2FoldChange - 1.96 * effect_df$lfcSE, na.rm = TRUE)
  y_max_eff <- max(effect_df$log2FoldChange + 1.96 * effect_df$lfcSE, na.rm = TRUE)
  y_rng_eff <- y_max_eff - y_min_eff
  off_eff <- if (is.finite(y_rng_eff) && y_rng_eff > 0) 0.08 * y_rng_eff else 0.3
  
  # -------------------------------------------------------
  # 4.7 Panel A — Mean relative abundance
  # -------------------------------------------------------
  pA <- ggplot() +
    geom_rect(
      data = bg_df,
      aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = is_alt),
      inherit.aes = FALSE
    ) +
    scale_fill_manual(values = c(`TRUE` = BAND_LIGHT, `FALSE` = "white"), guide = "none") +
    ggnewscale::new_scale_fill() +
    geom_col(
      data = bars_df,
      aes(x = tax_label, y = mean_pct, fill = Condition),
      position = position_dodge(width = 0.8),
      width = 0.75,
      colour = "black",
      linewidth = 0.55
    ) +
    geom_errorbar(
      data = bars_df,
      aes(x = tax_label, ymin = lo_pct, ymax = hi_pct, group = Condition),
      position = position_dodge(width = 0.8),
      width = 0.22,
      colour = "black",
      linewidth = 0.45
    ) +
    scale_fill_manual(values = c("Healthy" = COL_HEALTHY, "Diseased" = COL_DISEASED)) +
    scale_y_continuous(
      limits = c(0, y_limit_bar),
      breaks = pretty_breaks(n = 4)
    ) +
    labs(
      x = NULL,
      y = "Mean proportion (%)",
      title = paste0(month_name, " · ", compartment_name, " — top significant ASVs"),
      subtitle = paste0("Top ", top_n_per_direction, " Healthy + top ", top_n_per_direction, " Diseased")
    ) +
    coord_flip(clip = "off") +
    theme_big_labels() +
    theme(legend.position = "top")
  
  # -------------------------------------------------------
  # 4.8 Panel B — log2FoldChange
  # -------------------------------------------------------
  pB <- ggplot() +
    geom_rect(
      data = bg_df,
      aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = is_alt),
      inherit.aes = FALSE
    ) +
    scale_fill_manual(values = c(`TRUE` = BAND_LIGHT, `FALSE` = "white"), guide = "none") +
    ggnewscale::new_scale_fill() +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_errorbar(
      data = effect_df,
      aes(
        x = tax_label,
        ymin = log2FoldChange - 1.96 * lfcSE,
        ymax = log2FoldChange + 1.96 * lfcSE
      ),
      width = 0.2,
      colour = "black",
      linewidth = 0.45
    ) +
    geom_point(
      data = effect_df,
      aes(x = tax_label, y = log2FoldChange, fill = Direction_fill),
      shape = 21,
      colour = "black",
      stroke = 0.55,
      size = DOT_SIZE
    ) +
    scale_fill_manual(values = c("Healthy" = COL_HEALTHY, "Diseased" = COL_DISEASED)) +
    geom_text(
      data = effect_df,
      aes(x = tax_label, y = y_max_eff + off_eff, label = p_label),
      hjust = 0,
      vjust = 0.5,
      size = 3.1
    ) +
    scale_y_continuous(
      limits = c(y_min_eff, y_max_eff + 3 * off_eff),
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    labs(
      x = NULL,
      y = "log2FoldChange (Diseased vs Healthy)"
    ) +
    coord_flip(clip = "off") +
    theme_big_labels() +
    theme(
      legend.position = "top",
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      plot.title = element_blank(),
      plot.subtitle = element_blank()
    )
  
  # -------------------------------------------------------
  # 4.9 Combine
  # -------------------------------------------------------
  combined <- wrap_plots(
    pA, pB,
    nrow = 1,
    widths = c(2, 1),
    guides = "collect"
  )
  
  fig_height <- max(6.5, 0.24 * length(asv_order_plot) + 2.2)
  
  ggsave(
    filename = file.path(outdir_fig_asv, paste0("FIG_DESeq2_ASV_", short_label, ".png")),
    plot = combined,
    width = 11,
    height = fig_height,
    dpi = 300
  )
  
  ggsave(
    filename = file.path(outdir_fig_asv, paste0("FIG_DESeq2_ASV_", short_label, ".pdf")),
    plot = combined,
    width = 11,
    height = fig_height
  )
  
  cat("\n====================================================\n")
  cat("SELECTED ASVs FOR:", short_label, "\n")
  cat("====================================================\n")
  print(
    plot_df %>%
      select(ASV_hash, tax_label, log2FoldChange, padj, Direction) %>%
      arrange(desc(log2FoldChange), padj),
    n = Inf
  )
  
  return(combined)
}

############################################################
# 5) BUILD ALL 6 FIGURES
############################################################

plot_list_asv <- list()

for (i in seq_len(nrow(comparisons_tbl))) {
  row_i <- comparisons_tbl[i, ]
  
  fig_i <- make_deseq_asv_figure(
    month_name = as.character(row_i$Month),
    compartment_name = as.character(row_i$Compartment),
    healthy_label = as.character(row_i$Group_Healthy),
    diseased_label = as.character(row_i$Group_Diseased),
    short_label = as.character(row_i$short_label),
    top_n_per_direction = TOP_N_PER_DIRECTION
  )
  
  plot_list_asv[[as.character(row_i$short_label)]] <- fig_i
  print(fig_i)
}

############################################################
# 6) COMBINED PANEL OF THE 6 COMPARISONS
############################################################

panel_all_asv <- wrap_plots(plot_list_asv, ncol = 1)

ggsave(
  filename = file.path(outdir_fig_asv, "FIG_DESeq2_ASV_PANEL_all_comparisons.png"),
  plot = panel_all_asv,
  width = 11,
  height = 38,
  dpi = 300
)

ggsave(
  filename = file.path(outdir_fig_asv, "FIG_DESeq2_ASV_PANEL_all_comparisons.pdf"),
  plot = panel_all_asv,
  width = 11,
  height = 38
)

cat("\n====================================================\n")
cat("FILES WRITTEN TO:\n")
cat("====================================================\n")
print(normalizePath(outdir_fig_asv))







cat("\n=== 16S SUMMARY COUNTS ===\n")
print(sig_counts16, n = Inf, width = Inf)

cat("\n=== 16S TOP 10 PER CONTRAST ===\n")
deseq16 %>%
  group_by(Month, Compartment, result_summary) %>%
  arrange(padj, desc(abs(log2FoldChange))) %>%
  slice_head(n = 10) %>%
  ungroup() %>%
  select(Month, Compartment, result_summary, Genus, Phylum, log2FoldChange, padj) %>%
  arrange(Month, Compartment, result_summary, padj) %>%
  print(n = Inf, width = Inf)

cat("\n=== 16S RECURRENT HEALTHY GENERA n>=2 ===\n")
deseq16 %>%
  filter(result_summary == "Healthy significantly enriched") %>%
  group_by(Genus, Phylum) %>%
  summarise(
    n_contrasts = n(),
    compartments = paste(unique(Compartment), collapse = "; "),
    months = paste(unique(Month), collapse = "; "),
    min_padj = min(padj, na.rm = TRUE),
    max_abs_log2FC = max(abs(log2FoldChange), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(n_contrasts >= 2) %>%
  arrange(desc(n_contrasts), min_padj) %>%
  print(n = Inf, width = Inf)

cat("\n=== 16S RECURRENT DISEASED GENERA n>=2 ===\n")
deseq16 %>%
  filter(result_summary == "Diseased significantly enriched") %>%
  group_by(Genus, Phylum) %>%
  summarise(
    n_contrasts = n(),
    compartments = paste(unique(Compartment), collapse = "; "),
    months = paste(unique(Month), collapse = "; "),
    min_padj = min(padj, na.rm = TRUE),
    max_abs_log2FC = max(abs(log2FoldChange), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(n_contrasts >= 2) %>%
  arrange(desc(n_contrasts), min_padj) %>%
  print(n = Inf, width = Inf)


####################################################################################################################################################################################
####################################################################################################################################################################################
####################################################################################################################################################################################
####################################################################################################################################################################################






############################################################
# TABLE S4 — 16S BETA DIVERSITY MASTER TABLE
# Combina global + compartment + UniFrac
############################################################

setwd("path/to/QIIME2_files_to_R_16S")

library(readr)
library(dplyr)
library(tidyr)

outdir_beta_tables <- "02_beta_diversity_16S/tables"

beta16_global  <- read_csv(file.path(outdir_beta_tables, "01_beta_global_summary_July_October.csv"), show_col_types = FALSE)
beta16_comp    <- read_csv(file.path(outdir_beta_tables, "02_beta_by_compartment_summary_July_October.csv"), show_col_types = FALSE)
beta16_unifrac <- read_csv(file.path(outdir_beta_tables, "03_beta_unifrac_by_compartment_July_October.csv"), show_col_types = FALSE)

############################################################
# 1) GLOBAL BRAY / JACCARD
############################################################

global_tbl <- beta16_global %>%
  transmute(
    Analysis_level = "Global",
    Month,
    Compartment = NA_character_,
    Metric,
    n_samples,
    n_Healthy = NA_real_,
    n_Diseased = NA_real_,
    PCoA_axis1_var = round(PCoA_axis1_var, 2),
    PCoA_axis2_var = round(PCoA_axis2_var, 2),
    PERMANOVA_R2 = round(PERMANOVA_R2, 3),
    PERMANOVA_F = round(PERMANOVA_F, 3),
    PERMANOVA_p = signif(PERMANOVA_p, 3),
    ANOSIM_R = round(ANOSIM_R, 3),
    ANOSIM_p = signif(ANOSIM_p, 3),
    betadisper_F = round(betadisper_F, 3),
    betadisper_p = signif(betadisper_p, 3)
  )

comp_bray <- beta16_comp %>%
  transmute(
    Analysis_level = "Within compartment",
    Month,
    Compartment,
    Metric = "Bray-Curtis",
    n_samples,
    n_Healthy,
    n_Diseased,
    PCoA_axis1_var = round(Bray_PCoA_axis1_var, 2),
    PCoA_axis2_var = round(Bray_PCoA_axis2_var, 2),
    PERMANOVA_R2 = round(Bray_PERMANOVA_R2, 3),
    PERMANOVA_F = round(Bray_PERMANOVA_F, 3),
    PERMANOVA_p = signif(Bray_PERMANOVA_p, 3),
    ANOSIM_R = round(Bray_ANOSIM_R, 3),
    ANOSIM_p = signif(Bray_ANOSIM_p, 3),
    betadisper_F = round(Bray_betadisper_F, 3),
    betadisper_p = signif(Bray_betadisper_p, 3)
  )

comp_jaccard <- beta16_comp %>%
  transmute(
    Analysis_level = "Within compartment",
    Month,
    Compartment,
    Metric = "Jaccard",
    n_samples,
    n_Healthy,
    n_Diseased,
    PCoA_axis1_var = round(Jaccard_PCoA_axis1_var, 2),
    PCoA_axis2_var = round(Jaccard_PCoA_axis2_var, 2),
    PERMANOVA_R2 = round(Jaccard_PERMANOVA_R2, 3),
    PERMANOVA_F = round(Jaccard_PERMANOVA_F, 3),
    PERMANOVA_p = signif(Jaccard_PERMANOVA_p, 3),
    ANOSIM_R = round(Jaccard_ANOSIM_R, 3),
    ANOSIM_p = signif(Jaccard_ANOSIM_p, 3),
    betadisper_F = round(Jaccard_betadisper_F, 3),
    betadisper_p = signif(Jaccard_betadisper_p, 3)
  )

unifrac_tbl <- beta16_unifrac %>%
  transmute(
    Analysis_level = "Within compartment",
    Month,
    Compartment,
    Metric,
    n_samples = NA_real_,
    n_Healthy = NA_real_,
    n_Diseased = NA_real_,
    PCoA_axis1_var = round(PCoA_axis1_var, 2),
    PCoA_axis2_var = round(PCoA_axis2_var, 2),
    PERMANOVA_R2 = round(PERMANOVA_R2, 3),
    PERMANOVA_F = round(PERMANOVA_F, 3),
    PERMANOVA_p = signif(PERMANOVA_p, 3),
    ANOSIM_R = round(ANOSIM_R, 3),
    ANOSIM_p = signif(ANOSIM_p, 3),
    betadisper_F = round(betadisper_F, 3),
    betadisper_p = signif(betadisper_p, 3)
  )

beta16_master <- bind_rows(
  global_tbl,
  comp_bray,
  comp_jaccard,
  unifrac_tbl
) %>%
  arrange(
    Month,
    factor(Analysis_level, levels = c("Global", "Within compartment")),
    factor(Compartment, levels = c(NA, "Rhizosphere", "Endosphere", "Phyllosphere")),
    factor(Metric, levels = c("Bray-Curtis", "Jaccard", "Weighted UniFrac", "Unweighted UniFrac"))
  ) %>%
  mutate(across(everything(), ~ ifelse(is.na(.x), "-", as.character(.x))))

write_csv(
  beta16_master,
  file.path(outdir_beta_tables, "Table_S4_beta_diversity_16S_master.csv")
)

print(beta16_master, n = Inf, width = Inf)













############################################################
# SECTION 3B — TREE SUBSAMPLING EFFECT CHECK (16S)
#
# Objetivos:
# 1) Evaluar si las submuestras de un mismo árbol se agrupan entre sí
# 2) Generar PCoAs Bray-Curtis coloreados por Health_status y por árbol
# 3) Testar el efecto de Tree dentro de cada Month x Compartment x Health_status
# 4) Comparar distancias Within-tree vs Between-tree dentro de cada grupo
# 5) Agregar las submuestras a nivel árbol como análisis de robustez
#
# Nota:
# Este apartado NO sustituye al beta diversity principal.
# Sirve como control/sensibilidad para justificar el diseño de submuestreo.
############################################################

suppressPackageStartupMessages({
  library(phyloseq)
  library(tidyverse)
  library(vegan)
  library(ape)
  library(readr)
  library(ggplot2)
  library(patchwork)
})

############################################################
# 0) OUTPUT FOLDER
############################################################

outdir_treecheck <- "02B_tree_subsampling_check_16S"
outdir_treecheck_figs <- file.path(outdir_treecheck, "figures")
outdir_treecheck_tables <- file.path(outdir_treecheck, "tables")

dir.create(outdir_treecheck, recursive = TRUE, showWarnings = FALSE)
dir.create(outdir_treecheck_figs, recursive = TRUE, showWarnings = FALSE)
dir.create(outdir_treecheck_tables, recursive = TRUE, showWarnings = FALSE)


############################################################
# 0B) LOAD PHYLOSEQ IF ps0 DOES NOT EXIST
############################################################

# Si ejecutas este apartado de forma independiente, necesitamos cargar ps0.
# Si vienes ejecutando todo el script desde arriba, ps0 ya existirá y no se recarga.

if (!exists("ps0")) {
  
  cat("\nps0 no existe en la sesión. Cargando Phyloseq_16S.rds...\n")
  
  physeq_path <- "Phyloseq_16S.rds"
  
  if (!file.exists(physeq_path)) {
    stop("No encuentro el archivo: ", physeq_path,
         "\nComprueba que estás en el directorio correcto con getwd().")
  }
  
  ps0 <- readRDS(physeq_path)
  ps0 <- prune_samples(sample_sums(ps0) > 0, ps0)
  ps0 <- prune_taxa(taxa_sums(ps0) > 0, ps0)
}

cat("\nObjeto ps0 cargado correctamente:\n")
cat("Samples:", nsamples(ps0), "\n")
cat("Taxa:", ntaxa(ps0), "\n")
############################################################
# 1) CHECK REQUIRED METADATA COLUMNS
############################################################

# Usamos ps0 ya cargado en SECTION 3.
# Si ejecutas este bloque por separado, descomenta estas dos líneas:
# ps0 <- readRDS("Phyloseq_16S.rds")
# ps0 <- prune_samples(sample_sums(ps0) > 0, ps0)

meta_tree <- as(sample_data(ps0), "data.frame")
meta_tree$SampleID <- rownames(meta_tree)

required_cols_tree <- c("Month", "Compartment", "Health_status", "Trees")
missing_cols_tree <- setdiff(required_cols_tree, colnames(meta_tree))

if (length(missing_cols_tree) > 0) {
  stop(
    "Faltan columnas requeridas en sample_data para el análisis de submuestreo: ",
    paste(missing_cols_tree, collapse = ", ")
  )
}

# Orden de factores
meta_tree$Month <- factor(as.character(meta_tree$Month), levels = c("July", "October"))
meta_tree$Compartment <- factor(
  as.character(meta_tree$Compartment),
  levels = c("Bulk soil", "Rhizosphere", "Endosphere", "Phyllosphere")
)
meta_tree$Health_status <- factor(
  as.character(meta_tree$Health_status),
  levels = c(
    "Bulk soil",
    "Healthy Rhizosphere", "Diseased Rhizosphere",
    "Healthy Endosphere",  "Diseased Endosphere",
    "Healthy Phyllosphere","Diseased Phyllosphere"
  )
)
meta_tree$Trees <- factor(as.character(meta_tree$Trees))

sample_data(ps0) <- sample_data(meta_tree)

cat("\n====================================================\n")
cat("TREE SUBSAMPLING CHECK — BASIC CHECKS\n")
cat("====================================================\n")
cat("Metadata columns:\n")
print(colnames(meta_tree))

cat("\nTrees per Month x Compartment x Health_status:\n")
print(
  meta_tree %>%
    filter(Compartment != "Bulk soil") %>%
    group_by(Month, Compartment, Health_status) %>%
    summarise(
      n_samples = n(),
      n_trees = n_distinct(Trees),
      trees = paste(sort(unique(as.character(Trees))), collapse = "; "),
      .groups = "drop"
    ),
  n = Inf,
  width = Inf
)

############################################################
# 2) HELPER FUNCTIONS
############################################################

get_otu_samples_x_taxa <- function(ps) {
  otu <- as(otu_table(ps), "matrix")
  if (taxa_are_rows(ps)) {
    otu <- t(otu)
  }
  otu
}

make_bray_distance <- function(ps) {
  otu <- get_otu_samples_x_taxa(ps)
  vegan::vegdist(otu, method = "bray")
}

pcoa_labels_from_dist <- function(D) {
  pcoa_obj <- ape::pcoa(D)
  vals <- pcoa_obj$values$Relative_eig
  list(
    pcoa = pcoa_obj,
    xlab = paste0("PCoA1 (", round(100 * vals[1], 1), "%)"),
    ylab = paste0("PCoA2 (", round(100 * vals[2], 1), "%)")
  )
}

make_health_binary_from_health_status <- function(x) {
  case_when(
    grepl("^Healthy", as.character(x)) ~ "Healthy",
    grepl("^Diseased", as.character(x)) ~ "Diseased",
    TRUE ~ as.character(x)
  )
}

############################################################
# 3) PCoA BY MONTH AND COMPARTMENT — SAMPLE LEVEL
############################################################

run_treecheck_pcoa_one <- function(ps_in, month_name, compartment_name) {
  
  md <- as(sample_data(ps_in), "data.frame")
  
  keep_ids <- rownames(md)[
    as.character(md$Month) == month_name &
      as.character(md$Compartment) == compartment_name
  ]
  
  ps_sub <- prune_samples(keep_ids, ps_in)
  ps_sub <- prune_samples(sample_sums(ps_sub) > 0, ps_sub)
  ps_sub <- prune_taxa(taxa_sums(ps_sub) > 0, ps_sub)
  
  md_sub <- as(sample_data(ps_sub), "data.frame")
  md_sub$SampleID <- rownames(md_sub)
  md_sub$Health_binary <- factor(
    make_health_binary_from_health_status(md_sub$Health_status),
    levels = c("Healthy", "Diseased")
  )
  md_sub$Trees <- factor(as.character(md_sub$Trees))
  
  if (nsamples(ps_sub) < 4 || length(unique(md_sub$Trees)) < 2) {
    warning("Skipping ", month_name, " - ", compartment_name, ": not enough samples or trees.")
    return(NULL)
  }
  
  D_bray <- make_bray_distance(ps_sub)
  pcoa_info <- pcoa_labels_from_dist(D_bray)
  pcoa_obj <- pcoa_info$pcoa
  
  pcoa_df <- as.data.frame(pcoa_obj$vectors[, 1:2, drop = FALSE])
  pcoa_df$SampleID <- rownames(pcoa_df)
  pcoa_df <- left_join(pcoa_df, md_sub, by = "SampleID")
  
  # Plot 1: color = Health, shape = Tree
  p_health_tree <- ggplot(
    pcoa_df,
    aes(x = Axis.1, y = Axis.2, color = Health_binary, shape = Trees)
  ) +
    geom_point(size = 3.5, alpha = 0.9) +
    labs(
      title = paste0("16S Bray-Curtis PCoA — ", month_name, " · ", compartment_name),
      subtitle = "Colour = health status; shape = tree identity",
      x = pcoa_info$xlab,
      y = pcoa_info$ylab,
      color = "Health status",
      shape = "Tree"
    ) +
    scale_color_manual(values = c("Healthy" = "#4B9EA7", "Diseased" = "#801723")) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold")
    )
  
  # Plot 2: color = Tree, facet = Health
  p_tree_facet <- ggplot(
    pcoa_df,
    aes(x = Axis.1, y = Axis.2, color = Trees)
  ) +
    geom_point(size = 3.5, alpha = 0.9) +
    facet_wrap(~ Health_binary) +
    labs(
      title = paste0("16S Bray-Curtis PCoA — ", month_name, " · ", compartment_name),
      subtitle = "Colour = tree identity; faceted by health status",
      x = pcoa_info$xlab,
      y = pcoa_info$ylab,
      color = "Tree"
    ) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold")
    )
  
  prefix <- paste0(
    "16S_",
    month_name,
    "_",
    gsub(" ", "_", compartment_name)
  )
  
  ggsave(
    file.path(outdir_treecheck_figs, paste0(prefix, "_PCoA_health_shape_tree.png")),
    plot = p_health_tree,
    width = 7.5,
    height = 6,
    dpi = 600
  )
  
  ggsave(
    file.path(outdir_treecheck_figs, paste0(prefix, "_PCoA_tree_facet_health.png")),
    plot = p_tree_facet,
    width = 8,
    height = 5.5,
    dpi = 600
  )
  
  list(
    month = month_name,
    compartment = compartment_name,
    ps = ps_sub,
    metadata = md_sub,
    distance = D_bray,
    pcoa_df = pcoa_df,
    p_health_tree = p_health_tree,
    p_tree_facet = p_tree_facet
  )
}

months_treecheck <- c("July", "October")
compartments_treecheck <- c("Rhizosphere", "Endosphere", "Phyllosphere")

treecheck_pcoa_results <- list()

for (mo in months_treecheck) {
  for (co in compartments_treecheck) {
    key <- paste(mo, co, sep = "_")
    treecheck_pcoa_results[[key]] <- run_treecheck_pcoa_one(ps0, mo, co)
  }
}

############################################################
# 4) TEST TREE EFFECT WITHIN EACH MONTH x COMPARTMENT x HEALTH_STATUS
############################################################

test_tree_effect_one_group <- function(ps_in, month_name, compartment_name, health_status_name) {
  
  md <- as(sample_data(ps_in), "data.frame")
  
  keep_ids <- rownames(md)[
    as.character(md$Month) == month_name &
      as.character(md$Compartment) == compartment_name &
      as.character(md$Health_status) == health_status_name
  ]
  
  ps_sub <- prune_samples(keep_ids, ps_in)
  ps_sub <- prune_samples(sample_sums(ps_sub) > 0, ps_sub)
  ps_sub <- prune_taxa(taxa_sums(ps_sub) > 0, ps_sub)
  
  md_sub <- as(sample_data(ps_sub), "data.frame")
  md_sub$SampleID <- rownames(md_sub)
  md_sub$Trees <- factor(as.character(md_sub$Trees))
  
  if (nsamples(ps_sub) < 4 || length(unique(md_sub$Trees)) < 2) {
    return(tibble(
      Month = month_name,
      Compartment = compartment_name,
      Health_status = health_status_name,
      n_samples = nsamples(ps_sub),
      n_trees = length(unique(md_sub$Trees)),
      PERMANOVA_R2_Tree = NA_real_,
      PERMANOVA_F_Tree = NA_real_,
      PERMANOVA_p_Tree = NA_real_,
      Comment = "Not enough samples or tree levels"
    ))
  }
  
  D_bray <- make_bray_distance(ps_sub)
  ad <- vegan::adonis2(D_bray ~ Trees, data = md_sub, permutations = 999)
  
  tibble(
    Month = month_name,
    Compartment = compartment_name,
    Health_status = health_status_name,
    n_samples = nsamples(ps_sub),
    n_trees = length(unique(md_sub$Trees)),
    PERMANOVA_R2_Tree = as.numeric(ad$R2[1]),
    PERMANOVA_F_Tree = as.numeric(ad$F[1]),
    PERMANOVA_p_Tree = as.numeric(ad$`Pr(>F)`[1]),
    Comment = case_when(
      is.na(as.numeric(ad$`Pr(>F)`[1])) ~ "Not tested",
      as.numeric(ad$`Pr(>F)`[1]) < 0.05 ~ "Tree effect detected within this group",
      TRUE ~ "No significant tree effect within this group"
    )
  )
}

health_levels_treecheck <- c(
  "Healthy Rhizosphere", "Diseased Rhizosphere",
  "Healthy Endosphere", "Diseased Endosphere",
  "Healthy Phyllosphere", "Diseased Phyllosphere"
)

tree_effect_results <- purrr::map_dfr(
  months_treecheck,
  function(mo) {
    purrr::map_dfr(
      compartments_treecheck,
      function(co) {
        relevant_health <- health_levels_treecheck[
          grepl(co, health_levels_treecheck)
        ]
        purrr::map_dfr(
          relevant_health,
          function(hs) test_tree_effect_one_group(ps0, mo, co, hs)
        )
      }
    )
  }
)

write_csv(
  tree_effect_results,
  file.path(outdir_treecheck_tables, "01_tree_effect_within_month_compartment_health.csv")
)

############################################################
# 5) WITHIN-TREE VS BETWEEN-TREE DISTANCES
############################################################

compare_within_between_tree_distances_one_group <- function(ps_in, month_name, compartment_name, health_status_name) {
  
  md <- as(sample_data(ps_in), "data.frame")
  
  keep_ids <- rownames(md)[
    as.character(md$Month) == month_name &
      as.character(md$Compartment) == compartment_name &
      as.character(md$Health_status) == health_status_name
  ]
  
  ps_sub <- prune_samples(keep_ids, ps_in)
  ps_sub <- prune_samples(sample_sums(ps_sub) > 0, ps_sub)
  ps_sub <- prune_taxa(taxa_sums(ps_sub) > 0, ps_sub)
  
  md_sub <- as(sample_data(ps_sub), "data.frame")
  md_sub$SampleID <- rownames(md_sub)
  md_sub$Trees <- factor(as.character(md_sub$Trees))
  
  if (nsamples(ps_sub) < 4 || length(unique(md_sub$Trees)) < 2) {
    return(tibble(
      Month = month_name,
      Compartment = compartment_name,
      Health_status = health_status_name,
      Pair_type = NA_character_,
      n_pairs = NA_integer_,
      mean_Bray_distance = NA_real_,
      median_Bray_distance = NA_real_,
      sd_Bray_distance = NA_real_,
      Wilcoxon_p_within_vs_between = NA_real_,
      Comment = "Not enough samples or tree levels"
    ))
  }
  
  D_bray <- make_bray_distance(ps_sub)
  Dmat <- as.matrix(D_bray)
  
  sample_ids <- rownames(md_sub)
  pairs <- combn(sample_ids, 2, simplify = FALSE)
  
  pair_df <- purrr::map_dfr(
    pairs,
    function(pair_i) {
      s1 <- pair_i[1]
      s2 <- pair_i[2]
      
      tree1 <- as.character(md_sub$Trees[md_sub$SampleID == s1])
      tree2 <- as.character(md_sub$Trees[md_sub$SampleID == s2])
      
      tibble(
        Sample1 = s1,
        Sample2 = s2,
        Tree1 = tree1,
        Tree2 = tree2,
        Bray_distance = Dmat[s1, s2],
        Pair_type = ifelse(tree1 == tree2, "Within_tree", "Between_trees_same_group")
      )
    }
  )
  
  pval <- tryCatch(
    wilcox.test(Bray_distance ~ Pair_type, data = pair_df, exact = FALSE)$p.value,
    error = function(e) NA_real_
  )
  
  pair_df %>%
    group_by(Pair_type) %>%
    summarise(
      Month = month_name,
      Compartment = compartment_name,
      Health_status = health_status_name,
      n_pairs = n(),
      mean_Bray_distance = mean(Bray_distance, na.rm = TRUE),
      median_Bray_distance = median(Bray_distance, na.rm = TRUE),
      sd_Bray_distance = sd(Bray_distance, na.rm = TRUE),
      Wilcoxon_p_within_vs_between = pval,
      Comment = case_when(
        is.na(pval) ~ "Test not available",
        pval < 0.05 ~ "Within-tree and between-tree distances differ",
        TRUE ~ "No significant distance difference"
      ),
      .groups = "drop"
    ) %>%
    select(
      Month,
      Compartment,
      Health_status,
      Pair_type,
      n_pairs,
      mean_Bray_distance,
      median_Bray_distance,
      sd_Bray_distance,
      Wilcoxon_p_within_vs_between,
      Comment
    )
}

within_between_distance_results <- purrr::map_dfr(
  months_treecheck,
  function(mo) {
    purrr::map_dfr(
      compartments_treecheck,
      function(co) {
        relevant_health <- health_levels_treecheck[
          grepl(co, health_levels_treecheck)
        ]
        purrr::map_dfr(
          relevant_health,
          function(hs) compare_within_between_tree_distances_one_group(ps0, mo, co, hs)
        )
      }
    )
  }
)

write_csv(
  within_between_distance_results,
  file.path(outdir_treecheck_tables, "02_within_tree_vs_between_tree_distances.csv")
)

############################################################
# 6) TREE-LEVEL AGGREGATION ANALYSIS
#
# Este análisis suma las cuentas de las submuestras del mismo árbol
# dentro de Month x Compartment x Health_status x Trees.
# Así cada árbol pasa a ser una unidad.
############################################################

aggregate_phyloseq_by_tree <- function(ps_in) {
  
  otu <- get_otu_samples_x_taxa(ps_in)
  md <- as(sample_data(ps_in), "data.frame")
  md$SampleID <- rownames(md)
  
  md$TreeGroup <- paste(
    md$Month,
    md$Compartment,
    md$Health_status,
    md$Trees,
    sep = "__"
  )
  
  # Agregar conteos por grupo de árbol
  otu_agg_samples_x_taxa <- rowsum(otu, group = md$TreeGroup, reorder = FALSE)
  
  md_agg <- md %>%
    group_by(TreeGroup) %>%
    summarise(
      Month = first(Month),
      Compartment = first(Compartment),
      Health_status = first(Health_status),
      Trees = first(Trees),
      n_subsamples = n(),
      .groups = "drop"
    ) %>%
    as.data.frame()
  
  rownames(md_agg) <- md_agg$TreeGroup
  
  otu_agg_taxa_x_samples <- t(otu_agg_samples_x_taxa)
  
  ps_agg <- phyloseq(
    otu_table(otu_agg_taxa_x_samples, taxa_are_rows = TRUE),
    sample_data(md_agg),
    tax_table(ps_in)
  )
  
  if (!is.null(phy_tree(ps_in, errorIfNULL = FALSE))) {
    phy_tree(ps_agg) <- phy_tree(ps_in)
  }
  
  ps_agg <- prune_samples(sample_sums(ps_agg) > 0, ps_agg)
  ps_agg <- prune_taxa(taxa_sums(ps_agg) > 0, ps_agg)
  
  ps_agg
}

ps0_tree_agg <- aggregate_phyloseq_by_tree(ps0)

saveRDS(
  ps0_tree_agg,
  file.path(outdir_treecheck, "Phyloseq_16S_aggregated_by_tree.rds")
)

############################################################
# 7) PCoA AND PERMANOVA AT TREE-LEVEL
############################################################

run_treelevel_beta_one <- function(ps_in, month_name, compartment_name) {
  
  md <- as(sample_data(ps_in), "data.frame")
  
  keep_ids <- rownames(md)[
    as.character(md$Month) == month_name &
      as.character(md$Compartment) == compartment_name
  ]
  
  ps_sub <- prune_samples(keep_ids, ps_in)
  ps_sub <- prune_samples(sample_sums(ps_sub) > 0, ps_sub)
  ps_sub <- prune_taxa(taxa_sums(ps_sub) > 0, ps_sub)
  
  md_sub <- as(sample_data(ps_sub), "data.frame")
  md_sub$SampleID <- rownames(md_sub)
  md_sub$Health_binary <- factor(
    make_health_binary_from_health_status(md_sub$Health_status),
    levels = c("Healthy", "Diseased")
  )
  
  if (nsamples(ps_sub) < 4 || length(unique(md_sub$Health_binary)) < 2) {
    return(tibble(
      Month = month_name,
      Compartment = compartment_name,
      Analysis_level = "Tree-level aggregated",
      n_tree_units = nsamples(ps_sub),
      n_Healthy = sum(md_sub$Health_binary == "Healthy", na.rm = TRUE),
      n_Diseased = sum(md_sub$Health_binary == "Diseased", na.rm = TRUE),
      Bray_PERMANOVA_R2_Health = NA_real_,
      Bray_PERMANOVA_F_Health = NA_real_,
      Bray_PERMANOVA_p_Health = NA_real_,
      Comment = "Not enough tree-level units"
    ))
  }
  
  D_bray <- make_bray_distance(ps_sub)
  pcoa_info <- pcoa_labels_from_dist(D_bray)
  pcoa_obj <- pcoa_info$pcoa
  
  pcoa_df <- as.data.frame(pcoa_obj$vectors[, 1:2, drop = FALSE])
  pcoa_df$SampleID <- rownames(pcoa_df)
  pcoa_df <- left_join(pcoa_df, md_sub, by = "SampleID")
  
  p_treelevel <- ggplot(
    pcoa_df,
    aes(x = Axis.1, y = Axis.2, color = Health_binary, label = Trees)
  ) +
    geom_point(size = 4, alpha = 0.9) +
    geom_text(vjust = -0.9, size = 3.5, show.legend = FALSE) +
    labs(
      title = paste0("16S tree-level Bray-Curtis PCoA — ", month_name, " · ", compartment_name),
      subtitle = "Subsamples aggregated by tree",
      x = pcoa_info$xlab,
      y = pcoa_info$ylab,
      color = "Health status"
    ) +
    scale_color_manual(values = c("Healthy" = "#4B9EA7", "Diseased" = "#801723")) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold")
    )
  
  prefix <- paste0(
    "16S_",
    month_name,
    "_",
    gsub(" ", "_", compartment_name),
    "_tree_level"
  )
  
  ggsave(
    file.path(outdir_treecheck_figs, paste0(prefix, "_PCoA_Bray.png")),
    plot = p_treelevel,
    width = 7,
    height = 6,
    dpi = 600
  )
  
  ad <- vegan::adonis2(D_bray ~ Health_binary, data = md_sub, permutations = 999)
  
  tibble(
    Month = month_name,
    Compartment = compartment_name,
    Analysis_level = "Tree-level aggregated",
    n_tree_units = nsamples(ps_sub),
    n_Healthy = sum(md_sub$Health_binary == "Healthy", na.rm = TRUE),
    n_Diseased = sum(md_sub$Health_binary == "Diseased", na.rm = TRUE),
    Bray_PERMANOVA_R2_Health = as.numeric(ad$R2[1]),
    Bray_PERMANOVA_F_Health = as.numeric(ad$F[1]),
    Bray_PERMANOVA_p_Health = as.numeric(ad$`Pr(>F)`[1]),
    Comment = "Health effect tested after aggregating subsamples by tree"
  )
}

treelevel_health_results <- purrr::map_dfr(
  months_treecheck,
  function(mo) {
    purrr::map_dfr(
      compartments_treecheck,
      function(co) run_treelevel_beta_one(ps0_tree_agg, mo, co)
    )
  }
)

write_csv(
  treelevel_health_results,
  file.path(outdir_treecheck_tables, "03_tree_level_aggregated_health_effect.csv")
)

############################################################
# 8) PRINT SUMMARY TO CONSOLE
############################################################

cat("\n====================================================\n")
cat("TREE EFFECT WITHIN MONTH x COMPARTMENT x HEALTH_STATUS\n")
cat("====================================================\n")
print(tree_effect_results, n = Inf, width = Inf)

cat("\n====================================================\n")
cat("WITHIN-TREE VS BETWEEN-TREE DISTANCES\n")
cat("====================================================\n")
print(within_between_distance_results, n = Inf, width = Inf)

cat("\n====================================================\n")
cat("TREE-LEVEL AGGREGATED HEALTH EFFECT\n")
cat("====================================================\n")
print(treelevel_health_results, n = Inf, width = Inf)

cat("\n====================================================\n")
cat("FILES WRITTEN TO:\n")
cat("====================================================\n")
print(normalizePath(outdir_treecheck))
print(normalizePath(outdir_treecheck_figs))
print(normalizePath(outdir_treecheck_tables))






setwd("path/to/QIIME2_files_to_R_16S")

library(readr)
library(dplyr)

beta16S_comp <- read_csv(
  "02_beta_diversity_16S/tables/02_beta_by_compartment_summary_July_October.csv",
  show_col_types = FALSE
)

cat("\nColumn names:\n")
print(colnames(beta16S_comp))

beta16S_comp %>%
  select(
    Month,
    Compartment,
    n_samples,
    n_Healthy,
    n_Diseased,
    contains("PERMANOVA"),
    contains("betadisper")
  ) %>%
  print(n = Inf, width = Inf)








########### importante para el subsampling de arboles
############# Exploratorio
# ============================================================
# TREE-LEVEL ROBUSTNESS ANALYSIS - 16S
# Merge 4 subsamples per tree into one tree-level profile
# ============================================================

setwd("path/to/QIIME2_files_to_R_16S")

library(phyloseq)
library(vegan)
library(dplyr)
library(tibble)
library(readr)
library(ggplot2)
library(ape)

# ------------------------------------------------------------
# 1. Output folders
# ------------------------------------------------------------

OUTDIR <- "02C_tree_level_aggregation_16S"
FIGDIR <- file.path(OUTDIR, "figures")
TABDIR <- file.path(OUTDIR, "tables")

dir.create(OUTDIR, showWarnings = FALSE)
dir.create(FIGDIR, showWarnings = FALSE)
dir.create(TABDIR, showWarnings = FALSE)

# ------------------------------------------------------------
# 2. Load phyloseq object
# ------------------------------------------------------------

ps0 <- readRDS("Phyloseq_16S.rds")

cat("\nLoaded phyloseq object:\n")
print(ps0)

md0 <- as(sample_data(ps0), "data.frame")
md0$SampleID <- rownames(md0)

cat("\nMetadata columns:\n")
print(colnames(md0))

# ------------------------------------------------------------
# 3. Define health variable
# ------------------------------------------------------------

if ("Health_status_2" %in% colnames(md0)) {
  md0$Health_binary <- as.character(md0$Health_status_2)
} else if ("Health_status" %in% colnames(md0)) {
  md0$Health_binary <- case_when(
    grepl("Healthy", as.character(md0$Health_status), ignore.case = TRUE) ~ "Healthy",
    grepl("Diseased", as.character(md0$Health_status), ignore.case = TRUE) ~ "Diseased",
    TRUE ~ NA_character_
  )
} else {
  stop("No Health_status_2 or Health_status column found.")
}

required_cols <- c("Month", "Compartment", "Trees", "Health_binary")
missing_cols <- setdiff(required_cols, colnames(md0))

if (length(missing_cols) > 0) {
  stop(paste("Missing required metadata columns:", paste(missing_cols, collapse = ", ")))
}

# Keep only plant-associated compartments and Healthy/Diseased samples
md0 <- md0 %>%
  filter(
    Compartment %in% c("Rhizosphere", "Endosphere", "Phyllosphere"),
    Health_binary %in% c("Healthy", "Diseased")
  )

ps0 <- prune_samples(md0$SampleID, ps0)

cat("\nSamples retained for tree-level aggregation:\n")
print(table(md0$Month, md0$Compartment, md0$Health_binary))

# ------------------------------------------------------------
# 4. Create tree-level group ID
# ------------------------------------------------------------

md0 <- md0 %>%
  mutate(
    Tree_level_ID = paste(Month, Compartment, Health_binary, Trees, sep = "__")
  )

md0_for_phyloseq <- as.data.frame(md0)
rownames(md0_for_phyloseq) <- md0_for_phyloseq$SampleID
sample_data(ps0) <- sample_data(md0_for_phyloseq)

group_counts <- md0 %>%
  count(Month, Compartment, Health_binary, Trees, Tree_level_ID, name = "n_subsamples") %>%
  arrange(Month, Compartment, Health_binary, Trees)

write_csv(
  group_counts,
  file.path(TABDIR, "01_original_subsamples_per_tree_level_group_16S.csv")
)

cat("\nSummary of n_subsamples per tree-level group:\n")
utils::write.table(
  as.data.frame(table(as.data.frame(group_counts)$n_subsamples)),
  row.names = FALSE,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

cat("\nFirst rows of group_counts:\n")
utils::write.table(
  head(as.data.frame(group_counts), 30),
  row.names = FALSE,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

# ------------------------------------------------------------
# 5. Aggregate ASV counts by tree-level group - fast version
# ------------------------------------------------------------

otu <- as(otu_table(ps0), "matrix")

# Ensure samples are rows and ASVs are columns
if (taxa_are_rows(ps0)) {
  otu_samples <- t(otu)
} else {
  otu_samples <- otu
}

md_for_otu <- md0[match(rownames(otu_samples), md0$SampleID), ]

if (any(is.na(md_for_otu$Tree_level_ID))) {
  stop("Some samples in the OTU table do not have Tree_level_ID.")
}

otu_tree <- rowsum(
  otu_samples,
  group = md_for_otu$Tree_level_ID,
  reorder = FALSE
)

otu_tree <- otu_tree[, colSums(otu_tree) > 0, drop = FALSE]

otu_tree_phy <- otu_table(t(otu_tree), taxa_are_rows = TRUE)

cat("\nOriginal sample-level OTU table:\n")
cat(nrow(otu_samples), "samples x", ncol(otu_samples), "ASVs\n")

cat("\nAggregated tree-level OTU table:\n")
cat(nrow(otu_tree), "tree-level samples x", ncol(otu_tree), "ASVs\n")

# ------------------------------------------------------------
# 6. Create tree-level metadata and phyloseq object
# ------------------------------------------------------------

md_tree <- md0 %>%
  group_by(Tree_level_ID) %>%
  summarise(
    Month = first(Month),
    Compartment = first(Compartment),
    Health_binary = first(Health_binary),
    Trees = first(Trees),
    n_subsamples = n(),
    original_samples = paste(SampleID, collapse = "; "),
    .groups = "drop"
  ) %>%
  column_to_rownames("Tree_level_ID")

sample_data_tree <- sample_data(md_tree)

if (!is.null(tax_table(ps0, errorIfNULL = FALSE))) {
  tax_tree <- tax_table(ps0)[taxa_names(otu_tree_phy), ]
  ps_tree <- phyloseq(otu_tree_phy, tax_tree, sample_data_tree)
} else {
  ps_tree <- phyloseq(otu_tree_phy, sample_data_tree)
}

saveRDS(ps_tree, file.path(OUTDIR, "Phyloseq_16S_tree_level_aggregated.rds"))

write_csv(
  md_tree %>% rownames_to_column("Tree_level_ID"),
  file.path(TABDIR, "02_tree_level_metadata_16S.csv")
)

cat("\nTree-level phyloseq object:\n")
print(ps_tree)

cat("\nTree-level sample counts:\n")
print(table(sample_data(ps_tree)$Month, sample_data(ps_tree)$Compartment, sample_data(ps_tree)$Health_binary))

# ------------------------------------------------------------
# 7. Function: Bray PCoA + PERMANOVA + betadisper
# ------------------------------------------------------------

run_tree_level_beta <- function(ps_input, compartment_name) {
  
  md_all <- as(sample_data(ps_input), "data.frame")
  keep_ids <- rownames(md_all)[md_all$Compartment == compartment_name]
  
  ps_comp <- prune_samples(keep_ids, ps_input)
  ps_comp <- prune_taxa(taxa_sums(ps_comp) > 0, ps_comp)
  
  md <- as(sample_data(ps_comp), "data.frame")
  md$Tree_level_ID <- rownames(md)
  
  otu_mat <- as(otu_table(ps_comp), "matrix")
  if (taxa_are_rows(ps_comp)) {
    otu_mat <- t(otu_mat)
  }
  
  keep_samples <- rowSums(otu_mat) > 0
  otu_mat <- otu_mat[keep_samples, , drop = FALSE]
  md <- md[rownames(otu_mat), , drop = FALSE]
  
  dist_bray <- vegan::vegdist(otu_mat, method = "bray")
  
  # PCoA
  pcoa_res <- ape::pcoa(dist_bray)
  
  coords <- as.data.frame(pcoa_res$vectors[, 1:2])
  colnames(coords) <- c("PCoA1", "PCoA2")
  coords$Tree_level_ID <- rownames(coords)
  
  md_join <- md %>%
    select(Tree_level_ID, Month, Compartment, Health_binary, Trees, n_subsamples)
  
  coords <- coords %>%
    left_join(md_join, by = "Tree_level_ID")
  
  axis1 <- round(pcoa_res$values$Relative_eig[1] * 100, 1)
  axis2 <- round(pcoa_res$values$Relative_eig[2] * 100, 1)
  
  coords$PCoA1_var <- axis1
  coords$PCoA2_var <- axis2
  
  write_csv(
    coords,
    file.path(TABDIR, paste0("03_PCoA_coordinates_tree_level_", compartment_name, "_16S.csv"))
  )
  
  # PCoA without ellipses
  p_no_ellipse <- ggplot(coords, aes(x = PCoA1, y = PCoA2, color = Health_binary, shape = Month)) +
    geom_point(size = 4, alpha = 0.9) +
    geom_text(aes(label = Trees), vjust = -0.8, size = 3, show.legend = FALSE) +
    labs(
      title = paste0("Tree-level Bray-Curtis PCoA - ", compartment_name, " - 16S"),
      x = paste0("PCoA1 (", axis1, "%)"),
      y = paste0("PCoA2 (", axis2, "%)"),
      color = "Health status",
      shape = "Month"
    ) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5)
    )
  
  ggsave(
    file.path(FIGDIR, paste0("PCoA_tree_level_Bray_", compartment_name, "_16S.pdf")),
    p_no_ellipse,
    width = 7,
    height = 5
  )
  
  ggsave(
    file.path(FIGDIR, paste0("PCoA_tree_level_Bray_", compartment_name, "_16S.png")),
    p_no_ellipse,
    width = 7,
    height = 5,
    dpi = 300
  )
  
  # PCoA with health-status ellipses
  p_ellipse <- ggplot(coords, aes(x = PCoA1, y = PCoA2, color = Health_binary, shape = Month)) +
    geom_point(size = 4, alpha = 0.9) +
    geom_text(aes(label = Trees), vjust = -0.8, size = 3, show.legend = FALSE) +
    stat_ellipse(
      aes(group = Health_binary),
      type = "t",
      linetype = "dashed",
      linewidth = 0.6,
      show.legend = FALSE
    ) +
    labs(
      title = paste0("Tree-level Bray-Curtis PCoA with health ellipses - ", compartment_name, " - 16S"),
      x = paste0("PCoA1 (", axis1, "%)"),
      y = paste0("PCoA2 (", axis2, "%)"),
      color = "Health status",
      shape = "Month"
    ) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5)
    )
  
  ggsave(
    file.path(FIGDIR, paste0("PCoA_tree_level_Bray_", compartment_name, "_16S_with_health_ellipses.pdf")),
    p_ellipse,
    width = 7,
    height = 5
  )
  
  ggsave(
    file.path(FIGDIR, paste0("PCoA_tree_level_Bray_", compartment_name, "_16S_with_health_ellipses.png")),
    p_ellipse,
    width = 7,
    height = 5,
    dpi = 300
  )
  
  # Combined PERMANOVA: global model
  perm_combined <- adonis2(dist_bray ~ Month + Health_binary, data = md, permutations = 999)
  
  perm_combined_df <- as.data.frame(perm_combined) %>%
    rownames_to_column("term") %>%
    mutate(
      Compartment = compartment_name,
      Analysis = "Combined_July_October",
      Distance = "Bray"
    ) %>%
    relocate(Analysis, Compartment, Distance, term)
  
  # Combined PERMANOVA: term-wise model
  perm_terms <- adonis2(
    dist_bray ~ Month + Health_binary,
    data = md,
    permutations = 999,
    by = "terms"
  )
  
  perm_terms_df <- as.data.frame(perm_terms) %>%
    rownames_to_column("term") %>%
    mutate(
      Compartment = compartment_name,
      Analysis = "Tree_level_combined_by_terms",
      Distance = "Bray"
    ) %>%
    relocate(Analysis, Compartment, Distance, term)
  
  # Combined PERMANOVA with interaction
  perm_interaction <- adonis2(dist_bray ~ Month * Health_binary, data = md, permutations = 999)
  
  perm_interaction_df <- as.data.frame(perm_interaction) %>%
    rownames_to_column("term") %>%
    mutate(
      Compartment = compartment_name,
      Analysis = "Combined_July_October_interaction",
      Distance = "Bray"
    ) %>%
    relocate(Analysis, Compartment, Distance, term)
  
  # PERMANOVA and betadisper by month
  by_month_list <- list()
  betadisper_list <- list()
  
  for (m in sort(unique(md$Month))) {
    
    sample_ids_m <- rownames(md)[md$Month == m]
    
    otu_m <- otu_mat[sample_ids_m, , drop = FALSE]
    md_m <- md[sample_ids_m, , drop = FALSE]
    
    if (length(unique(md_m$Health_binary)) < 2) next
    
    dist_m <- vegan::vegdist(otu_m, method = "bray")
    
    perm_m <- adonis2(dist_m ~ Health_binary, data = md_m, permutations = 999)
    
    perm_m_df <- as.data.frame(perm_m) %>%
      rownames_to_column("term") %>%
      mutate(
        Analysis = "By_month",
        Month = m,
        Compartment = compartment_name,
        Distance = "Bray"
      ) %>%
      relocate(Analysis, Month, Compartment, Distance, term)
    
    by_month_list[[m]] <- perm_m_df
    
    bd <- betadisper(dist_m, md_m$Health_binary)
    bd_perm <- permutest(bd, permutations = 999)
    
    bd_dist <- data.frame(
      Tree_level_ID = names(bd$distances),
      distance_to_centroid = as.numeric(bd$distances)
    ) %>%
      left_join(
        md_m %>%
          select(Tree_level_ID, Month, Compartment, Health_binary, Trees),
        by = "Tree_level_ID"
      )
    
    bd_summary <- bd_dist %>%
      group_by(Month, Compartment, Health_binary) %>%
      summarise(
        n = n(),
        mean_distance_to_centroid = mean(distance_to_centroid),
        sd_distance_to_centroid = sd(distance_to_centroid),
        .groups = "drop"
      )
    
    bd_result <- data.frame(
      Month = m,
      Compartment = compartment_name,
      Distance = "Bray",
      betadisper_F = bd_perm$tab[1, "F"],
      betadisper_p = bd_perm$tab[1, "Pr(>F)"]
    )
    
    betadisper_list[[m]] <- list(
      result = bd_result,
      distances = bd_dist,
      summary = bd_summary
    )
  }
  
  return(list(
    pcoa = coords,
    permanova_combined = perm_combined_df,
    permanova_by_terms = perm_terms_df,
    permanova_interaction = perm_interaction_df,
    permanova_by_month = bind_rows(by_month_list),
    betadisper_results = bind_rows(lapply(betadisper_list, `[[`, "result")),
    betadisper_distances = bind_rows(lapply(betadisper_list, `[[`, "distances")),
    betadisper_summary = bind_rows(lapply(betadisper_list, `[[`, "summary"))
  ))
}

# ------------------------------------------------------------
# 8. Run analyses by compartment
# ------------------------------------------------------------

compartments <- c("Rhizosphere", "Endosphere", "Phyllosphere")

results <- lapply(compartments, function(x) run_tree_level_beta(ps_tree, x))
names(results) <- compartments

permanova_combined_all <- bind_rows(lapply(results, `[[`, "permanova_combined"))
permanova_by_terms_all <- bind_rows(lapply(results, `[[`, "permanova_by_terms"))
permanova_interaction_all <- bind_rows(lapply(results, `[[`, "permanova_interaction"))
permanova_by_month_all <- bind_rows(lapply(results, `[[`, "permanova_by_month"))
betadisper_results_all <- bind_rows(lapply(results, `[[`, "betadisper_results"))
betadisper_distances_all <- bind_rows(lapply(results, `[[`, "betadisper_distances"))
betadisper_summary_all <- bind_rows(lapply(results, `[[`, "betadisper_summary"))

# ------------------------------------------------------------
# 9. Save output tables
# ------------------------------------------------------------

write_csv(
  permanova_combined_all,
  file.path(TABDIR, "04_PERMANOVA_tree_level_combined_Month_plus_Health_16S.csv")
)

write_csv(
  permanova_interaction_all,
  file.path(TABDIR, "05_PERMANOVA_tree_level_combined_Month_times_Health_16S.csv")
)

write_csv(
  permanova_by_month_all,
  file.path(TABDIR, "06_PERMANOVA_tree_level_by_month_Health_16S.csv")
)

write_csv(
  betadisper_results_all,
  file.path(TABDIR, "07_betadisper_tree_level_by_month_Health_16S.csv")
)

write_csv(
  betadisper_distances_all,
  file.path(TABDIR, "08_betadisper_tree_level_distances_to_centroid_16S.csv")
)

write_csv(
  betadisper_summary_all,
  file.path(TABDIR, "09_betadisper_tree_level_summary_by_health_16S.csv")
)

write_csv(
  permanova_by_terms_all,
  file.path(TABDIR, "10_PERMANOVA_tree_level_combined_by_terms_16S.csv")
)

# ------------------------------------------------------------
# 10. Print key outputs robustly
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("TREE-LEVEL PERMANOVA: combined July + October\n")
cat("Model: Bray ~ Month + Health_binary\n")
cat("============================================================\n")
utils::write.table(
  as.data.frame(permanova_combined_all),
  row.names = FALSE,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

cat("\n============================================================\n")
cat("TREE-LEVEL PERMANOVA: combined by terms\n")
cat("Model: Bray ~ Month + Health_binary, by = 'terms'\n")
cat("============================================================\n")
utils::write.table(
  as.data.frame(permanova_by_terms_all),
  row.names = FALSE,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

cat("\n============================================================\n")
cat("TREE-LEVEL PERMANOVA: Month * Health_binary interaction\n")
cat("============================================================\n")
utils::write.table(
  as.data.frame(permanova_interaction_all),
  row.names = FALSE,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

cat("\n============================================================\n")
cat("TREE-LEVEL PERMANOVA: by month\n")
cat("Model: Bray ~ Health_binary\n")
cat("============================================================\n")
utils::write.table(
  as.data.frame(permanova_by_month_all),
  row.names = FALSE,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

cat("\n============================================================\n")
cat("TREE-LEVEL BETADISPER: by month\n")
cat("============================================================\n")
utils::write.table(
  as.data.frame(betadisper_results_all),
  row.names = FALSE,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

cat("\n============================================================\n")
cat("TREE-LEVEL DISTANCE TO CENTROID SUMMARY\n")
cat("============================================================\n")
utils::write.table(
  as.data.frame(betadisper_summary_all),
  row.names = FALSE,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

cat("\nDONE. Outputs saved in:\n")
cat(OUTDIR, "\n")
