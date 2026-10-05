
############################################################
# SECTION 1 — BASIC METRICS (ITS, OCTOBER ONLY)
#
# Objetivos:
# 1) Resumir profundidad de secuenciación por muestra y por grupo
# 2) Calcular ASVs detectados por compartimento y estado de salud
# 3) Calcular géneros detectados por compartimento y estado de salud
# 4) Exportar tablas limpias para escritura del manuscrito
############################################################

setwd("path/to/QIIME2_files_to_R_ITS")

suppressPackageStartupMessages({
  library(phyloseq)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(stringr)
})

############################################################
# 0) INPUTS AND OUTPUT
############################################################

physeq_path <- "Phyloseq_ITS.rds"

outdir_basic <- "00_basic_metrics_ITS"
dir.create(outdir_basic, recursive = TRUE, showWarnings = FALSE)

############################################################
# 1) LOAD PHYLOSEQ
############################################################

ps <- readRDS(physeq_path)
ps <- prune_samples(sample_sums(ps) > 0, ps)

sample_df <- as(sample_data(ps), "data.frame") %>%
  rownames_to_column("SampleID") %>%
  mutate(
    Reads = sample_sums(ps),
    Month = as.character(Month),
    Compartment = as.character(Compartment),
    Health_status = as.character(Health_status)
  )

cat("\n====================================================\n")
cat("PHYLOSEQ LOADED\n")
cat("====================================================\n")
print(ps)

############################################################
# 2) BASIC CHECKS
############################################################

cat("\n====================================================\n")
cat("BASIC SAMPLE CHECKS\n")
cat("====================================================\n")
cat("\nMonth table:\n")
print(table(sample_df$Month, useNA = "ifany"))

cat("\nCompartment table:\n")
print(table(sample_df$Compartment, useNA = "ifany"))

cat("\nHealth_status table:\n")
print(table(sample_df$Health_status, useNA = "ifany"))

############################################################
# 3) GLOBAL TOTALS
############################################################

global_totals <- tibble(
  n_samples = nsamples(ps),
  total_reads = sum(sample_sums(ps)),
  mean_reads = mean(sample_sums(ps)),
  median_reads = median(sample_sums(ps)),
  min_reads = min(sample_sums(ps)),
  max_reads = max(sample_sums(ps)),
  n_ASVs_total = ntaxa(ps)
)

write_csv(
  global_totals,
  file.path(outdir_basic, "01_global_totals_ITS.csv")
)

############################################################
# 4) HELPER FUNCTIONS
############################################################

count_detected_asvs <- function(sample_ids, otu_mat) {
  submat <- otu_mat[, sample_ids, drop = FALSE]
  sum(rowSums(submat > 0) > 0)
}

count_detected_genera <- function(sample_ids, otu_mat, tax_df) {
  submat <- otu_mat[, sample_ids, drop = FALSE]
  present_asvs <- rownames(submat)[rowSums(submat > 0) > 0]
  
  tax_df %>%
    filter(ASV_ID %in% present_asvs) %>%
    filter(!is.na(Genus), Genus != "") %>%
    summarise(n_genera = n_distinct(Genus)) %>%
    pull(n_genera)
}

############################################################
# 5) PREPARE OTU + TAX TABLES
############################################################

otu_mat <- as(otu_table(ps), "matrix")
if (!taxa_are_rows(ps)) {
  otu_mat <- t(otu_mat)
}

tax_df <- as.data.frame(tax_table(ps)) %>%
  rownames_to_column("ASV_ID") %>%
  mutate(
    Genus = as.character(Genus),
    Genus = ifelse(Genus == "" | is.na(Genus), NA, Genus)
  )

############################################################
# 6) SUMMARY BY COMPARTMENT
############################################################

summary_compartment <- sample_df %>%
  group_by(Compartment) %>%
  summarise(
    n_samples = n(),
    total_reads = sum(Reads, na.rm = TRUE),
    mean_reads = mean(Reads, na.rm = TRUE),
    median_reads = median(Reads, na.rm = TRUE),
    min_reads = min(Reads, na.rm = TRUE),
    max_reads = max(Reads, na.rm = TRUE),
    sd_reads = sd(Reads, na.rm = TRUE),
    sample_ids = list(SampleID),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    n_ASVs_detected = count_detected_asvs(sample_ids, otu_mat),
    n_genera_detected = count_detected_genera(sample_ids, otu_mat, tax_df)
  ) %>%
  ungroup() %>%
  select(-sample_ids) %>%
  arrange(Compartment)

write_csv(
  summary_compartment,
  file.path(outdir_basic, "02_summary_by_compartment_ITS.csv")
)

############################################################
# 7) SUMMARY BY COMPARTMENT + HEALTH STATUS
############################################################

summary_group <- sample_df %>%
  group_by(Compartment, Health_status) %>%
  summarise(
    n_samples = n(),
    total_reads = sum(Reads, na.rm = TRUE),
    mean_reads = mean(Reads, na.rm = TRUE),
    median_reads = median(Reads, na.rm = TRUE),
    min_reads = min(Reads, na.rm = TRUE),
    max_reads = max(Reads, na.rm = TRUE),
    sd_reads = sd(Reads, na.rm = TRUE),
    sample_ids = list(SampleID),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    n_ASVs_detected = count_detected_asvs(sample_ids, otu_mat),
    n_genera_detected = count_detected_genera(sample_ids, otu_mat, tax_df)
  ) %>%
  ungroup() %>%
  select(-sample_ids) %>%
  arrange(Compartment, Health_status)

write_csv(
  summary_group,
  file.path(outdir_basic, "03_summary_by_compartment_health_ITS.csv")
)

############################################################
# 8) PER-SAMPLE METRICS
############################################################

per_sample_basic_metrics <- sample_df %>%
  rowwise() %>%
  mutate(
    n_ASVs_detected = count_detected_asvs(SampleID, otu_mat),
    n_genera_detected = count_detected_genera(SampleID, otu_mat, tax_df)
  ) %>%
  ungroup() %>%
  select(
    SampleID,
    Month,
    Compartment,
    Health_status,
    Reads,
    n_ASVs_detected,
    n_genera_detected
  ) %>%
  arrange(Compartment, Health_status, SampleID)

write_csv(
  per_sample_basic_metrics,
  file.path(outdir_basic, "04_per_sample_basic_metrics_ITS.csv")
)

############################################################
# 9) PER-SAMPLE SUMMARY BY GROUP
############################################################

per_sample_summary_by_group <- per_sample_basic_metrics %>%
  group_by(Compartment, Health_status) %>%
  summarise(
    n_samples = n(),
    reads_mean = mean(Reads),
    reads_median = median(Reads),
    reads_min = min(Reads),
    reads_max = max(Reads),
    reads_sd = sd(Reads),
    ASVs_mean = mean(n_ASVs_detected),
    ASVs_median = median(n_ASVs_detected),
    ASVs_min = min(n_ASVs_detected),
    ASVs_max = max(n_ASVs_detected),
    ASVs_sd = sd(n_ASVs_detected),
    genera_mean = mean(n_genera_detected),
    genera_median = median(n_genera_detected),
    genera_min = min(n_genera_detected),
    genera_max = max(n_genera_detected),
    genera_sd = sd(n_genera_detected),
    .groups = "drop"
  ) %>%
  arrange(Compartment, Health_status)

write_csv(
  per_sample_summary_by_group,
  file.path(outdir_basic, "05_per_sample_summary_by_group_ITS.csv")
)

############################################################
# 10) OUTLIER SCREEN (Z-SCORE, DESCRIPTIVE ONLY)
############################################################

per_sample_outlier_flags <- per_sample_basic_metrics %>%
  group_by(Compartment, Health_status) %>%
  mutate(
    reads_z = ifelse(sd(Reads) == 0, 0, (Reads - mean(Reads)) / sd(Reads)),
    ASVs_z = ifelse(sd(n_ASVs_detected) == 0, 0, (n_ASVs_detected - mean(n_ASVs_detected)) / sd(n_ASVs_detected)),
    genera_z = ifelse(sd(n_genera_detected) == 0, 0, (n_genera_detected - mean(n_genera_detected)) / sd(n_genera_detected)),
    potential_outlier = ifelse(
      abs(reads_z) > 3 | abs(ASVs_z) > 3 | abs(genera_z) > 3,
      "YES",
      "NO"
    )
  ) %>%
  ungroup()

write_csv(
  per_sample_outlier_flags,
  file.path(outdir_basic, "06_per_sample_outlier_flags_ITS.csv")
)

############################################################
# 11) MANUSCRIPT-STYLE TABLE
############################################################

methods_table <- summary_group %>%
  mutate(
    reads_mean_range = paste0(
      round(mean_reads, 1), " (",
      min_reads, "–", max_reads, ")"
    )
  ) %>%
  select(
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

write_csv(
  methods_table,
  file.path(outdir_basic, "07_methods_results_ready_table_ITS.csv")
)

############################################################
# 12) CONSOLE OUTPUT
############################################################

cat("\n====================================================\n")
cat("GLOBAL TOTALS\n")
cat("====================================================\n")
print(global_totals)

cat("\n====================================================\n")
cat("SUMMARY BY COMPARTMENT\n")
cat("====================================================\n")
print(summary_compartment, n = Inf)

cat("\n====================================================\n")
cat("SUMMARY BY COMPARTMENT + HEALTH STATUS\n")
cat("====================================================\n")
print(summary_group, n = Inf)

cat("\n====================================================\n")
cat("PER-SAMPLE BASIC METRICS\n")
cat("====================================================\n")
print(as.data.frame(head(per_sample_basic_metrics, 20)))

cat("\n====================================================\n")
cat("PER-SAMPLE SUMMARY BY GROUP\n")
cat("====================================================\n")
print(per_sample_summary_by_group, n = Inf)

cat("\n====================================================\n")
cat("POTENTIAL OUTLIERS\n")
cat("====================================================\n")
print(
  per_sample_outlier_flags %>%
    filter(potential_outlier == "YES"),
  n = Inf
)

cat("\n====================================================\n")
cat("MANUSCRIPT-STYLE SENTENCES BY COMPARTMENT\n")
cat("====================================================\n")
for (i in seq_len(nrow(summary_compartment))) {
  row_i <- summary_compartment[i, ]
  cat(
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




###############################################################################################################################################################
###############################################################################################################################################################
###############################################################################################################################################################
###############################################################################################################################################################
###############################################################################################################################################################

############################################################
# SECTION 2 — ALPHA DIVERSITY (ITS, OCTOBER ONLY)
#
# Objetivos:
# 1) Leer alpha diversity desde el .qzv de rarefacción ITS
# 2) Unir metadata externo (metadata_ITS_analysis.tsv)
# 3) Usar solo las métricas realmente presentes en el .qzv
# 4) Resumir valores por compartimento y estado de salud
# 5) Hacer contrastes Healthy vs Diseased dentro de cada compartimento
# 6) Generar tablas limpias para Results / Supplementary
############################################################


setwd("path/to/QIIME2_files_to_R_ITS")

suppressPackageStartupMessages({
  library(tidyverse)
  library(multcompView)
  library(patchwork)
  library(readr)
})

############################################################
# 0) INPUTS AND OUTPUT
############################################################

qzv_file <- "alpha-rarefaction_ITS_maxdepth39000.qzv"
metadata_file <- "metadata_ITS_analysis.tsv"
target_depth <- 39000
month_to_do <- "October"

outdir_alpha <- "01_alpha_diversity_ITS"
dir.create(outdir_alpha, recursive = TRUE, showWarnings = FALSE)

hs_levels <- c(
  "Bulk soil",
  "Healthy Rhizosphere", "Diseased Rhizosphere",
  "Healthy Endosphere",  "Diseased Endosphere",
  "Healthy Phyllosphere","Diseased Phyllosphere"
)

hs_colors <- c(
  "Bulk soil"              = "#896054",
  "Healthy Rhizosphere"    = "#4B9EA7",
  "Diseased Rhizosphere"   = "#801723",
  "Healthy Endosphere"     = "#6B9863",
  "Diseased Endosphere"    = "#E7C313",
  "Healthy Phyllosphere"   = "#4A5263",
  "Diseased Phyllosphere"  = "#6B4EA7"
)

# Requested metrics
metrics_requested <- c(
  "Observed features"   = "observed_features\\.csv$",
  "Chao1"               = "chao1\\.csv$",
  "Shannon's Diversity" = "shannon\\.csv$"
)

############################################################
# 1) BASIC CHECKS
############################################################

cat("\n====================================================\n")
cat("WORKING DIRECTORY\n")
cat("====================================================\n")
print(getwd())

cat("\n====================================================\n")
cat("INPUT FILES EXIST?\n")
cat("====================================================\n")
print(file.exists(qzv_file))
print(file.exists(metadata_file))

if (!file.exists(qzv_file)) {
  stop("Missing qzv file: ", qzv_file)
}
if (!file.exists(metadata_file)) {
  stop("Missing metadata file: ", metadata_file)
}

zip_listing <- tryCatch(
  unzip(qzv_file, list = TRUE),
  error = function(e) e
)

cat("\n====================================================\n")
cat("ZIP LISTING TEST\n")
cat("====================================================\n")
print(zip_listing)

if (inherits(zip_listing, "error")) {
  stop("R cannot read the .qzv as a zip file.")
}

############################################################
# 2) DETECT AVAILABLE METRICS INSIDE QZV
############################################################

available_metric_names <- zip_listing$Name

metrics_available <- metrics_requested[
  vapply(
    metrics_requested,
    function(pat) any(grepl(pat, available_metric_names)),
    logical(1)
  )
]

metrics_missing <- setdiff(names(metrics_requested), names(metrics_available))

cat("\n====================================================\n")
cat("METRICS REQUESTED\n")
cat("====================================================\n")
print(names(metrics_requested))

cat("\n====================================================\n")
cat("METRICS AVAILABLE IN THIS QZV\n")
cat("====================================================\n")
print(names(metrics_available))

cat("\n====================================================\n")
cat("METRICS MISSING IN THIS QZV\n")
cat("====================================================\n")
print(metrics_missing)

if (length(metrics_available) == 0) {
  stop("No requested alpha metrics were found inside the qzv.")
}

############################################################
# 3) LOAD EXTERNAL METADATA — ROBUST
############################################################

metadata_external_raw <- read.delim(
  metadata_file,
  sep = "\t",
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

cat("\n====================================================\n")
cat("RAW EXTERNAL METADATA COLUMNS\n")
cat("====================================================\n")
print(colnames(metadata_external_raw))

possible_id_cols <- intersect(
  c("sample-id", "SampleID", "sample_id", "#SampleID"),
  colnames(metadata_external_raw)
)

if (length(possible_id_cols) == 0) {
  stop(
    paste0(
      "Could not identify the sample ID column in metadata file.\n",
      "Available columns are: ",
      paste(colnames(metadata_external_raw), collapse = ", ")
    )
  )
}

id_col <- possible_id_cols[1]

metadata_external <- metadata_external_raw %>%
  as_tibble() %>%
  mutate(sample_id = .data[[id_col]]) %>%
  mutate(
    sample_id = as.character(sample_id),
    Month = as.character(Month),
    Compartment = as.character(Compartment),
    Health_status = as.character(Health_status)
  )

cat("\n====================================================\n")
cat("EXTERNAL METADATA LOADED\n")
cat("====================================================\n")
print(dim(metadata_external))
print(colnames(metadata_external))
cat("\nUsing sample ID column:\n")
print(id_col)

############################################################
# 4) HELPERS
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

label_stacked <- function(x) gsub(" ", "\n", x)

sig_mark <- function(p) ifelse(
  is.na(p), NA_character_,
  ifelse(p < 0.001, "***",
         ifelse(p < 0.01,  "**",
                ifelse(p < 0.05,  "*", "ns")))
)

digits_by_metric <- function(metric) dplyr::case_when(
  grepl("Shannon", metric, ignore.case = TRUE)        ~ 2L,
  grepl("Observed", metric, ignore.case = TRUE)       ~ 0L,
  grepl("Chao1", metric, ignore.case = TRUE)          ~ 0L,
  TRUE ~ 2L
)

read_metric_csv <- function(qzv, filename_regex) {
  tmpdir <- file.path(tempdir(), "alpha_raref_tmp_its")
  
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
    stop("No file matching pattern found inside qzv: ", filename_regex)
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
# 5) BUILD DATA FOR ONE METRIC
############################################################

df_for_metric_october <- function(metric_regex) {
  wide <- read_metric_csv(qzv_file, metric_regex)
  
  cat("\nColumns detected in qzv metric table:\n")
  print(names(wide))
  
  if (!"sample_id" %in% names(wide)) {
    stop("Missing required qzv column: sample_id")
  }
  
  depth_iter_cols <- grep("^depth_[0-9]+_iter_[0-9]+$", names(wide), value = TRUE)
  
  if (length(depth_iter_cols) == 0) {
    stop("No depth_*_iter_* columns found.")
  }
  
  long <- wide %>%
    mutate(sample_id = as.character(sample_id)) %>%
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
    left_join(
      metadata_external %>%
        select(sample_id, Name, Month, Compartment, Health_status),
      by = "sample_id"
    )
  
  if (any(is.na(long$Month))) {
    missing_ids <- long %>%
      filter(is.na(Month)) %>%
      distinct(sample_id) %>%
      pull(sample_id)
    
    stop(
      paste0(
        "Some sample_id values from qzv could not be matched to metadata_ITS_analysis.tsv.\n",
        "Examples: ",
        paste(head(missing_ids, 10), collapse = ", ")
      )
    )
  }
  
  long <- long %>%
    filter(tolower(Month) == tolower(month_to_do))
  
  if (nrow(long) == 0) {
    stop("No rows found for October after metadata join.")
  }
  
  available_depths <- sort(unique(long$depth))
  depth_use <- available_depths[which.min(abs(available_depths - target_depth))]
  
  cat("\nAvailable depths for October:\n")
  print(available_depths)
  cat("Target depth:", target_depth, "\n")
  cat("Depth used:", depth_use, "\n")
  
  alpha <- long %>%
    filter(depth == depth_use) %>%
    group_by(sample_id, Name, Health_status, Compartment, Month) %>%
    summarise(alpha_value = mean(value, na.rm = TRUE), .groups = "drop") %>%
    mutate(
      Health_status = factor(Health_status, levels = hs_levels),
      alpha_value   = safe_num(alpha_value)
    ) %>%
    filter(Compartment %in% c("Bulk soil", "Rhizosphere", "Endosphere", "Phyllosphere")) %>%
    select(
      sample_id,
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

for (mtr in names(metrics_available)) {
  cat("\n========================================\n")
  cat("PROCESSING:", mtr, "\n")
  cat("========================================\n")
  
  tmp <- df_for_metric_october(metrics_available[[mtr]]) %>%
    mutate(Metric = mtr, .before = 2)
  
  alpha_values_list[[mtr]] <- tmp
}

alpha_values_long <- bind_rows(alpha_values_list) %>%
  relocate(Metric, Month, sample_id, Compartment, Health_status, alpha_value)

write_csv(alpha_values_long, file.path(outdir_alpha, "01_alpha_values_long_ITS.csv"))

############################################################
# 7) SUMMARY TABLE
############################################################

summarise_alpha <- function(df) {
  df %>%
    group_by(Health_status) %>%
    summarise(
      n       = dplyr::n(),
      mean    = mean(alpha_value, na.rm = TRUE),
      sd      = if (n() > 1) sd(alpha_value, na.rm = TRUE) else 0,
      se      = if (n() > 1) sd / sqrt(n()) else 0,
      tval    = if (n() > 1) qt(0.975, df = n() - 1) else 0,
      ci95_lo = mean - tval * se,
      ci95_hi = mean + tval * se,
      median  = median(alpha_value, na.rm = TRUE),
      q1      = quantile(alpha_value, 0.25, na.rm = TRUE),
      q3      = quantile(alpha_value, 0.75, na.rm = TRUE),
      min     = min(alpha_value, na.rm = TRUE),
      max     = max(alpha_value, na.rm = TRUE),
      .groups = "drop"
    )
}

alpha_summary_by_group <- alpha_values_long %>%
  group_by(Metric) %>%
  group_modify(~ summarise_alpha(.x)) %>%
  ungroup() %>%
  mutate(Health_status = factor(Health_status, levels = hs_levels)) %>%
  arrange(Metric, Health_status)

write_csv(
  alpha_summary_by_group,
  file.path(outdir_alpha, "02_alpha_summary_by_group_ITS.csv")
)

############################################################
# 8) ANOVA + TUKEY
############################################################

tidy_tukey <- function(df) {
  fit <- aov(alpha_value ~ Health_status, data = df)
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

pairwise_all_min <- alpha_values_long %>%
  group_by(Metric) %>%
  group_modify(~ tidy_tukey(.x)) %>%
  ungroup() %>%
  separate(
    contrast,
    into = c("group1", "group2"),
    sep = "-",
    remove = FALSE
  ) %>%
  mutate(
    Compartment_pair = case_when(
      str_detect(group1, "Rhizosphere")  & str_detect(group2, "Rhizosphere")  ~ "Rhizosphere",
      str_detect(group1, "Endosphere")   & str_detect(group2, "Endosphere")   ~ "Endosphere",
      str_detect(group1, "Phyllosphere") & str_detect(group2, "Phyllosphere") ~ "Phyllosphere",
      TRUE ~ "Mixed"
    ),
    dig = digits_by_metric(Metric),
    estimate  = round(estimate,  dig),
    conf.low  = round(conf.low,  dig),
    conf.high = round(conf.high, dig),
    sig = sig_mark(p_adj)
  ) %>%
  select(Metric, Compartment_pair, contrast, estimate, conf.low, conf.high, p_adj, sig) %>%
  arrange(
    Metric,
    factor(Compartment_pair, levels = c("Rhizosphere", "Endosphere", "Phyllosphere", "Mixed")),
    contrast
  )

write_csv(
  pairwise_all_min,
  file.path(outdir_alpha, "03_alpha_pairwise_Tukey_min_ITS.csv")
)

############################################################
# 9) WITHIN-COMPARTMENT HEALTHY VS DISEASED ONLY
############################################################

target_contrasts <- c(
  "Diseased Rhizosphere-Healthy Rhizosphere",
  "Diseased Endosphere-Healthy Endosphere",
  "Diseased Phyllosphere-Healthy Phyllosphere"
)

alpha_within_compartment_tests <- pairwise_all_min %>%
  filter(contrast %in% target_contrasts) %>%
  mutate(
    Month = "October",
    Compartment = case_when(
      str_detect(contrast, "Rhizosphere")  ~ "Rhizosphere",
      str_detect(contrast, "Endosphere")   ~ "Endosphere",
      str_detect(contrast, "Phyllosphere") ~ "Phyllosphere"
    )
  ) %>%
  select(Metric, Month, Compartment, contrast, estimate, conf.low, conf.high, p_adj, sig)

write_csv(
  alpha_within_compartment_tests,
  file.path(outdir_alpha, "04_alpha_within_compartment_tests_ITS.csv")
)

############################################################
# 10) WILCOX VALIDATION
############################################################

alpha_within_compartment_wilcox <- alpha_values_long %>%
  filter(Health_status != "Bulk soil") %>%
  mutate(
    Compartment = as.character(Compartment),
    Health_status = as.character(Health_status)
  ) %>%
  group_by(Metric, Compartment) %>%
  group_modify(~{
    df <- .x
    groups <- unique(df$Health_status)
    
    if (length(groups) != 2) {
      return(tibble(
        group1 = NA_character_,
        group2 = NA_character_,
        wilcox_p = NA_real_
      ))
    }
    
    wt <- wilcox.test(alpha_value ~ Health_status, data = df, exact = FALSE)
    
    tibble(
      group1 = groups[1],
      group2 = groups[2],
      wilcox_p = wt$p.value
    )
  }) %>%
  ungroup() %>%
  mutate(Month = "October") %>%
  select(Metric, Month, Compartment, group1, group2, wilcox_p) %>%
  arrange(Metric, Compartment)

write_csv(
  alpha_within_compartment_wilcox,
  file.path(outdir_alpha, "05_alpha_within_compartment_wilcox_ITS.csv")
)

############################################################
# 11) KEY CONTRASTS TABLE
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
  arrange(Compartment, Metric)

write_csv(
  alpha_key_contrasts,
  file.path(outdir_alpha, "06_alpha_key_contrasts_tukey_wilcox_ITS.csv")
)

alpha_key_contrasts_compact <- alpha_key_contrasts %>%
  mutate(
    estimate_round = round(estimate, digits_by_metric(Metric)),
    conf.low_round = round(conf.low, digits_by_metric(Metric)),
    conf.high_round = round(conf.high, digits_by_metric(Metric)),
    p_adj_round = signif(p_adj, 3),
    wilcox_p_round = signif(wilcox_p, 3),
    estimate_ci = paste0(
      estimate_round, " [", conf.low_round, ", ", conf.high_round, "]"
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
  arrange(Compartment, Metric)

write_csv(
  alpha_key_contrasts_compact,
  file.path(outdir_alpha, "07_alpha_key_contrasts_compact_for_paper_ITS.csv")
)

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
  arrange(Compartment, Metric)

write_csv(
  alpha_key_contrasts_writing,
  file.path(outdir_alpha, "08_alpha_key_contrasts_writing_table_ITS.csv")
)

############################################################
# 12) TABLE FOR WRITING
############################################################

alpha_table_for_writing <- alpha_summary_by_group %>%
  mutate(
    mean_se = paste0(round(mean, 2), " ± ", round(se, 2))
  ) %>%
  select(Metric, Health_status, n, mean, se, mean_se, median, ci95_lo, ci95_hi)

write_csv(
  alpha_table_for_writing,
  file.path(outdir_alpha, "09_alpha_table_for_writing_ITS.csv")
)

############################################################
# 13) PLOT FUNCTION
############################################################

plot_violin_with_SE_Tukey <- function(metric_name) {
  dfm <- alpha_values_long %>%
    filter(Metric == metric_name) %>%
    mutate(Health_status = factor(Health_status, levels = hs_levels))
  
  sumdf <- dfm %>%
    group_by(Health_status) %>%
    summarise(
      n    = n(),
      mean = mean(alpha_value, na.rm = TRUE),
      sd   = if (n() > 1) sd(alpha_value, na.rm = TRUE) else 0,
      se   = if (n() > 1) sd / sqrt(n()) else 0,
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
      subtitle = "October — violins by condition; mean ± SE; Tukey letters"
    ) +
    theme_bw() +
    theme(
      legend.position = "none",
      axis.title.y = element_text(size = 16, face = "bold"),
      axis.text.x = element_text(size = 12, lineheight = 0.9),
      panel.grid = element_blank(),
      plot.margin = margin(10, 10, 20, 10, "pt")
    )
}

available_metric_names_clean <- unique(alpha_values_long$Metric)

plot_list <- lapply(available_metric_names_clean, plot_violin_with_SE_Tukey)

p_oct <- wrap_plots(plot_list, ncol = 1) +
  patchwork::plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")")

print(p_oct)

ggsave(
  file.path(outdir_alpha, "10_alpha_panel_ITS.png"),
  plot = p_oct,
  width = 14, height = 5 * length(plot_list), dpi = 600
)

ggsave(
  file.path(outdir_alpha, "10_alpha_panel_ITS.tiff"),
  plot = p_oct,
  width = 14, height = 5 * length(plot_list), dpi = 600,
  compression = "lzw"
)




############################################################
# 14) CONSOLE OUTPUT
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
cat("ALPHA KEY CONTRASTS — WRITING TABLE\n")
cat("====================================================\n")
print(alpha_key_contrasts_writing, n = Inf, width = Inf)

cat("\n====================================================\n")
cat("FILES WRITTEN TO:\n")
cat("====================================================\n")
print(normalizePath(outdir_alpha))


############################################################
# FIX — FINAL TABLES FOR ITS ALPHA (EXPLICIT VERSION)
############################################################

alpha_key_contrasts_compact <- alpha_key_contrasts %>%
  mutate(
    estimate_round = round(estimate, digits_by_metric(Metric)),
    conf_low_round = round(conf.low, digits_by_metric(Metric)),
    conf_high_round = round(conf.high, digits_by_metric(Metric)),
    tukey_p_round = signif(p_adj, 3),
    wilcox_p_round = signif(wilcox_p, 3),
    Tukey_estimate_CI = paste0(
      estimate_round, " [", conf_low_round, ", ", conf_high_round, "]"
    )
  ) %>%
  transmute(
    Metric,
    Month,
    Compartment,
    contrast,
    direction,
    Tukey_estimate_CI,
    Tukey_p = tukey_p_round,
    tukey_sig,
    Wilcox_p = wilcox_p_round,
    wilcox_sig
  ) %>%
  arrange(Compartment, Metric)

write_csv(
  alpha_key_contrasts_compact,
  file.path(outdir_alpha, "07_alpha_key_contrasts_compact_for_paper_ITS.csv")
)

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
  transmute(
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
  arrange(Compartment, Metric)

write_csv(
  alpha_key_contrasts_writing,
  file.path(outdir_alpha, "08_alpha_key_contrasts_writing_table_ITS.csv")
)

cat("\n====================================================\n")
cat("ALPHA KEY CONTRASTS — COMPACT PAPER TABLE\n")
cat("====================================================\n")
print(alpha_key_contrasts_compact, n = Inf, width = Inf)

cat("\n====================================================\n")
cat("ALPHA KEY CONTRASTS — WRITING TABLE\n")
cat("====================================================\n")
print(alpha_key_contrasts_writing, n = Inf, width = Inf)





setwd("path/to/QIIME2_files_to_R_ITS")

library(readr)
library(dplyr)
library(stringr)

summary_ITS <- read_csv("01_alpha_diversity_ITS/09_alpha_table_for_writing_ITS.csv", show_col_types = FALSE) %>%
  mutate(
    Month = "October",
    Compartment = case_when(
      str_detect(Health_status, "Rhizosphere") ~ "Rhizosphere",
      str_detect(Health_status, "Endosphere") ~ "Endosphere",
      str_detect(Health_status, "Phyllosphere") ~ "Phyllosphere",
      Health_status == "Bulk soil" ~ "Bulk soil"
    )
  )

tests_ITS <- read_csv("01_alpha_diversity_ITS/08_alpha_key_contrasts_writing_table_ITS.csv", show_col_types = FALSE)

alphaITS_master <- summary_ITS %>%
  left_join(
    tests_ITS,
    by = c("Metric", "Month", "Compartment")
  )

write_csv(alphaITS_master, "01_alpha_diversity_ITS/Table_S3_alpha_diversity_ITS_master.csv")

print(alphaITS_master, n = Inf, width = Inf)



###############################################################################################################################################################
###############################################################################################################################################################
###############################################################################################################################################################
###############################################################################################################################################################
###############################################################################################################################################################
############################################################
# SECTION 3 — BETA DIVERSITY (ITS, OCTOBER ONLY)
#
# Objetivos:
# 1) Calcular beta diversity no filogenética para ITS
# 2) Hacer análisis global con Bray-Curtis y Jaccard
# 3) Hacer análisis por compartimento: Healthy vs Diseased
# 4) Exportar tablas y figuras base para Results
# 5) Dibujar elipses en los PCoA
############################################################



setwd("path/to/QIIME2_files_to_R_ITS")

suppressPackageStartupMessages({
  library(phyloseq)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(vegan)
  library(ape)
  library(ggplot2)
  library(readr)
  library(patchwork)
})

############################################################
# 0) INPUTS AND OUTPUT
############################################################

physeq_path <- "Phyloseq_ITS.rds"

outdir_beta <- "02_beta_diversity_ITS"
dir.create(outdir_beta, recursive = TRUE, showWarnings = FALSE)

hs_levels <- c(
  "Bulk soil",
  "Healthy Rhizosphere", "Diseased Rhizosphere",
  "Healthy Endosphere",  "Diseased Endosphere",
  "Healthy Phyllosphere","Diseased Phyllosphere"
)

hs_colors <- c(
  "Bulk soil"              = "#896054",
  "Healthy Rhizosphere"    = "#4B9EA7",
  "Diseased Rhizosphere"   = "#801723",
  "Healthy Endosphere"     = "#6B9863",
  "Diseased Endosphere"    = "#E7C313",
  "Healthy Phyllosphere"   = "#4A5263",
  "Diseased Phyllosphere"  = "#6B4EA7"
)

############################################################
# 1) LOAD PHYLOSEQ
############################################################

ps0 <- readRDS(physeq_path)
ps0 <- prune_samples(sample_sums(ps0) > 0, ps0)
ps0 <- prune_taxa(taxa_sums(ps0) > 0, ps0)

md0 <- as(sample_data(ps0), "data.frame") %>%
  rownames_to_column("SampleID") %>%
  mutate(
    Month = as.character(Month),
    Compartment = as.character(Compartment),
    Health_status = as.character(Health_status)
  )

cat("\n====================================================\n")
cat("PHYLOSEQ LOADED\n")
cat("====================================================\n")
print(ps0)

cat("\nMonth table:\n")
print(table(md0$Month, useNA = "ifany"))

cat("\nCompartment table:\n")
print(table(md0$Compartment, useNA = "ifany"))

cat("\nHealth_status table:\n")
print(table(md0$Health_status, useNA = "ifany"))

############################################################
# 2) HELPER FUNCTIONS
############################################################

get_otu_samples_rows <- function(ps) {
  mat <- as(otu_table(ps), "matrix")
  if (taxa_are_rows(ps)) mat <- t(mat)
  mat
}

make_bray <- function(ps) {
  mat <- get_otu_samples_rows(ps)
  vegan::vegdist(mat, method = "bray")
}

make_jaccard <- function(ps) {
  mat <- get_otu_samples_rows(ps)
  mat_pa <- ifelse(mat > 0, 1, 0)
  vegan::vegdist(mat_pa, method = "jaccard", binary = TRUE)
}

run_pcoa_df <- function(dist_obj, md_df, metric_name) {
  pcoa_res <- ape::pcoa(as.matrix(dist_obj))
  
  coords <- as.data.frame(pcoa_res$vectors[, 1:2, drop = FALSE]) %>%
    rownames_to_column("SampleID") %>%
    left_join(md_df, by = "SampleID") %>%
    mutate(Metric = metric_name)
  
  colnames(coords)[2:3] <- c("PCoA1", "PCoA2")
  
  var1 <- pcoa_res$values$Relative_eig[1] * 100
  var2 <- pcoa_res$values$Relative_eig[2] * 100
  
  list(coords = coords, var1 = var1, var2 = var2)
}

safe_anosim <- function(dist_obj, grouping) {
  set.seed(123)
  vegan::anosim(dist_obj, grouping = grouping, permutations = 999)
}

safe_betadisper_test <- function(dist_obj, grouping) {
  bd <- vegan::betadisper(dist_obj, grouping)
  an <- anova(bd)
  list(bd = bd, anova = an)
}

pairwise_distance_summary <- function(dist_obj, metadata_df, group_var = "Health_status") {
  mat <- as.matrix(dist_obj)
  samp <- rownames(mat)
  
  pairs <- t(combn(samp, 2))
  
  pairs_df <- tibble(
    Sample1 = pairs[, 1],
    Sample2 = pairs[, 2]
  ) %>%
    mutate(
      distance = purrr::map2_dbl(Sample1, Sample2, ~ mat[.x, .y]),
      group1 = metadata_df[[group_var]][match(Sample1, metadata_df$SampleID)],
      group2 = metadata_df[[group_var]][match(Sample2, metadata_df$SampleID)]
    ) %>%
    mutate(
      comparison_type = case_when(
        group1 == group2 ~ "Intra-group",
        TRUE ~ "Inter-group"
      )
    )
  
  pairs_df %>%
    group_by(comparison_type) %>%
    summarise(
      mean_distance = mean(distance, na.rm = TRUE),
      sd_distance = sd(distance, na.rm = TRUE),
      .groups = "drop"
    )
}

classify_beta_result <- function(permanova_p, betadisper_p) {
  case_when(
    is.na(permanova_p) ~ "Not tested",
    permanova_p < 0.05 & !is.na(betadisper_p) & betadisper_p < 0.05 ~ "Significant, but dispersion differs",
    permanova_p < 0.05 ~ "Significant centroid shift",
    TRUE ~ "No significant difference"
  )
}

plot_global_pcoa <- function(df, metric_name, var1, var2) {
  ggplot(df, aes(x = PCoA1, y = PCoA2, color = Health_status, fill = Health_status)) +
    stat_ellipse(
      aes(group = Health_status),
      type = "norm",
      geom = "polygon",
      alpha = 0.12,
      color = NA,
      level = 0.68,
      show.legend = FALSE
    ) +
    stat_ellipse(
      aes(group = Health_status),
      type = "norm",
      linewidth = 0.8,
      level = 0.68,
      show.legend = FALSE
    ) +
    geom_point(size = 3, alpha = 0.9) +
    scale_color_manual(values = hs_colors, breaks = hs_levels, drop = FALSE) +
    scale_fill_manual(values = hs_colors, breaks = hs_levels, drop = FALSE) +
    labs(
      title = paste("ITS beta diversity —", metric_name),
      subtitle = "October only",
      x = paste0("PCoA1 (", round(var1, 1), "%)"),
      y = paste0("PCoA2 (", round(var2, 1), "%)")
    ) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      legend.title = element_blank(),
      plot.title = element_text(face = "bold")
    )
}

plot_compartment_pcoa <- function(df, metric_name, comp_name, var1, var2) {
  ggplot(df, aes(x = PCoA1, y = PCoA2, color = Health_status, fill = Health_status)) +
    stat_ellipse(
      aes(group = Health_status),
      type = "norm",
      geom = "polygon",
      alpha = 0.15,
      color = NA,
      level = 0.68,
      show.legend = FALSE
    ) +
    stat_ellipse(
      aes(group = Health_status),
      type = "norm",
      linewidth = 0.9,
      level = 0.68,
      show.legend = FALSE
    ) +
    geom_point(size = 3, alpha = 0.9) +
    scale_color_manual(values = hs_colors, breaks = hs_levels, drop = FALSE) +
    scale_fill_manual(values = hs_colors, breaks = hs_levels, drop = FALSE) +
    labs(
      title = paste(comp_name, "—", metric_name),
      x = paste0("PCoA1 (", round(var1, 1), "%)"),
      y = paste0("PCoA2 (", round(var2, 1), "%)")
    ) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      legend.title = element_blank(),
      plot.title = element_text(face = "bold")
    )
}

############################################################
# 3) GLOBAL ANALYSIS
############################################################

run_global_beta <- function(ps_in) {
  md <- as(sample_data(ps_in), "data.frame") %>%
    rownames_to_column("SampleID")
  
  dist_bray <- make_bray(ps_in)
  dist_jacc <- make_jaccard(ps_in)
  
  pcoa_bray <- run_pcoa_df(dist_bray, md, "Bray-Curtis")
  pcoa_jacc <- run_pcoa_df(dist_jacc, md, "Jaccard")
  
  set.seed(123)
  ad_bray <- vegan::adonis2(dist_bray ~ Health_status, data = md, permutations = 999)
  
  set.seed(123)
  ad_jacc <- vegan::adonis2(dist_jacc ~ Health_status, data = md, permutations = 999)
  
  an_bray <- safe_anosim(dist_bray, md$Health_status)
  an_jacc <- safe_anosim(dist_jacc, md$Health_status)
  
  bd_bray <- safe_betadisper_test(dist_bray, md$Health_status)
  bd_jacc <- safe_betadisper_test(dist_jacc, md$Health_status)
  
  global_summary <- tibble(
    Month = "October",
    Metric = c("Bray-Curtis", "Jaccard"),
    n_samples = nsamples(ps_in),
    n_taxa = ntaxa(ps_in),
    PCoA_axis1_var = c(pcoa_bray$var1, pcoa_jacc$var1),
    PCoA_axis2_var = c(pcoa_bray$var2, pcoa_jacc$var2),
    PERMANOVA_R2 = c(ad_bray$R2[1], ad_jacc$R2[1]),
    PERMANOVA_F = c(ad_bray$F[1], ad_jacc$F[1]),
    PERMANOVA_p = c(ad_bray$`Pr(>F)`[1], ad_jacc$`Pr(>F)`[1]),
    ANOSIM_R = c(as.numeric(an_bray$statistic), as.numeric(an_jacc$statistic)),
    ANOSIM_p = c(an_bray$signif, an_jacc$signif),
    betadisper_F = c(bd_bray$anova$`F value`[1], bd_jacc$anova$`F value`[1]),
    betadisper_p = c(bd_bray$anova$`Pr(>F)`[1], bd_jacc$anova$`Pr(>F)`[1])
  )
  
  list(
    summary = global_summary,
    pcoa_bray = pcoa_bray$coords,
    pcoa_jacc = pcoa_jacc$coords
  )
}

global_beta <- run_global_beta(ps0)
global_beta_summary <- global_beta$summary

write_csv(
  global_beta_summary,
  file.path(outdir_beta, "01_global_beta_summary_ITS.csv")
)

write_csv(
  global_beta$pcoa_bray,
  file.path(outdir_beta, "02_global_pcoa_bray_ITS.csv")
)

write_csv(
  global_beta$pcoa_jacc,
  file.path(outdir_beta, "03_global_pcoa_jaccard_ITS.csv")
)

############################################################
# 4) GLOBAL PLOTS
############################################################

p_global_bray <- plot_global_pcoa(
  global_beta$pcoa_bray,
  "Bray-Curtis",
  global_beta_summary$PCoA_axis1_var[global_beta_summary$Metric == "Bray-Curtis"],
  global_beta_summary$PCoA_axis2_var[global_beta_summary$Metric == "Bray-Curtis"]
)

p_global_jacc <- plot_global_pcoa(
  global_beta$pcoa_jacc,
  "Jaccard",
  global_beta_summary$PCoA_axis1_var[global_beta_summary$Metric == "Jaccard"],
  global_beta_summary$PCoA_axis2_var[global_beta_summary$Metric == "Jaccard"]
)

p_global <- p_global_bray / p_global_jacc
print(p_global)

ggsave(
  file.path(outdir_beta, "04_global_beta_pcoa_ITS.png"),
  plot = p_global,
  width = 10, height = 12, dpi = 600
)

############################################################
# 5) WITHIN-COMPARTMENT ANALYSIS
############################################################

prepare_subset_for_compartment <- function(ps_in, compartment_name, healthy_label, diseased_label) {
  md <- as(sample_data(ps_in), "data.frame")
  
  keep_ids <- rownames(md)[
    as.character(md$Compartment) == compartment_name &
      as.character(md$Health_status) %in% c(healthy_label, diseased_label)
  ]
  
  ps_sub <- prune_samples(keep_ids, ps_in)
  ps_sub <- prune_samples(sample_sums(ps_sub) > 0, ps_sub)
  ps_sub <- prune_taxa(taxa_sums(ps_sub) > 0, ps_sub)
  ps_sub
}

run_beta_one_compartment <- function(ps_in, compartment_name, healthy_label, diseased_label) {
  ps_sub <- prepare_subset_for_compartment(ps_in, compartment_name, healthy_label, diseased_label)
  
  md <- as(sample_data(ps_sub), "data.frame") %>%
    rownames_to_column("SampleID")
  
  dist_bray <- make_bray(ps_sub)
  dist_jacc <- make_jaccard(ps_sub)
  
  pcoa_bray <- run_pcoa_df(dist_bray, md, "Bray-Curtis")
  pcoa_jacc <- run_pcoa_df(dist_jacc, md, "Jaccard")
  
  set.seed(123)
  ad_bray <- vegan::adonis2(dist_bray ~ Health_status, data = md, permutations = 999)
  
  set.seed(123)
  ad_jacc <- vegan::adonis2(dist_jacc ~ Health_status, data = md, permutations = 999)
  
  an_bray <- safe_anosim(dist_bray, md$Health_status)
  an_jacc <- safe_anosim(dist_jacc, md$Health_status)
  
  bd_bray <- safe_betadisper_test(dist_bray, md$Health_status)
  bd_jacc <- safe_betadisper_test(dist_jacc, md$Health_status)
  
  pair_bray <- pairwise_distance_summary(dist_bray, md)
  pair_jacc <- pairwise_distance_summary(dist_jacc, md)
  
  summary_tbl <- tibble(
    Month = "October",
    Compartment = compartment_name,
    n_samples = nsamples(ps_sub),
    n_Healthy = sum(md$Health_status == healthy_label),
    n_Diseased = sum(md$Health_status == diseased_label),
    
    Bray_PCoA_axis1_var = pcoa_bray$var1,
    Bray_PCoA_axis2_var = pcoa_bray$var2,
    Bray_PERMANOVA_R2 = ad_bray$R2[1],
    Bray_PERMANOVA_F = ad_bray$F[1],
    Bray_PERMANOVA_p = ad_bray$`Pr(>F)`[1],
    Bray_ANOSIM_R = as.numeric(an_bray$statistic),
    Bray_ANOSIM_p = an_bray$signif,
    Bray_betadisper_F = bd_bray$anova$`F value`[1],
    Bray_betadisper_p = bd_bray$anova$`Pr(>F)`[1],
    Bray_mean_Intra = pair_bray$mean_distance[pair_bray$comparison_type == "Intra-group"],
    Bray_sd_Intra = pair_bray$sd_distance[pair_bray$comparison_type == "Intra-group"],
    Bray_mean_Inter = pair_bray$mean_distance[pair_bray$comparison_type == "Inter-group"],
    Bray_sd_Inter = pair_bray$sd_distance[pair_bray$comparison_type == "Inter-group"],
    
    Jaccard_PCoA_axis1_var = pcoa_jacc$var1,
    Jaccard_PCoA_axis2_var = pcoa_jacc$var2,
    Jaccard_PERMANOVA_R2 = ad_jacc$R2[1],
    Jaccard_PERMANOVA_F = ad_jacc$F[1],
    Jaccard_PERMANOVA_p = ad_jacc$`Pr(>F)`[1],
    Jaccard_ANOSIM_R = as.numeric(an_jacc$statistic),
    Jaccard_ANOSIM_p = an_jacc$signif,
    Jaccard_betadisper_F = bd_jacc$anova$`F value`[1],
    Jaccard_betadisper_p = bd_jacc$anova$`Pr(>F)`[1],
    Jaccard_mean_Intra = pair_jacc$mean_distance[pair_jacc$comparison_type == "Intra-group"],
    Jaccard_sd_Intra = pair_jacc$sd_distance[pair_jacc$comparison_type == "Intra-group"],
    Jaccard_mean_Inter = pair_jacc$mean_distance[pair_jacc$comparison_type == "Inter-group"],
    Jaccard_sd_Inter = pair_jacc$sd_distance[pair_jacc$comparison_type == "Inter-group"]
  )
  
  list(
    summary = summary_tbl,
    pcoa_bray = pcoa_bray$coords,
    pcoa_jacc = pcoa_jacc$coords
  )
}

comp_list <- list(
  list(comp = "Rhizosphere", healthy = "Healthy Rhizosphere", diseased = "Diseased Rhizosphere"),
  list(comp = "Endosphere", healthy = "Healthy Endosphere", diseased = "Diseased Endosphere"),
  list(comp = "Phyllosphere", healthy = "Healthy Phyllosphere", diseased = "Diseased Phyllosphere")
)

beta_results <- lapply(comp_list, function(x) {
  run_beta_one_compartment(ps0, x$comp, x$healthy, x$diseased)
})

beta_by_compartment_summary <- bind_rows(lapply(beta_results, `[[`, "summary"))

write_csv(
  beta_by_compartment_summary,
  file.path(outdir_beta, "05_beta_by_compartment_summary_ITS.csv")
)

############################################################
# 6) KEY TABLE FOR PAPER
############################################################

beta_key_table <- beta_by_compartment_summary %>%
  transmute(
    Month,
    Compartment,
    
    Bray_result = classify_beta_result(Bray_PERMANOVA_p, Bray_betadisper_p),
    Bray_R2 = Bray_PERMANOVA_R2,
    Bray_F = Bray_PERMANOVA_F,
    Bray_p = Bray_PERMANOVA_p,
    Bray_ANOSIM_R = Bray_ANOSIM_R,
    Bray_betadisper_p = Bray_betadisper_p,
    
    Jaccard_result = classify_beta_result(Jaccard_PERMANOVA_p, Jaccard_betadisper_p),
    Jaccard_R2 = Jaccard_PERMANOVA_R2,
    Jaccard_F = Jaccard_PERMANOVA_F,
    Jaccard_p = Jaccard_PERMANOVA_p,
    Jaccard_ANOSIM_R = Jaccard_ANOSIM_R,
    Jaccard_betadisper_p = Jaccard_betadisper_p
  ) %>%
  arrange(Compartment)

write_csv(
  beta_key_table,
  file.path(outdir_beta, "06_beta_key_table_for_paper_ITS.csv")
)

############################################################
# 7) COMPARTMENT PLOTS
############################################################

for (i in seq_along(beta_results)) {
  res_i <- beta_results[[i]]
  comp_i <- res_i$summary$Compartment[1]
  
  p1 <- plot_compartment_pcoa(
    res_i$pcoa_bray,
    "Bray-Curtis",
    comp_i,
    res_i$summary$Bray_PCoA_axis1_var[1],
    res_i$summary$Bray_PCoA_axis2_var[1]
  )
  
  p2 <- plot_compartment_pcoa(
    res_i$pcoa_jacc,
    "Jaccard",
    comp_i,
    res_i$summary$Jaccard_PCoA_axis1_var[1],
    res_i$summary$Jaccard_PCoA_axis2_var[1]
  )
  
  p_comp <- p1 / p2
  print(p_comp)
  
  ggsave(
    file.path(outdir_beta, paste0("07_pcoa_", comp_i, "_ITS.png")),
    plot = p_comp,
    width = 9, height = 11, dpi = 600
  )
}

############################################################
# 8) CONSOLE OUTPUT
############################################################

cat("\n====================================================\n")
cat("GLOBAL BETA SUMMARY\n")
cat("====================================================\n")
print(global_beta_summary, n = Inf, width = Inf)

cat("\n====================================================\n")
cat("BETA SUMMARY BY COMPARTMENT\n")
cat("====================================================\n")
print(beta_by_compartment_summary, n = Inf, width = Inf)

cat("\n====================================================\n")
cat("BETA KEY TABLE FOR PAPER\n")
cat("====================================================\n")
print(beta_key_table, n = Inf, width = Inf)

cat("\n====================================================\n")
cat("FILES WRITTEN TO:\n")
cat("====================================================\n")
print(normalizePath(outdir_beta))




setwd("path/to/QIIME2_files_to_R_ITS")

library(readr)
library(dplyr)

betaITS_global <- read_csv("02_beta_diversity_ITS/01_global_beta_summary_ITS.csv", show_col_types = FALSE)
betaITS_comp   <- read_csv("02_beta_diversity_ITS/05_beta_by_compartment_summary_ITS.csv", show_col_types = FALSE)
betaITS_key    <- read_csv("02_beta_diversity_ITS/06_beta_key_table_for_paper_ITS.csv", show_col_types = FALSE)

write_csv(betaITS_global, "02_beta_diversity_ITS/Table_S7_beta_ITS_global_master.csv")
write_csv(betaITS_comp, "02_beta_diversity_ITS/Table_S8_beta_ITS_compartment_master.csv")

print(betaITS_global, n = Inf, width = Inf)
print(betaITS_key, n = Inf, width = Inf)



###############################################################################################################################################################
###############################################################################################################################################################
###############################################################################################################################################################
###############################################################################################################################################################
###############################################################################################################################################################


############################################################
# SECTION 4 — TAXONOMY (ITS, OCTOBER ONLY)
#
# Objetivos:
# 1) Resumir taxonomía ITS a niveles altos
# 2) Generar tablas de abundancia relativa por grupo
# 3) Extraer taxa diana (Rosellinia, R. necatrix,
#    Botryosphaeriaceae y Neofusicoccum)
# 4) Graficar abundancia relativa de esos taxa diana
############################################################



setwd("path/to/QIIME2_files_to_R_ITS")

suppressPackageStartupMessages({
  library(phyloseq)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(readr)
  library(stringr)
  library(forcats)
})

############################################################
# 0) INPUTS AND OUTPUT
############################################################

physeq_path <- "Phyloseq_ITS.rds"

outdir_tax <- "03_taxonomy_ITS"
dir.create(outdir_tax, recursive = TRUE, showWarnings = FALSE)

outdir_tax_tables <- file.path(outdir_tax, "tables")
dir.create(outdir_tax_tables, recursive = TRUE, showWarnings = FALSE)

outdir_tax_figs <- file.path(outdir_tax, "figures")
dir.create(outdir_tax_figs, recursive = TRUE, showWarnings = FALSE)

hs_levels <- c(
  "Bulk soil",
  "Healthy Rhizosphere", "Diseased Rhizosphere",
  "Healthy Endosphere",  "Diseased Endosphere",
  "Healthy Phyllosphere","Diseased Phyllosphere"
)

hs_colors <- c(
  "Bulk soil"              = "#896054",
  "Healthy Rhizosphere"    = "#4B9EA7",
  "Diseased Rhizosphere"   = "#801723",
  "Healthy Endosphere"     = "#6B9863",
  "Diseased Endosphere"    = "#E7C313",
  "Healthy Phyllosphere"   = "#4A5263",
  "Diseased Phyllosphere"  = "#6B4EA7"
)

############################################################
# 1) LOAD PHYLOSEQ
############################################################

ps0 <- readRDS(physeq_path)
ps0 <- prune_samples(sample_sums(ps0) > 0, ps0)
ps0 <- prune_taxa(taxa_sums(ps0) > 0, ps0)

md0 <- as(sample_data(ps0), "data.frame") %>%
  rownames_to_column("SampleID") %>%
  mutate(
    Month = as.character(Month),
    Compartment = as.character(Compartment),
    Health_status = factor(as.character(Health_status), levels = hs_levels)
  )

otu_mat <- as(otu_table(ps0), "matrix")
if (taxa_are_rows(ps0)) {
  otu_mat <- otu_mat
} else {
  otu_mat <- t(otu_mat)
}

tax_df <- as.data.frame(tax_table(ps0)) %>%
  rownames_to_column("ASV_ID") %>%
  as_tibble() %>%
  mutate(
    Kingdom = as.character(Kingdom),
    Phylum = as.character(Phylum),
    Class = as.character(Class),
    Order = as.character(Order),
    Family = as.character(Family),
    Genus = as.character(Genus),
    Species = as.character(Species)
  )

cat("\n====================================================\n")
cat("PHYLOSEQ LOADED\n")
cat("====================================================\n")
print(ps0)

############################################################
# 2) HELPER FUNCTIONS
############################################################

clean_tax_label <- function(x, fallback) {
  x <- as.character(x)
  x[is.na(x) | x == ""] <- fallback
  x
}

make_tax_rel_long <- function(rank_name) {
  rank_sym <- rlang::sym(rank_name)
  
  rel_mat <- apply(otu_mat, 2, function(x) {
    s <- sum(x)
    if (s == 0) rep(0, length(x)) else x / s
  })
  
  if (is.vector(rel_mat)) {
    rel_mat <- matrix(rel_mat, ncol = 1)
    colnames(rel_mat) <- colnames(otu_mat)
    rownames(rel_mat) <- rownames(otu_mat)
  }
  
  rel_df <- as.data.frame(rel_mat) %>%
    rownames_to_column("ASV_ID") %>%
    pivot_longer(
      cols = -ASV_ID,
      names_to = "SampleID",
      values_to = "Relative_abundance"
    ) %>%
    left_join(tax_df, by = "ASV_ID") %>%
    left_join(md0, by = "SampleID") %>%
    mutate(
      Taxon = clean_tax_label(!!rank_sym, paste0("Unclassified_", rank_name))
    ) %>%
    group_by(SampleID, Month, Compartment, Health_status, Taxon) %>%
    summarise(
      Relative_abundance = sum(Relative_abundance, na.rm = TRUE),
      .groups = "drop"
    )
  
  rel_df
}

summarise_tax_group <- function(rel_long) {
  rel_long %>%
    group_by(Month, Compartment, Health_status, Taxon) %>%
    summarise(
      mean_rel_abundance = mean(Relative_abundance, na.rm = TRUE),
      sd_rel_abundance = sd(Relative_abundance, na.rm = TRUE),
      se_rel_abundance = sd_rel_abundance / sqrt(dplyr::n()),
      prevalence = sum(Relative_abundance > 0, na.rm = TRUE),
      prevalence_pct = 100 * prevalence / dplyr::n(),
      mean_rel_abundance_pct = 100 * mean_rel_abundance,
      .groups = "drop"
    )
}

save_plot_both <- function(plot_obj, base_filename, width = 8, height = 5) {
  ggsave(
    file.path(outdir_tax_figs, paste0(base_filename, ".png")),
    plot = plot_obj, width = width, height = height, dpi = 600
  )
  ggsave(
    file.path(outdir_tax_figs, paste0(base_filename, ".tiff")),
    plot = plot_obj, width = width, height = height, dpi = 600, compression = "lzw"
  )
}

plot_target_abundance <- function(df, title_text, y_text = "Relative abundance (%)") {
  ggplot(df, aes(x = Health_status, y = mean_rel_abundance_pct, fill = Health_status)) +
    geom_col(width = 0.75, color = "black", linewidth = 0.3) +
    geom_errorbar(
      aes(ymin = pmax(mean_rel_abundance_pct - se_rel_abundance_pct, 0),
          ymax = mean_rel_abundance_pct + se_rel_abundance_pct),
      width = 0.15,
      linewidth = 0.5
    ) +
    geom_jitter(
      data = df %>% select(Health_status, SampleID, Relative_abundance_pct),
      aes(x = Health_status, y = Relative_abundance_pct),
      width = 0.12,
      height = 0,
      size = 2,
      shape = 21,
      color = "black",
      inherit.aes = FALSE
    ) +
    scale_fill_manual(values = hs_colors, drop = FALSE) +
    labs(
      title = title_text,
      x = NULL,
      y = y_text
    ) +
    theme_bw() +
    theme(
      legend.position = "none",
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 20, hjust = 1),
      plot.title = element_text(face = "bold")
    )
}

############################################################
# 3) GLOBAL TAXONOMIC TABLES
############################################################

phylum_long  <- make_tax_rel_long("Phylum")
family_long  <- make_tax_rel_long("Family")
genus_long   <- make_tax_rel_long("Genus")
species_long <- make_tax_rel_long("Species")

phylum_summary  <- summarise_tax_group(phylum_long)
family_summary  <- summarise_tax_group(family_long)
genus_summary   <- summarise_tax_group(genus_long)
species_summary <- summarise_tax_group(species_long)

write_csv(phylum_long,  file.path(outdir_tax_tables, "01_phylum_relative_by_sample_ITS.csv"))
write_csv(family_long,  file.path(outdir_tax_tables, "02_family_relative_by_sample_ITS.csv"))
write_csv(genus_long,   file.path(outdir_tax_tables, "03_genus_relative_by_sample_ITS.csv"))
write_csv(species_long, file.path(outdir_tax_tables, "04_species_relative_by_sample_ITS.csv"))

write_csv(phylum_summary,  file.path(outdir_tax_tables, "05_phylum_group_summary_ITS.csv"))
write_csv(family_summary,  file.path(outdir_tax_tables, "06_family_group_summary_ITS.csv"))
write_csv(genus_summary,   file.path(outdir_tax_tables, "07_genus_group_summary_ITS.csv"))
write_csv(species_summary, file.path(outdir_tax_tables, "08_species_group_summary_ITS.csv"))

############################################################
# 4) TOP TAXA TABLES
############################################################

top15_phyla <- phylum_summary %>%
  group_by(Taxon) %>%
  summarise(global_mean = mean(mean_rel_abundance, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(global_mean)) %>%
  slice_head(n = 15) %>%
  pull(Taxon)

top15_phyla

top20_families <- family_summary %>%
  group_by(Taxon) %>%
  summarise(global_mean = mean(mean_rel_abundance, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(global_mean)) %>%
  slice_head(n = 20) %>%
  pull(Taxon)

top20_genera <- genus_summary %>%
  group_by(Taxon) %>%
  summarise(global_mean = mean(mean_rel_abundance, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(global_mean)) %>%
  slice_head(n = 20) %>%
  pull(Taxon)
top20_genera

write_csv(
  tibble(top15_phyla = top15_phyla),
  file.path(outdir_tax_tables, "09_top15_phyla_ITS.csv")
)

write_csv(
  tibble(top20_families = top20_families),
  file.path(outdir_tax_tables, "10_top20_families_ITS.csv")
)

write_csv(
  tibble(top20_genera = top20_genera),
  file.path(outdir_tax_tables, "11_top20_genera_ITS.csv")
)

############################################################
# 5) TARGET TAXA — BUILD SAMPLE-LEVEL TABLES
############################################################

target_taxa_sample <- tax_df %>%
  transmute(
    ASV_ID,
    Family_clean  = clean_tax_label(Family, "Unclassified_Family"),
    Genus_clean   = clean_tax_label(Genus, "Unclassified_Genus"),
    Species_clean = clean_tax_label(Species, "Unclassified_Species")
  ) %>%
  left_join(
    as.data.frame(otu_mat) %>%
      rownames_to_column("ASV_ID") %>%
      pivot_longer(cols = -ASV_ID, names_to = "SampleID", values_to = "Reads"),
    by = "ASV_ID"
  ) %>%
  left_join(md0, by = "SampleID") %>%
  group_by(SampleID, Month, Compartment, Health_status) %>%
  mutate(
    Total_reads_sample = sum(Reads, na.rm = TRUE)
  ) %>%
  ungroup()

# Rosellinia genus — rhizosphere + endosphere
rosellinia_genus_sample <- target_taxa_sample %>%
  filter(
    Compartment %in% c("Rhizosphere", "Endosphere"),
    Health_status %in% c(
      "Healthy Rhizosphere", "Diseased Rhizosphere",
      "Healthy Endosphere", "Diseased Endosphere"
    )
  ) %>%
  mutate(Target_reads = ifelse(Genus_clean == "Rosellinia", Reads, 0)) %>%
  group_by(SampleID, Month, Compartment, Health_status, Total_reads_sample) %>%
  summarise(Target_reads = sum(Target_reads, na.rm = TRUE), .groups = "drop") %>%
  mutate(Relative_abundance_pct = 100 * Target_reads / Total_reads_sample)

# Rosellinia necatrix species — rhizosphere + endosphere
rosellinia_necatrix_sample <- target_taxa_sample %>%
  filter(
    Compartment %in% c("Rhizosphere", "Endosphere"),
    Health_status %in% c(
      "Healthy Rhizosphere", "Diseased Rhizosphere",
      "Healthy Endosphere", "Diseased Endosphere"
    )
  ) %>%
  mutate(
    Target_reads = ifelse(
      str_to_lower(Species_clean) == "rosellinia_necatrix",
      Reads, 0
    )
  ) %>%
  group_by(SampleID, Month, Compartment, Health_status, Total_reads_sample) %>%
  summarise(Target_reads = sum(Target_reads, na.rm = TRUE), .groups = "drop") %>%
  mutate(Relative_abundance_pct = 100 * Target_reads / Total_reads_sample)

# Botryosphaeriaceae family — phyllosphere
botryo_family_sample <- target_taxa_sample %>%
  filter(
    Compartment == "Phyllosphere",
    Health_status %in% c("Healthy Phyllosphere", "Diseased Phyllosphere")
  ) %>%
  mutate(Target_reads = ifelse(Family_clean == "Botryosphaeriaceae", Reads, 0)) %>%
  group_by(SampleID, Month, Compartment, Health_status, Total_reads_sample) %>%
  summarise(Target_reads = sum(Target_reads, na.rm = TRUE), .groups = "drop") %>%
  mutate(Relative_abundance_pct = 100 * Target_reads / Total_reads_sample)

# Neofusicoccum genus — phyllosphere
neofusicoccum_sample <- target_taxa_sample %>%
  filter(
    Compartment == "Phyllosphere",
    Health_status %in% c("Healthy Phyllosphere", "Diseased Phyllosphere")
  ) %>%
  mutate(Target_reads = ifelse(Genus_clean == "Neofusicoccum", Reads, 0)) %>%
  group_by(SampleID, Month, Compartment, Health_status, Total_reads_sample) %>%
  summarise(Target_reads = sum(Target_reads, na.rm = TRUE), .groups = "drop") %>%
  mutate(Relative_abundance_pct = 100 * Target_reads / Total_reads_sample)

write_csv(rosellinia_genus_sample,      file.path(outdir_tax_tables, "12_rosellinia_genus_sample_level_ITS.csv"))
write_csv(rosellinia_necatrix_sample,   file.path(outdir_tax_tables, "13_rosellinia_necatrix_sample_level_ITS.csv"))
write_csv(botryo_family_sample,         file.path(outdir_tax_tables, "14_botryosphaeriaceae_sample_level_ITS.csv"))
write_csv(neofusicoccum_sample,         file.path(outdir_tax_tables, "15_neofusicoccum_sample_level_ITS.csv"))

############################################################
# 6) TARGET TAXA — GROUP SUMMARIES
############################################################

summarise_target <- function(df) {
  df %>%
    group_by(Month, Compartment, Health_status) %>%
    summarise(
      n = dplyr::n(),
      mean_rel_abundance_pct = mean(Relative_abundance_pct, na.rm = TRUE),
      sd_rel_abundance_pct = sd(Relative_abundance_pct, na.rm = TRUE),
      se_rel_abundance_pct = sd_rel_abundance_pct / sqrt(n),
      prevalence = sum(Relative_abundance_pct > 0, na.rm = TRUE),
      prevalence_pct = 100 * prevalence / n,
      .groups = "drop"
    )
}

rosellinia_genus_summary    <- summarise_target(rosellinia_genus_sample)
rosellinia_necatrix_summary <- summarise_target(rosellinia_necatrix_sample)
botryo_family_summary       <- summarise_target(botryo_family_sample)
neofusicoccum_summary       <- summarise_target(neofusicoccum_sample)

write_csv(rosellinia_genus_summary,    file.path(outdir_tax_tables, "16_rosellinia_genus_group_summary_ITS.csv"))
write_csv(rosellinia_necatrix_summary, file.path(outdir_tax_tables, "17_rosellinia_necatrix_group_summary_ITS.csv"))
write_csv(botryo_family_summary,       file.path(outdir_tax_tables, "18_botryosphaeriaceae_group_summary_ITS.csv"))
write_csv(neofusicoccum_summary,       file.path(outdir_tax_tables, "19_neofusicoccum_group_summary_ITS.csv"))

############################################################
# 7) TARGET TAXA — PLOTS
############################################################

rosellinia_genus_plot_df <- rosellinia_genus_summary %>%
  left_join(
    rosellinia_genus_sample %>%
      select(SampleID, Health_status, Relative_abundance_pct),
    by = "Health_status"
  )

rosellinia_necatrix_plot_df <- rosellinia_necatrix_summary %>%
  left_join(
    rosellinia_necatrix_sample %>%
      select(SampleID, Health_status, Relative_abundance_pct),
    by = "Health_status"
  )

botryo_family_plot_df <- botryo_family_summary %>%
  left_join(
    botryo_family_sample %>%
      select(SampleID, Health_status, Relative_abundance_pct),
    by = "Health_status"
  )

neofusicoccum_plot_df <- neofusicoccum_summary %>%
  left_join(
    neofusicoccum_sample %>%
      select(SampleID, Health_status, Relative_abundance_pct),
    by = "Health_status"
  )

p_ros_genus <- plot_target_abundance(
  rosellinia_genus_plot_df,
  "Rosellinia relative abundance (%) — Rhizosphere and Endosphere"
)

p_ros_nec <- plot_target_abundance(
  rosellinia_necatrix_plot_df,
  "Rosellinia necatrix relative abundance (%) — Rhizosphere and Endosphere"
)

p_botryo <- plot_target_abundance(
  botryo_family_plot_df,
  "Botryosphaeriaceae relative abundance (%) — Phyllosphere"
)

p_neofu <- plot_target_abundance(
  neofusicoccum_plot_df,
  "Neofusicoccum relative abundance (%) — Phyllosphere"
)

print(p_ros_genus)
print(p_ros_nec)
print(p_botryo)
print(p_neofu)

save_plot_both(p_ros_genus, "20_plot_rosellinia_genus_ITS", width = 9, height = 5.5)
save_plot_both(p_ros_nec,   "21_plot_rosellinia_necatrix_ITS", width = 9, height = 5.5)
save_plot_both(p_botryo,    "22_plot_botryosphaeriaceae_ITS", width = 7, height = 5.5)
save_plot_both(p_neofu,     "23_plot_neofusicoccum_ITS", width = 7, height = 5.5)

############################################################
# 8) CONSOLE OUTPUT
############################################################

cat("\n====================================================\n")
cat("TOP 15 PHYLA\n")
cat("====================================================\n")
print(top15_phyla)

cat("\n====================================================\n")
cat("TOP 20 FAMILIES\n")
cat("====================================================\n")
print(top20_families)

cat("\n====================================================\n")
cat("TOP 20 GENERA\n")
cat("====================================================\n")
print(top20_genera)

cat("\n====================================================\n")
cat("ROSELLINIA GENUS — GROUP SUMMARY\n")
cat("====================================================\n")
print(rosellinia_genus_summary, n = Inf, width = Inf)

cat("\n====================================================\n")
cat("ROSELLINIA NECATRIX — GROUP SUMMARY\n")
cat("====================================================\n")
print(rosellinia_necatrix_summary, n = Inf, width = Inf)

cat("\n====================================================\n")
cat("BOTRYOSPHAERIACEAE — GROUP SUMMARY\n")
cat("====================================================\n")
print(botryo_family_summary, n = Inf, width = Inf)

cat("\n====================================================\n")
cat("NEOFUSICOCCUM — GROUP SUMMARY\n")
cat("====================================================\n")
print(neofusicoccum_summary, n = Inf, width = Inf)

cat("\n====================================================\n")
cat("FILES WRITTEN TO:\n")
cat("====================================================\n")
print(normalizePath(outdir_tax))





setwd("path/to/QIIME2_files_to_R_ITS")

library(readr)
library(dplyr)

phylumITS <- read_csv("03_taxonomy_ITS/tables/05_phylum_group_summary_ITS.csv", show_col_types = FALSE)
familyITS <- read_csv("03_taxonomy_ITS/tables/06_family_group_summary_ITS.csv", show_col_types = FALSE)
genusITS  <- read_csv("03_taxonomy_ITS/tables/07_genus_group_summary_ITS.csv", show_col_types = FALSE)

cat("\n=== ITS TOP PHYLA BY GROUP ===\n")
phylumITS %>%
  group_by(Compartment, Health_status) %>%
  slice_max(mean_rel_abundance_pct, n = 5, with_ties = FALSE) %>%
  arrange(Compartment, Health_status, desc(mean_rel_abundance_pct)) %>%
  print(n = Inf, width = Inf)

cat("\n=== ITS TOP FAMILIES BY GROUP ===\n")
familyITS %>%
  group_by(Compartment, Health_status) %>%
  slice_max(mean_rel_abundance_pct, n = 8, with_ties = FALSE) %>%
  arrange(Compartment, Health_status, desc(mean_rel_abundance_pct)) %>%
  print(n = Inf, width = Inf)

cat("\n=== ITS TOP GENERA BY GROUP ===\n")
genusITS %>%
  group_by(Compartment, Health_status) %>%
  slice_max(mean_rel_abundance_pct, n = 8, with_ties = FALSE) %>%
  arrange(Compartment, Health_status, desc(mean_rel_abundance_pct)) %>%
  print(n = Inf, width = Inf)











###############################################################################################################################################################
###############################################################################################################################################################
###############################################################################################################################################################
###############################################################################################################################################################
###############################################################################################################################################################

############################################################
# SECTION 5 — DIFFERENTIAL ABUNDANCE (ITS, OCTOBER ONLY)
# DESeq2 at GENUS level only
#
# Objetivos:
# 1) Construir una matriz de counts por género a partir de ASVs
# 2) Hacer DESeq2 solo dentro de cada compartimento
# 3) Comparar Healthy vs Diseased en:
#    - Rhizosphere
#    - Endosphere
#    - Phyllosphere
# 4) Exportar tablas completas y significativas
############################################################



setwd("path/to/QIIME2_files_to_R_ITS")

suppressPackageStartupMessages({
  library(phyloseq)
  library(DESeq2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(stringr)
})

############################################################
# 0) INPUTS AND OUTPUT
############################################################

physeq_path <- "Phyloseq_ITS.rds"

outdir_deseq <- "04_differential_abundance_ITS_DESeq2_genus"
dir.create(outdir_deseq, recursive = TRUE, showWarnings = FALSE)

min_total_counts_genus <- 10
alpha_sig <- 0.05

############################################################
# 1) LOAD PHYLOSEQ
############################################################

ps0 <- readRDS(physeq_path)
ps0 <- prune_samples(sample_sums(ps0) > 0, ps0)
ps0 <- prune_taxa(taxa_sums(ps0) > 0, ps0)

md0 <- as(sample_data(ps0), "data.frame") %>%
  rownames_to_column("SampleID") %>%
  mutate(
    Month = as.character(Month),
    Compartment = as.character(Compartment),
    Health_status = as.character(Health_status)
  )

cat("\n====================================================\n")
cat("PHYLOSEQ LOADED\n")
cat("====================================================\n")
print(ps0)

############################################################
# 2) PREPARE COUNT MATRIX + TAXONOMY
############################################################

count_mat_asv <- as(otu_table(ps0), "matrix")
if (!taxa_are_rows(ps0)) {
  count_mat_asv <- t(count_mat_asv)
}

tax_df <- as.data.frame(tax_table(ps0)) %>%
  rownames_to_column("ASV_ID") %>%
  as_tibble() %>%
  mutate(
    Kingdom = as.character(Kingdom),
    Phylum  = as.character(Phylum),
    Class   = as.character(Class),
    Order   = as.character(Order),
    Family  = as.character(Family),
    Genus   = as.character(Genus),
    Species = as.character(Species)
  )

cat("\n====================================================\n")
cat("ASV COUNT MATRIX DIMENSIONS\n")
cat("====================================================\n")
print(dim(count_mat_asv))

############################################################
# 3) BUILD EXPLICIT GENUS COUNT MATRIX
############################################################

# Mantener solo ASVs con género anotado
tax_df_genus <- tax_df %>%
  filter(!is.na(Genus), Genus != "")

if (nrow(tax_df_genus) == 0) {
  stop("No ASVs with assigned genus were found in the ITS taxonomy table.")
}

count_df_long <- as.data.frame(count_mat_asv) %>%
  rownames_to_column("ASV_ID") %>%
  pivot_longer(
    cols = -ASV_ID,
    names_to = "SampleID",
    values_to = "Count"
  ) %>%
  inner_join(
    tax_df_genus %>%
      select(ASV_ID, Kingdom, Phylum, Class, Order, Family, Genus, Species),
    by = "ASV_ID"
  )

# Agregar explícitamente por género
genus_count_long <- count_df_long %>%
  group_by(SampleID, Kingdom, Phylum, Class, Order, Family, Genus) %>%
  summarise(
    Count = sum(Count, na.rm = TRUE),
    .groups = "drop"
  )

# Tabla taxonómica estable por género
pick_first_nonempty <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & x != ""]
  if (length(x) == 0) return(NA_character_)
  x[1]
}

genus_tax_table <- genus_count_long %>%
  group_by(Genus) %>%
  summarise(
    Kingdom = pick_first_nonempty(Kingdom),
    Phylum  = pick_first_nonempty(Phylum),
    Class   = pick_first_nonempty(Class),
    Order   = pick_first_nonempty(Order),
    Family  = pick_first_nonempty(Family),
    .groups = "drop"
  )

# Matriz de counts género x muestra
genus_count_mat <- genus_count_long %>%
  select(Genus, SampleID, Count) %>%
  pivot_wider(
    names_from = SampleID,
    values_from = Count,
    values_fill = 0
  ) %>%
  as.data.frame()

rownames(genus_count_mat) <- genus_count_mat$Genus
genus_count_mat$Genus <- NULL
genus_count_mat <- as.matrix(genus_count_mat)

mode(genus_count_mat) <- "numeric"

cat("\n====================================================\n")
cat("GENUS COUNT MATRIX DIMENSIONS\n")
cat("====================================================\n")
print(dim(genus_count_mat))

############################################################
# 4) HELPER FUNCTIONS
############################################################

prepare_subset_for_deseq <- function(
    md,
    compartment_name,
    healthy_label,
    diseased_label
) {
  keep_ids <- md$SampleID[
    md$Compartment == compartment_name &
      md$Health_status %in% c(healthy_label, diseased_label)
  ]
  
  md_sub <- md %>%
    filter(SampleID %in% keep_ids) %>%
    mutate(
      Condition = ifelse(
        Health_status == healthy_label,
        "Healthy",
        "Diseased"
      ),
      Condition = factor(Condition, levels = c("Healthy", "Diseased"))
    )
  
  md_sub
}

run_deseq2_one_comparison <- function(
    count_matrix_genus,
    genus_tax_df,
    md,
    compartment_name,
    healthy_label,
    diseased_label,
    min_total_counts = 10,
    alpha = 0.05
) {
  md_sub <- prepare_subset_for_deseq(
    md = md,
    compartment_name = compartment_name,
    healthy_label = healthy_label,
    diseased_label = diseased_label
  )
  
  sample_ids <- md_sub$SampleID
  
  count_sub <- count_matrix_genus[, sample_ids, drop = FALSE]
  
  keep_genera <- rowSums(count_sub) >= min_total_counts
  count_sub <- count_sub[keep_genera, , drop = FALSE]
  
  if (nrow(count_sub) == 0) {
    stop(
      paste0(
        "No genera passed the minimum count filter in compartment: ",
        compartment_name
      )
    )
  }
  
  coldata <- md_sub %>%
    as.data.frame()
  rownames(coldata) <- coldata$SampleID
  coldata <- coldata[colnames(count_sub), , drop = FALSE]
  
  dds <- DESeqDataSetFromMatrix(
    countData = round(count_sub),
    colData = coldata,
    design = ~ Condition
  )
  
  dds <- dds[rowSums(counts(dds)) > 0, ]
  dds <- DESeq(dds, fitType = "parametric")
  
  res <- results(
    dds,
    contrast = c("Condition", "Diseased", "Healthy"),
    alpha = alpha
  )
  
  res_df <- as.data.frame(res) %>%
    rownames_to_column("Genus") %>%
    as_tibble() %>%
    left_join(genus_tax_df, by = "Genus") %>%
    mutate(
      Taxonomic_level = "Genus",
      Month = "October",
      Compartment = compartment_name,
      Comparison = paste0(healthy_label, " vs ", diseased_label),
      Direction = case_when(
        !is.na(padj) & padj < alpha & log2FoldChange > 0 ~ "Diseased > Healthy",
        !is.na(padj) & padj < alpha & log2FoldChange < 0 ~ "Healthy > Diseased",
        TRUE ~ "No significant difference"
      )
    ) %>%
    select(
      Taxonomic_level,
      Month,
      Compartment,
      Comparison,
      Genus,
      Kingdom,
      Phylum,
      Class,
      Order,
      Family,
      baseMean,
      log2FoldChange,
      lfcSE,
      stat,
      pvalue,
      padj,
      Direction
    ) %>%
    arrange(padj, desc(abs(log2FoldChange)))
  
  sig_res_df <- res_df %>%
    filter(!is.na(padj), padj < alpha) %>%
    arrange(padj, desc(abs(log2FoldChange)))
  
  list(
    dds = dds,
    full = res_df,
    sig = sig_res_df,
    md_sub = md_sub,
    count_sub = count_sub
  )
}

############################################################
# 5) DEFINE COMPARISONS
############################################################

comparisons_tbl <- tibble(
  Compartment = c("Rhizosphere", "Endosphere", "Phyllosphere"),
  Group_Healthy = c(
    "Healthy Rhizosphere",
    "Healthy Endosphere",
    "Healthy Phyllosphere"
  ),
  Group_Diseased = c(
    "Diseased Rhizosphere",
    "Diseased Endosphere",
    "Diseased Phyllosphere"
  ),
  short_label = c(
    "October_Rhizosphere",
    "October_Endosphere",
    "October_Phyllosphere"
  )
)

############################################################
# 6) RUN DESEQ2 FOR EACH COMPARISON
############################################################

all_full_results <- list()
all_sig_results <- list()
summary_counts <- list()

for (i in seq_len(nrow(comparisons_tbl))) {
  row_i <- comparisons_tbl[i, ]
  
  cat("\n====================================================\n")
  cat("RUNNING DESEQ2 — GENUS LEVEL:", row_i$short_label, "\n")
  cat("====================================================\n")
  
  res_i <- run_deseq2_one_comparison(
    count_matrix_genus = genus_count_mat,
    genus_tax_df = genus_tax_table,
    md = md0,
    compartment_name = row_i$Compartment,
    healthy_label = row_i$Group_Healthy,
    diseased_label = row_i$Group_Diseased,
    min_total_counts = min_total_counts_genus,
    alpha = alpha_sig
  )
  
  all_full_results[[row_i$short_label]] <- res_i$full
  all_sig_results[[row_i$short_label]] <- res_i$sig
  
  summary_counts[[row_i$short_label]] <- tibble(
    Month = "October",
    Compartment = row_i$Compartment,
    Comparison = paste0(row_i$Group_Healthy, " vs ", row_i$Group_Diseased),
    n_samples = nrow(res_i$md_sub),
    n_healthy = sum(res_i$md_sub$Condition == "Healthy"),
    n_diseased = sum(res_i$md_sub$Condition == "Diseased"),
    n_genera_tested = nrow(res_i$full),
    n_genera_significant = nrow(res_i$sig),
    n_diseased_higher = sum(res_i$sig$Direction == "Diseased > Healthy", na.rm = TRUE),
    n_healthy_higher = sum(res_i$sig$Direction == "Healthy > Diseased", na.rm = TRUE)
  )
  
  write_csv(
    res_i$full,
    file.path(outdir_deseq, paste0("01_DESeq2_Genus_", row_i$short_label, "_full.csv"))
  )
  
  write_csv(
    res_i$sig,
    file.path(outdir_deseq, paste0("02_DESeq2_Genus_", row_i$short_label, "_significant.csv"))
  )
  
  sig_compact <- res_i$sig %>%
    mutate(
      baseMean = round(baseMean, 3),
      log2FoldChange = round(log2FoldChange, 3),
      lfcSE = round(lfcSE, 3),
      stat = round(stat, 3),
      pvalue = signif(pvalue, 3),
      padj = signif(padj, 3)
    )
  
  write_csv(
    sig_compact,
    file.path(outdir_deseq, paste0("03_DESeq2_Genus_", row_i$short_label, "_significant_compact.csv"))
  )
}

############################################################
# 7) MERGE ALL RESULTS
############################################################

deseq_genus_full_all <- bind_rows(all_full_results)
deseq_genus_sig_all <- bind_rows(all_sig_results)
deseq_genus_summary <- bind_rows(summary_counts)

write_csv(
  deseq_genus_full_all,
  file.path(outdir_deseq, "04_DESeq2_Genus_all_comparisons_full.csv")
)

write_csv(
  deseq_genus_sig_all,
  file.path(outdir_deseq, "05_DESeq2_Genus_all_comparisons_significant.csv")
)

write_csv(
  deseq_genus_summary,
  file.path(outdir_deseq, "06_DESeq2_Genus_summary_counts.csv")
)

############################################################
# 8) OPTIONAL WRITING TABLE
############################################################

deseq_genus_writing <- deseq_genus_sig_all %>%
  transmute(
    Month,
    Compartment,
    Comparison,
    Genus,
    Family,
    Order,
    Class,
    Phylum,
    log2FoldChange = round(log2FoldChange, 3),
    padj = signif(padj, 3),
    Direction
  ) %>%
  arrange(Compartment, padj, desc(abs(log2FoldChange)))

write_csv(
  deseq_genus_writing,
  file.path(outdir_deseq, "07_DESeq2_Genus_writing_table.csv")
)

############################################################
# 9) CONSOLE OUTPUT
############################################################

cat("\n====================================================\n")
cat("DESEQ2 GENUS SUMMARY COUNTS\n")
cat("====================================================\n")
print(deseq_genus_summary, n = Inf, width = Inf)

cat("\n====================================================\n")
cat("DESEQ2 GENUS SIGNIFICANT RESULTS — HEAD\n")
cat("====================================================\n")
print(head(deseq_genus_sig_all, 30), width = Inf)

cat("\n====================================================\n")
cat("FILES WRITTEN TO:\n")
cat("====================================================\n")
print(normalizePath(outdir_deseq))



############################################################
# ITS — DESeq2 GENUS FIGURES
# mismo estilo que en 16S
############################################################

rm(list = ls())
graphics.off()

setwd("path/to/QIIME2_files_to_R_ITS")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(stringr)
  library(forcats)
  library(patchwork)
})

############################################################
# 0) INPUT / OUTPUT
############################################################

indir_deseq <- "04_differential_abundance_ITS_DESeq2_genus"
outdir_fig_genus <- "05_DESeq2_figures_Genus_ITS"
dir.create(outdir_fig_genus, recursive = TRUE, showWarnings = FALSE)

TOP_N_PER_DIRECTION <- 15

comparisons_tbl <- tibble::tribble(
  ~Month,    ~Compartment,   ~Group_Healthy,            ~Group_Diseased,             ~short_label,
  "October", "Rhizosphere",  "Healthy Rhizosphere",     "Diseased Rhizosphere",      "October_Rhizosphere",
  "October", "Endosphere",   "Healthy Endosphere",      "Diseased Endosphere",       "October_Endosphere",
  "October", "Phyllosphere", "Healthy Phyllosphere",    "Diseased Phyllosphere",     "October_Phyllosphere"
)

############################################################
# 1) HELPER FUNCTIONS
############################################################

make_genus_label <- function(x) {
  x <- as.character(x)
  x[is.na(x) | x == ""] <- "Unclassified genus"
  x
}

save_plot_both <- function(plot_obj, base_filename, width = 11, height = 8) {
  ggsave(
    file.path(outdir_fig_genus, paste0(base_filename, ".png")),
    plot = plot_obj,
    width = width,
    height = height,
    dpi = 600
  )
  
  ggsave(
    file.path(outdir_fig_genus, paste0(base_filename, ".tiff")),
    plot = plot_obj,
    width = width,
    height = height,
    dpi = 600,
    compression = "lzw"
  )
}

read_one_deseq_file <- function(short_label) {
  f <- file.path(
    indir_deseq,
    paste0("02_DESeq2_Genus_", short_label, "_significant.csv")
  )
  
  if (!file.exists(f)) {
    stop("File not found: ", f)
  }
  
  read_csv(f, show_col_types = FALSE)
}

make_deseq_genus_figure <- function(
    month_name,
    compartment_name,
    healthy_label,
    diseased_label,
    short_label,
    top_n_per_direction = 15
) {
  
  df <- read_one_deseq_file(short_label) %>%
    mutate(
      Genus_label = make_genus_label(Genus),
      minus_log10_padj = -log10(padj),
      Direction = factor(
        Direction,
        levels = c("Healthy > Diseased", "Diseased > Healthy")
      )
    ) %>%
    filter(Direction %in% c("Healthy > Diseased", "Diseased > Healthy"))
  
  if (nrow(df) == 0) {
    return(
      ggplot() +
        annotate("text", x = 0, y = 0, label = paste("No significant genera:", short_label), size = 6) +
        theme_void()
    )
  }
  
  top_df <- df %>%
    group_by(Direction) %>%
    arrange(padj, desc(abs(log2FoldChange))) %>%
    slice_head(n = top_n_per_direction) %>%
    ungroup()
  
  top_df <- top_df %>%
    mutate(
      plot_label = paste0(Genus_label),
      plot_label = forcats::fct_reorder(plot_label, log2FoldChange)
    )
  
  p_bar <- ggplot(
    top_df,
    aes(x = plot_label, y = log2FoldChange, fill = Direction)
  ) +
    geom_col(color = "black", linewidth = 0.25) +
    coord_flip() +
    scale_fill_manual(
      values = c(
        "Healthy > Diseased" = "#4B9EA7",
        "Diseased > Healthy" = "#801723"
      )
    ) +
    labs(
      title = paste0(month_name, " — ", compartment_name),
      subtitle = paste0(
        "Top ", top_n_per_direction, " enriched in healthy and top ",
        top_n_per_direction, " enriched in diseased"
      ),
      x = NULL,
      y = "log2 fold change (Diseased vs Healthy)"
    ) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      legend.title = element_blank(),
      plot.title = element_text(face = "bold"),
      axis.text.y = element_text(size = 10)
    )
  
  p_dot <- ggplot(
    top_df,
    aes(x = log2FoldChange, y = minus_log10_padj, color = Direction)
  ) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
    geom_point(size = 3, alpha = 0.9) +
    geom_text(
      aes(label = Genus_label),
      size = 3,
      vjust = -0.4,
      check_overlap = TRUE
    ) +
    scale_color_manual(
      values = c(
        "Healthy > Diseased" = "#4B9EA7",
        "Diseased > Healthy" = "#801723"
      )
    ) +
    labs(
      x = "log2 fold change (Diseased vs Healthy)",
      y = expression(-log[10](adjusted~p))
    ) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      legend.title = element_blank()
    )
  
  p_bar / p_dot
}

############################################################
# 2) MAKE ALL FIGURES
############################################################

plot_list <- list()

for (i in seq_len(nrow(comparisons_tbl))) {
  row_i <- comparisons_tbl[i, ]
  
  fig_i <- make_deseq_genus_figure(
    month_name = row_i$Month,
    compartment_name = row_i$Compartment,
    healthy_label = row_i$Group_Healthy,
    diseased_label = row_i$Group_Diseased,
    short_label = row_i$short_label,
    top_n_per_direction = TOP_N_PER_DIRECTION
  )
  
  plot_list[[row_i$short_label]] <- fig_i
  print(fig_i)
  
  save_plot_both(
    fig_i,
    paste0("ITS_DESeq2_Genus_", row_i$short_label),
    width = 11,
    height = 8
  )
}

############################################################
# 3) OPTIONAL: COMBINED PANEL
############################################################

combined_panel <- wrap_plots(plot_list, ncol = 1)
print(combined_panel)

save_plot_both(
  combined_panel,
  "ITS_DESeq2_Genus_all_comparisons_combined",
  width = 12,
  height = 8 * length(plot_list)
)

############################################################
# 4) OUTPUT CHECK
############################################################

cat("\n====================================================\n")
cat("FILES WRITTEN TO:\n")
cat("====================================================\n")
print(normalizePath(outdir_fig_genus))



setwd("path/to/QIIME2_files_to_R_ITS")

library(readr)
library(dplyr)

sig_countsITS <- read_csv(
  "04_differential_abundance_ITS_DESeq2_genus/06_DESeq2_Genus_summary_counts.csv",
  show_col_types = FALSE
)

deseqITS <- read_csv(
  "04_differential_abundance_ITS_DESeq2_genus/07_DESeq2_Genus_writing_table.csv",
  show_col_types = FALSE
)



cat("\n=== ITS SUMMARY COUNTS ===\n")
print(sig_countsITS, n = Inf, width = Inf)

cat("\n=== ITS TOP 10 PER COMPARTMENT / DIRECTION ===\n")
deseqITS %>%
  group_by(Compartment, Direction) %>%
  arrange(padj, desc(abs(log2FoldChange))) %>%
  slice_head(n = 10) %>%
  ungroup() %>%
  select(Month, Compartment, Direction, Genus, Family, Phylum, log2FoldChange, padj) %>%
  arrange(Compartment, Direction, padj) %>%
  print(n = Inf, width = Inf)

cat("\n=== ITS RECURRENT GENERA n>=2 COMPARTMENTS ===\n")
deseqITS %>%
  group_by(Direction, Genus, Family, Phylum) %>%
  summarise(
    n_compartments = n_distinct(Compartment),
    compartments = paste(unique(Compartment), collapse = "; "),
    min_padj = min(padj, na.rm = TRUE),
    max_abs_log2FC = max(abs(log2FoldChange), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(n_compartments >= 2) %>%
  arrange(Direction, desc(n_compartments), min_padj) %>%
  print(n = Inf, width = Inf)



















###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################







############################################################
# ITS — Botryosphaeriaceae abundance in phyllosphere
# Sample-level plot (Healthy vs Diseased)
############################################################

rm(list = ls())
graphics.off()

setwd("path/to/QIIME2_files_to_R_ITS")

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(forcats)
  library(stringr)
})

############################################################
# 0) INPUT / OUTPUT
############################################################

infile <- "03_taxonomy_ITS/tables/14_botryosphaeriaceae_sample_level_ITS.csv"

outdir <- "03_taxonomy_ITS/figures"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

hs_colors <- c(
  "Bulk soil"              = "#896054",
  "Healthy Rhizosphere"    = "#4B9EA7",
  "Diseased Rhizosphere"   = "#801723",
  "Healthy Endosphere"     = "#6B9863",
  "Diseased Endosphere"    = "#E7C313",
  "Healthy Phyllosphere"   = "#4A5263",
  "Diseased Phyllosphere"  = "#6B4EA7"
)

############################################################
# 1) LOAD DATA
############################################################

botryo <- read_csv(infile, show_col_types = FALSE) %>%
  mutate(
    Month = as.character(Month),
    Compartment = as.character(Compartment),
    Health_status = as.character(Health_status),
    SampleID = as.character(SampleID)
  ) %>%
  filter(
    Compartment == "Phyllosphere",
    Health_status %in% c("Healthy Phyllosphere", "Diseased Phyllosphere")
  )

cat("\n====================================================\n")
cat("BOTRYOSPHAERIACEAE SAMPLE-LEVEL TABLE\n")
cat("====================================================\n")
print(botryo, n = Inf, width = Inf)

############################################################
# 2) ORDER SAMPLES
############################################################

sample_order <- botryo %>%
  arrange(Health_status, SampleID) %>%
  pull(SampleID)

botryo <- botryo %>%
  mutate(
    Health_status = factor(
      Health_status,
      levels = c("Healthy Phyllosphere", "Diseased Phyllosphere")
    ),
    SampleID = factor(SampleID, levels = sample_order)
  )

############################################################
# 3) SUMMARY FOR MEAN ± SE
############################################################

sumdf <- botryo %>%
  group_by(Health_status) %>%
  summarise(
    n = n(),
    mean_abund = mean(Relative_abundance_pct, na.rm = TRUE),
    sd_abund = sd(Relative_abundance_pct, na.rm = TRUE),
    se_abund = sd_abund / sqrt(n),
    .groups = "drop"
  )

############################################################
# 4) PLOT 1 — SAMPLE-LEVEL DOT PLOT
############################################################

p_dots <- ggplot(
  botryo,
  aes(x = SampleID, y = Relative_abundance_pct, fill = Health_status)
) +
  geom_col(color = "black", linewidth = 0.25, width = 0.8) +
  scale_fill_manual(values = hs_colors, drop = FALSE) +
  labs(
    title = "Botryosphaeriaceae abundance across phyllosphere samples",
    subtitle = "ITS — October",
    x = NULL,
    y = "Relative abundance (%)"
  ) +
  facet_grid(. ~ Health_status, scales = "free_x", space = "free_x") +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold", size = 12),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 9),
    axis.title.y = element_text(size = 14, face = "bold"),
    plot.title = element_text(face = "bold")
  )

print(p_dots)

ggsave(
  file.path(outdir, "24_botryosphaeriaceae_phyllosphere_sample_level.png"),
  plot = p_dots,
  width = 12,
  height = 5.5,
  dpi = 600
)

ggsave(
  file.path(outdir, "24_botryosphaeriaceae_phyllosphere_sample_level.tiff"),
  plot = p_dots,
  width = 12,
  height = 5.5,
  dpi = 600,
  compression = "lzw"
)

############################################################
# 5) PLOT 2 — GROUPED DOT PLOT WITH MEAN ± SE
############################################################

p_group <- ggplot(
  botryo,
  aes(x = Health_status, y = Relative_abundance_pct, fill = Health_status)
) +
  geom_violin(trim = FALSE, linewidth = 0.6, alpha = 0.9) +
  geom_jitter(
    width = 0.12,
    size = 2.2,
    shape = 21,
    stroke = 0.4,
    color = "black",
    alpha = 0.9
  ) +
  geom_errorbar(
    data = sumdf,
    aes(
      x = Health_status,
      y = mean_abund,
      ymin = mean_abund - se_abund,
      ymax = mean_abund + se_abund
    ),
    width = 0.12,
    linewidth = 0.7,
    inherit.aes = FALSE
  ) +
  geom_point(
    data = sumdf,
    aes(x = Health_status, y = mean_abund),
    size = 3,
    shape = 23,
    fill = "white",
    color = "black",
    stroke = 0.6,
    inherit.aes = FALSE
  ) +
  scale_fill_manual(values = hs_colors, drop = FALSE) +
  labs(
    title = "Botryosphaeriaceae abundance in phyllosphere",
    subtitle = "ITS — October",
    x = NULL,
    y = "Relative abundance (%)"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 11),
    axis.title.y = element_text(size = 14, face = "bold"),
    plot.title = element_text(face = "bold")
  )

print(p_group)

ggsave(
  file.path(outdir, "25_botryosphaeriaceae_phyllosphere_grouped.png"),
  plot = p_group,
  width = 6,
  height = 5.5,
  dpi = 600
)

ggsave(
  file.path(outdir, "25_botryosphaeriaceae_phyllosphere_grouped.tiff"),
  plot = p_group,
  width = 6,
  height = 5.5,
  dpi = 600,
  compression = "lzw"
)

############################################################
# 6) EXPORT CLEAN TABLE FOR WRITING
############################################################

write_csv(
  botryo,
  file.path(outdir, "24a_botryosphaeriaceae_phyllosphere_sample_level_table.csv")
)

write_csv(
  sumdf,
  file.path(outdir, "25a_botryosphaeriaceae_phyllosphere_summary_table.csv")
)

cat("\n====================================================\n")
cat("SUMMARY\n")
cat("====================================================\n")
print(sumdf, n = Inf, width = Inf)

cat("\n====================================================\n")
cat("FILES WRITTEN TO:\n")
cat("====================================================\n")
print(normalizePath(outdir))









############################################################
# ITS — Neofusicoccum abundance in phyllosphere
# Sample-level plot (Healthy vs Diseased)
############################################################

rm(list = ls())
graphics.off()

setwd("path/to/QIIME2_files_to_R_ITS")

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(forcats)
  library(stringr)
})

############################################################
# 0) INPUT / OUTPUT
############################################################

infile <- "03_taxonomy_ITS/tables/15_neofusicoccum_sample_level_ITS.csv"

outdir <- "03_taxonomy_ITS/figures"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

hs_colors <- c(
  "Bulk soil"              = "#896054",
  "Healthy Rhizosphere"    = "#4B9EA7",
  "Diseased Rhizosphere"   = "#801723",
  "Healthy Endosphere"     = "#6B9863",
  "Diseased Endosphere"    = "#E7C313",
  "Healthy Phyllosphere"   = "#4A5263",
  "Diseased Phyllosphere"  = "#6B4EA7"
)

############################################################
# 1) LOAD DATA
############################################################

neofu <- read_csv(infile, show_col_types = FALSE) %>%
  mutate(
    Month = as.character(Month),
    Compartment = as.character(Compartment),
    Health_status = as.character(Health_status),
    SampleID = as.character(SampleID)
  ) %>%
  filter(
    Compartment == "Phyllosphere",
    Health_status %in% c("Healthy Phyllosphere", "Diseased Phyllosphere")
  )

cat("\n====================================================\n")
cat("NEOFUSICOCCUM SAMPLE-LEVEL TABLE\n")
cat("====================================================\n")
print(neofu, n = Inf, width = Inf)

############################################################
# 2) ORDER SAMPLES
############################################################

sample_order <- neofu %>%
  arrange(Health_status, SampleID) %>%
  pull(SampleID)

neofu <- neofu %>%
  mutate(
    Health_status = factor(
      Health_status,
      levels = c("Healthy Phyllosphere", "Diseased Phyllosphere")
    ),
    SampleID = factor(SampleID, levels = sample_order)
  )

############################################################
# 3) SUMMARY FOR MEAN ± SE
############################################################

sumdf <- neofu %>%
  group_by(Health_status) %>%
  summarise(
    n = n(),
    mean_abund = mean(Relative_abundance_pct, na.rm = TRUE),
    sd_abund = sd(Relative_abundance_pct, na.rm = TRUE),
    se_abund = sd_abund / sqrt(n),
    .groups = "drop"
  )

############################################################
# 4) PLOT 1 — SAMPLE-LEVEL DOT/BAR PLOT
############################################################

p_dots <- ggplot(
  neofu,
  aes(x = SampleID, y = Relative_abundance_pct, fill = Health_status)
) +
  geom_col(color = "black", linewidth = 0.25, width = 0.8) +
  scale_fill_manual(values = hs_colors, drop = FALSE) +
  labs(
    title = "Neofusicoccum abundance across phyllosphere samples",
    subtitle = "ITS — October",
    x = NULL,
    y = "Relative abundance (%)"
  ) +
  facet_grid(. ~ Health_status, scales = "free_x", space = "free_x") +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold", size = 12),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 9),
    axis.title.y = element_text(size = 14, face = "bold"),
    plot.title = element_text(face = "bold")
  )

print(p_dots)

ggsave(
  file.path(outdir, "26_neofusicoccum_phyllosphere_sample_level.png"),
  plot = p_dots,
  width = 12,
  height = 5.5,
  dpi = 600
)

ggsave(
  file.path(outdir, "26_neofusicoccum_phyllosphere_sample_level.tiff"),
  plot = p_dots,
  width = 12,
  height = 5.5,
  dpi = 600,
  compression = "lzw"
)

############################################################
# 5) PLOT 2 — GROUPED VIOLIN + DOTS + MEAN ± SE
############################################################

p_group <- ggplot(
  neofu,
  aes(x = Health_status, y = Relative_abundance_pct, fill = Health_status)
) +
  geom_violin(trim = FALSE, linewidth = 0.6, alpha = 0.9) +
  geom_jitter(
    width = 0.12,
    size = 2.2,
    shape = 21,
    stroke = 0.4,
    color = "black",
    alpha = 0.9
  ) +
  geom_errorbar(
    data = sumdf,
    aes(
      x = Health_status,
      y = mean_abund,
      ymin = mean_abund - se_abund,
      ymax = mean_abund + se_abund
    ),
    width = 0.12,
    linewidth = 0.7,
    inherit.aes = FALSE
  ) +
  geom_point(
    data = sumdf,
    aes(x = Health_status, y = mean_abund),
    size = 3,
    shape = 23,
    fill = "white",
    color = "black",
    stroke = 0.6,
    inherit.aes = FALSE
  ) +
  scale_fill_manual(values = hs_colors, drop = FALSE) +
  labs(
    title = "Neofusicoccum abundance in phyllosphere",
    subtitle = "ITS — October",
    x = NULL,
    y = "Relative abundance (%)"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 11),
    axis.title.y = element_text(size = 14, face = "bold"),
    plot.title = element_text(face = "bold")
  )

print(p_group)

ggsave(
  file.path(outdir, "27_neofusicoccum_phyllosphere_grouped.png"),
  plot = p_group,
  width = 6,
  height = 5.5,
  dpi = 600
)

ggsave(
  file.path(outdir, "27_neofusicoccum_phyllosphere_grouped.tiff"),
  plot = p_group,
  width = 6,
  height = 5.5,
  dpi = 600,
  compression = "lzw"
)

############################################################
# 6) EXPORT CLEAN TABLES
############################################################

write_csv(
  neofu,
  file.path(outdir, "26a_neofusicoccum_phyllosphere_sample_level_table.csv")
)

write_csv(
  sumdf,
  file.path(outdir, "27a_neofusicoccum_phyllosphere_summary_table.csv")
)

cat("\n====================================================\n")
cat("SUMMARY\n")
cat("====================================================\n")
print(sumdf, n = Inf, width = Inf)

cat("\n====================================================\n")
cat("FILES WRITTEN TO:\n")
cat("====================================================\n")
print(normalizePath(outdir))


















############################################################
# SECTION 3B — TREE SUBSAMPLING EFFECT CHECK (ITS, OCTOBER ONLY)
#
# Objetivos:
# 1) Evaluar si las submuestras de un mismo árbol se agrupan entre sí
# 2) Generar PCoAs Bray-Curtis coloreados por Health_status y por árbol
# 3) Testar el efecto de Tree dentro de cada Compartment x Health_status
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

outdir_treecheck <- "02B_tree_subsampling_check_ITS"
outdir_treecheck_figs <- file.path(outdir_treecheck, "figures")
outdir_treecheck_tables <- file.path(outdir_treecheck, "tables")

dir.create(outdir_treecheck, recursive = TRUE, showWarnings = FALSE)
dir.create(outdir_treecheck_figs, recursive = TRUE, showWarnings = FALSE)
dir.create(outdir_treecheck_tables, recursive = TRUE, showWarnings = FALSE)

############################################################
# 0B) LOAD PHYLOSEQ IF ps0 DOES NOT EXIST
############################################################

if (!exists("ps0")) {
  
  cat("\nps0 no existe en la sesión. Cargando Phyloseq_ITS.rds...\n")
  
  physeq_path <- "Phyloseq_ITS.rds"
  
  if (!file.exists(physeq_path)) {
    stop(
      "No encuentro el archivo: ", physeq_path,
      "\nComprueba que estás en el directorio correcto con getwd()."
    )
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

meta_tree$Month <- factor(as.character(meta_tree$Month), levels = c("October"))
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
cat("ITS TREE SUBSAMPLING CHECK — BASIC CHECKS\n")
cat("====================================================\n")
cat("Metadata columns:\n")
print(colnames(meta_tree))

cat("\nTrees per Compartment x Health_status:\n")
print(
  meta_tree %>%
    filter(Compartment != "Bulk soil") %>%
    group_by(Compartment, Health_status) %>%
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
# 3) PCoA BY COMPARTMENT — SAMPLE LEVEL
############################################################

run_treecheck_pcoa_one <- function(ps_in, compartment_name) {
  
  md <- as(sample_data(ps_in), "data.frame")
  
  keep_ids <- rownames(md)[
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
    warning("Skipping ", compartment_name, ": not enough samples or trees.")
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
      title = paste0("ITS Bray-Curtis PCoA — October · ", compartment_name),
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
      title = paste0("ITS Bray-Curtis PCoA — October · ", compartment_name),
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
    "ITS_October_",
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
    compartment = compartment_name,
    ps = ps_sub,
    metadata = md_sub,
    distance = D_bray,
    pcoa_df = pcoa_df,
    p_health_tree = p_health_tree,
    p_tree_facet = p_tree_facet
  )
}

compartments_treecheck <- c("Rhizosphere", "Endosphere", "Phyllosphere")

treecheck_pcoa_results <- list()

for (co in compartments_treecheck) {
  treecheck_pcoa_results[[co]] <- run_treecheck_pcoa_one(ps0, co)
}

############################################################
# 4) TEST TREE EFFECT WITHIN EACH COMPARTMENT x HEALTH_STATUS
############################################################

test_tree_effect_one_group <- function(ps_in, compartment_name, health_status_name) {
  
  md <- as(sample_data(ps_in), "data.frame")
  
  keep_ids <- rownames(md)[
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
      Month = "October",
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
    Month = "October",
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
  compartments_treecheck,
  function(co) {
    relevant_health <- health_levels_treecheck[
      grepl(co, health_levels_treecheck)
    ]
    purrr::map_dfr(
      relevant_health,
      function(hs) test_tree_effect_one_group(ps0, co, hs)
    )
  }
)

write_csv(
  tree_effect_results,
  file.path(outdir_treecheck_tables, "01_tree_effect_within_compartment_health_ITS.csv")
)

############################################################
# 5) WITHIN-TREE VS BETWEEN-TREE DISTANCES
############################################################

compare_within_between_tree_distances_one_group <- function(ps_in, compartment_name, health_status_name) {
  
  md <- as(sample_data(ps_in), "data.frame")
  
  keep_ids <- rownames(md)[
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
      Month = "October",
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
      Month = "October",
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
  compartments_treecheck,
  function(co) {
    relevant_health <- health_levels_treecheck[
      grepl(co, health_levels_treecheck)
    ]
    purrr::map_dfr(
      relevant_health,
      function(hs) compare_within_between_tree_distances_one_group(ps0, co, hs)
    )
  }
)

write_csv(
  within_between_distance_results,
  file.path(outdir_treecheck_tables, "02_within_tree_vs_between_tree_distances_ITS.csv")
)

############################################################
# 6) TREE-LEVEL AGGREGATION ANALYSIS
#
# Suma las cuentas de las submuestras del mismo árbol
# dentro de Compartment x Health_status x Trees.
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
  
  ps_agg <- prune_samples(sample_sums(ps_agg) > 0, ps_agg)
  ps_agg <- prune_taxa(taxa_sums(ps_agg) > 0, ps_agg)
  
  ps_agg
}

ps0_tree_agg <- aggregate_phyloseq_by_tree(ps0)

saveRDS(
  ps0_tree_agg,
  file.path(outdir_treecheck, "Phyloseq_ITS_aggregated_by_tree.rds")
)

############################################################
# 7) PCoA AND PERMANOVA AT TREE-LEVEL
############################################################

run_treelevel_beta_one <- function(ps_in, compartment_name) {
  
  md <- as(sample_data(ps_in), "data.frame")
  
  keep_ids <- rownames(md)[
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
      Month = "October",
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
      title = paste0("ITS tree-level Bray-Curtis PCoA — October · ", compartment_name),
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
    "ITS_October_",
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
    Month = "October",
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
  compartments_treecheck,
  function(co) run_treelevel_beta_one(ps0_tree_agg, co)
)

write_csv(
  treelevel_health_results,
  file.path(outdir_treecheck_tables, "03_tree_level_aggregated_health_effect_ITS.csv")
)

############################################################
# 8) PRINT SUMMARY TO CONSOLE
############################################################

cat("\n====================================================\n")
cat("ITS — TREE EFFECT WITHIN COMPARTMENT x HEALTH_STATUS\n")
cat("====================================================\n")
print(tree_effect_results, n = Inf, width = Inf)

cat("\n====================================================\n")
cat("ITS — WITHIN-TREE VS BETWEEN-TREE DISTANCES\n")
cat("====================================================\n")
print(within_between_distance_results, n = Inf, width = Inf)

cat("\n====================================================\n")
cat("ITS — TREE-LEVEL AGGREGATED HEALTH EFFECT\n")
cat("====================================================\n")
print(treelevel_health_results, n = Inf, width = Inf)

cat("\n====================================================\n")
cat("FILES WRITTEN TO:\n")
cat("====================================================\n")
print(normalizePath(outdir_treecheck))
print(normalizePath(outdir_treecheck_figs))
print(normalizePath(outdir_treecheck_tables))





setwd("path/to/QIIME2_files_to_R_ITS")

library(readr)
library(dplyr)

betaITS_comp <- read_csv(
  "02_beta_diversity_ITS/05_beta_by_compartment_summary_ITS.csv",
  show_col_types = FALSE
)

betaITS_comp %>%
  select(
    Month,
    Compartment,
    n_samples,
    n_Healthy,
    n_Diseased,
    Bray_PERMANOVA_R2,
    Bray_PERMANOVA_F,
    Bray_PERMANOVA_p,
    Bray_betadisper_F,
    Bray_betadisper_p,
    Jaccard_PERMANOVA_R2,
    Jaccard_PERMANOVA_F,
    Jaccard_PERMANOVA_p,
    Jaccard_betadisper_F,
    Jaccard_betadisper_p
  ) %>%
  print(n = Inf, width = Inf)










############# Exploratorio
# ============================================================
# TREE-LEVEL ROBUSTNESS ANALYSIS - ITS
# Merge subsamples per biological tree into one tree-level profile
# Corrected for ITS Trees codes:
# A1RO, A2RO, A3RO, A4RO -> Tree_real = A
# ============================================================

setwd("path/to/QIIME2_files_to_R_ITS")

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

OUTDIR <- "02C_tree_level_aggregation_ITS"
FIGDIR <- file.path(OUTDIR, "figures")
TABDIR <- file.path(OUTDIR, "tables")

dir.create(OUTDIR, showWarnings = FALSE)
dir.create(FIGDIR, showWarnings = FALSE)
dir.create(TABDIR, showWarnings = FALSE)

# ------------------------------------------------------------
# 2. Load phyloseq object
# ------------------------------------------------------------

ps0 <- readRDS("Phyloseq_ITS.rds")

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
# 4. Create real biological tree ID
# ------------------------------------------------------------

# ITS Trees values are sub-sample IDs, for example:
# A1RO, A2RO, A3RO, A4RO.
# The real biological tree is the first letter: A, B, C, D, E, F, G, H, I.

md0 <- md0 %>%
  mutate(
    Tree_real = sub("^([A-Z]).*$", "\\1", Trees),
    Tree_level_ID = paste(Month, Compartment, Health_binary, Tree_real, sep = "__")
  )

md0_for_phyloseq <- as.data.frame(md0)
rownames(md0_for_phyloseq) <- md0_for_phyloseq$SampleID
sample_data(ps0) <- sample_data(md0_for_phyloseq)

# Check number of subsamples per real tree-level group
group_counts <- md0 %>%
  count(Month, Compartment, Health_binary, Tree_real, Tree_level_ID, name = "n_subsamples") %>%
  arrange(Month, Compartment, Health_binary, Tree_real)

write_csv(
  group_counts,
  file.path(TABDIR, "01_original_subsamples_per_tree_level_group_ITS.csv")
)

cat("\nSummary of n_subsamples per real tree-level group:\n")
utils::write.table(
  as.data.frame(table(as.data.frame(group_counts)$n_subsamples)),
  row.names = FALSE,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

cat("\nFirst rows of group_counts:\n")
utils::write.table(
  head(as.data.frame(group_counts), 50),
  row.names = FALSE,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

# ------------------------------------------------------------
# 5. Aggregate ASV counts by real tree-level group
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
    Tree_real = first(Tree_real),
    n_subsamples = n(),
    original_Trees_values = paste(Trees, collapse = "; "),
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

saveRDS(ps_tree, file.path(OUTDIR, "Phyloseq_ITS_tree_level_aggregated.rds"))

write_csv(
  md_tree %>% rownames_to_column("Tree_level_ID"),
  file.path(TABDIR, "02_tree_level_metadata_ITS.csv")
)

cat("\nTree-level phyloseq object:\n")
print(ps_tree)

cat("\nTree-level sample counts:\n")
print(table(sample_data(ps_tree)$Month, sample_data(ps_tree)$Compartment, sample_data(ps_tree)$Health_binary))

# ------------------------------------------------------------
# 7. Function: Bray PCoA + PERMANOVA + betadisper
# ------------------------------------------------------------

run_tree_level_beta_ITS <- function(ps_input, compartment_name) {
  
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
    select(Tree_level_ID, Month, Compartment, Health_binary, Tree_real, n_subsamples)
  
  coords <- coords %>%
    left_join(md_join, by = "Tree_level_ID")
  
  axis1 <- round(pcoa_res$values$Relative_eig[1] * 100, 1)
  axis2 <- round(pcoa_res$values$Relative_eig[2] * 100, 1)
  
  coords$PCoA1_var <- axis1
  coords$PCoA2_var <- axis2
  
  write_csv(
    coords,
    file.path(TABDIR, paste0("03_PCoA_coordinates_tree_level_", compartment_name, "_ITS.csv"))
  )
  
  # PCoA without ellipses
  p_no_ellipse <- ggplot(coords, aes(x = PCoA1, y = PCoA2, color = Health_binary)) +
    geom_point(size = 4, alpha = 0.9) +
    geom_text(aes(label = Tree_real), vjust = -0.8, size = 3, show.legend = FALSE) +
    labs(
      title = paste0("Tree-level Bray-Curtis PCoA - ", compartment_name, " - ITS"),
      x = paste0("PCoA1 (", axis1, "%)"),
      y = paste0("PCoA2 (", axis2, "%)"),
      color = "Health status"
    ) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5)
    )
  
  ggsave(
    file.path(FIGDIR, paste0("PCoA_tree_level_Bray_", compartment_name, "_ITS.pdf")),
    p_no_ellipse,
    width = 7,
    height = 5
  )
  
  ggsave(
    file.path(FIGDIR, paste0("PCoA_tree_level_Bray_", compartment_name, "_ITS.png")),
    p_no_ellipse,
    width = 7,
    height = 5,
    dpi = 300
  )
  
  # PCoA with health-status ellipses
  p_ellipse <- ggplot(coords, aes(x = PCoA1, y = PCoA2, color = Health_binary)) +
    geom_point(size = 4, alpha = 0.9) +
    geom_text(aes(label = Tree_real), vjust = -0.8, size = 3, show.legend = FALSE) +
    stat_ellipse(
      aes(group = Health_binary),
      type = "t",
      linetype = "dashed",
      linewidth = 0.6,
      show.legend = FALSE
    ) +
    labs(
      title = paste0("Tree-level Bray-Curtis PCoA with health ellipses - ", compartment_name, " - ITS"),
      x = paste0("PCoA1 (", axis1, "%)"),
      y = paste0("PCoA2 (", axis2, "%)"),
      color = "Health status"
    ) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5)
    )
  
  ggsave(
    file.path(FIGDIR, paste0("PCoA_tree_level_Bray_", compartment_name, "_ITS_with_health_ellipses.pdf")),
    p_ellipse,
    width = 7,
    height = 5
  )
  
  ggsave(
    file.path(FIGDIR, paste0("PCoA_tree_level_Bray_", compartment_name, "_ITS_with_health_ellipses.png")),
    p_ellipse,
    width = 7,
    height = 5,
    dpi = 300
  )
  
  # PERMANOVA: Health
  perm_health <- adonis2(dist_bray ~ Health_binary, data = md, permutations = 999)
  
  perm_health_df <- as.data.frame(perm_health) %>%
    rownames_to_column("term") %>%
    mutate(
      Compartment = compartment_name,
      Analysis = "Tree_level_Health",
      Distance = "Bray"
    ) %>%
    relocate(Analysis, Compartment, Distance, term)
  
  # betadisper: Health
  bd <- betadisper(dist_bray, md$Health_binary)
  bd_perm <- permutest(bd, permutations = 999)
  
  bd_dist <- data.frame(
    Tree_level_ID = names(bd$distances),
    distance_to_centroid = as.numeric(bd$distances)
  ) %>%
    left_join(
      md %>%
        select(Tree_level_ID, Month, Compartment, Health_binary, Tree_real),
      by = "Tree_level_ID"
    )
  
  bd_summary <- bd_dist %>%
    group_by(Compartment, Health_binary) %>%
    summarise(
      n = n(),
      mean_distance_to_centroid = mean(distance_to_centroid),
      sd_distance_to_centroid = sd(distance_to_centroid),
      .groups = "drop"
    )
  
  bd_result <- data.frame(
    Compartment = compartment_name,
    Distance = "Bray",
    betadisper_F = bd_perm$tab[1, "F"],
    betadisper_p = bd_perm$tab[1, "Pr(>F)"]
  )
  
  return(list(
    pcoa = coords,
    permanova_health = perm_health_df,
    betadisper_results = bd_result,
    betadisper_distances = bd_dist,
    betadisper_summary = bd_summary
  ))
}

# ------------------------------------------------------------
# 8. Run analyses by compartment
# ------------------------------------------------------------

compartments <- c("Rhizosphere", "Endosphere", "Phyllosphere")

results <- lapply(compartments, function(x) run_tree_level_beta_ITS(ps_tree, x))
names(results) <- compartments

permanova_health_all <- bind_rows(lapply(results, `[[`, "permanova_health"))
betadisper_results_all <- bind_rows(lapply(results, `[[`, "betadisper_results"))
betadisper_distances_all <- bind_rows(lapply(results, `[[`, "betadisper_distances"))
betadisper_summary_all <- bind_rows(lapply(results, `[[`, "betadisper_summary"))

# ------------------------------------------------------------
# 9. Save output tables
# ------------------------------------------------------------

write_csv(
  permanova_health_all,
  file.path(TABDIR, "04_PERMANOVA_tree_level_Health_ITS.csv")
)

write_csv(
  betadisper_results_all,
  file.path(TABDIR, "05_betadisper_tree_level_Health_ITS.csv")
)

write_csv(
  betadisper_distances_all,
  file.path(TABDIR, "06_betadisper_tree_level_distances_to_centroid_ITS.csv")
)

write_csv(
  betadisper_summary_all,
  file.path(TABDIR, "07_betadisper_tree_level_summary_by_health_ITS.csv")
)

# ------------------------------------------------------------
# 10. Print key outputs robustly
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("ITS TREE-LEVEL PERMANOVA\n")
cat("Model: Bray ~ Health_binary\n")
cat("============================================================\n")
utils::write.table(
  as.data.frame(permanova_health_all),
  row.names = FALSE,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

cat("\n============================================================\n")
cat("ITS TREE-LEVEL BETADISPER\n")
cat("============================================================\n")
utils::write.table(
  as.data.frame(betadisper_results_all),
  row.names = FALSE,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

cat("\n============================================================\n")
cat("ITS TREE-LEVEL DISTANCE TO CENTROID SUMMARY\n")
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

