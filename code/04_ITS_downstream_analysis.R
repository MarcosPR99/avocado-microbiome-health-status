#!/usr/bin/env Rscript

###############################################################################
# 04_ITS_downstream_analysis.R
#
# Downstream analyses of ITS amplicon data (October) for:
# "Fungal diseases reshape avocado-associated multi-kingdom microbiomes through
#  compartment-specific restructuring and network rewiring"
#
# Outputs:
#   Table S2  - sequencing output and taxonomic recovery (ITS part)
#   Table S3  - alpha diversity (ITS part)
#   Table S4  - beta diversity, Bray-Curtis and Jaccard (ITS part)
#   Table S6  - taxonomic composition (ITS part)
#   Table S7  - differentially abundant genera, DESeq2 (ITS part)
#   Table S8  - R. necatrix-compatible ITS signals
#   Table S9  - branch dieback-associated fungal signals
#   Table S5 - tree-level robustness analysis (ITS part)
#
# Inputs (working directory; see REPRODUCIBILITY_README.md):
#   Phyloseq_ITS.rds
#   metadata_ITS_analysis.tsv
#   alpha-rarefaction_ITS_maxdepth39000.qzv
###############################################################################

# ----------------------------------------------------------------------------
# Configuration
# ----------------------------------------------------------------------------

physeq_path <- "Phyloseq_ITS.rds"
metadata_file <- "metadata_ITS_analysis.tsv"
alpha_qzv_file <- "alpha-rarefaction_ITS_maxdepth39000.qzv"
alpha_depth <- 39000
outdir <- "results_ITS"

# R. necatrix-compatible ASVs: candidate Xylariaceae/Xylariales ASVs were queried
# with BLASTn (BLAST+ v2.14.0); ASVs whose best hit was R. necatrix (accession
# PP485889) with >= 98% identity and 100% query coverage were retained.
rnecatrix_asvs <- c(
  "3bf02d83576f1a60930d9a073c4a7d7e",  # 100% identity, 100% coverage
  "8f63bbc73058beaaf9c53b38024c529b"   # 98.819% identity, 100% coverage
)

seed <- 123
n_permutations <- 999
deseq_alpha <- 0.05
min_total_count_genus <- 10

plant_compartments <- c("Rhizosphere", "Endosphere", "Phyllosphere")
health_levels <- c(
  "Bulk soil",
  "Healthy Rhizosphere", "Diseased Rhizosphere",
  "Healthy Endosphere", "Diseased Endosphere",
  "Healthy Phyllosphere", "Diseased Phyllosphere"
)

alpha_metrics <- c(
  "Observed features" = "observed_features\\.csv$",
  "Shannon's Diversity" = "shannon\\.csv$"
)

suppressPackageStartupMessages({
  library(phyloseq)
  library(tidyverse)
  library(vegan)
  library(ape)
  library(DESeq2)
})

# Check required inputs before starting computational analyses.
required_inputs <- c(physeq_path, alpha_qzv_file, metadata_file)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs)) {
  stop("Missing input files in working directory: ",
       paste(missing_inputs, collapse = ", "),
       "\nSee REPRODUCIBILITY_README.md for instructions.")
}

dir.create(outdir, showWarnings = FALSE)

# ----------------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------------

otu_taxa_rows <- function(ps) {
  m <- as(otu_table(ps), "matrix")
  if (!taxa_are_rows(ps)) m <- t(m)
  m
}

otu_samples_rows <- function(ps) t(otu_taxa_rows(ps))

subset_ps <- function(ps, sample_ids) {
  ps <- prune_samples(sample_ids, ps)
  ps <- prune_samples(sample_sums(ps) > 0, ps)
  prune_taxa(taxa_sums(ps) > 0, ps)
}

# Missing or empty taxonomic labels are treated as unassigned
assigned <- function(x) {
  x <- as.character(x)
  x[is.na(x) | x == ""] <- NA_character_
  x
}

health_binary <- function(x) {
  case_when(
    grepl("Healthy", as.character(x)) ~ "Healthy",
    grepl("Diseased", as.character(x)) ~ "Diseased",
    TRUE ~ NA_character_
  )
}

# PCoA variance, PERMANOVA, ANOSIM (999 permutations) and betadisper for one distance matrix
beta_stats <- function(dist_obj, grouping) {
  df <- data.frame(grouping = grouping)
  set.seed(seed)
  per <- adonis2(dist_obj ~ grouping, data = df, permutations = n_permutations)
  set.seed(seed)
  ano <- anosim(dist_obj, grouping, permutations = n_permutations)
  bd <- anova(betadisper(dist_obj, grouping))
  eig <- ape::pcoa(as.matrix(dist_obj))$values$Relative_eig
  tibble(
    PCoA_axis1_var = 100 * eig[1], PCoA_axis2_var = 100 * eig[2],
    PERMANOVA_R2 = per$R2[1], PERMANOVA_F = per$F[1], PERMANOVA_p = per$`Pr(>F)`[1],
    ANOSIM_R = as.numeric(ano$statistic), ANOSIM_p = ano$signif,
    betadisper_F = bd$`F value`[1], betadisper_p = bd$`Pr(>F)`[1]
  )
}

beta_rows <- function(ps_sub, analysis, compartment) {
  mat <- otu_samples_rows(ps_sub)
  grouping <- as.character(sample_data(ps_sub)$Health_status)
  dists <- list(
    "Bray-Curtis" = vegdist(mat, method = "bray"),
    "Jaccard" = vegdist((mat > 0) * 1, method = "jaccard", binary = TRUE)
  )
  imap_dfr(dists, function(d, metric) {
    tibble(Analysis = analysis, Month = "October", Compartment = compartment, Metric = metric,
           n_samples = nsamples(ps_sub)) %>%
      bind_cols(beta_stats(d, grouping))
  })
}

# ----------------------------------------------------------------------------
# Input data
# ----------------------------------------------------------------------------

ps_raw <- readRDS(physeq_path)
ps_raw <- prune_samples(sample_sums(ps_raw) > 0, ps_raw)
ps0 <- prune_taxa(taxa_sums(ps_raw) > 0, ps_raw)
message("ITS input: ", nsamples(ps0), " samples; ", ntaxa(ps0), " ASVs")
# Expected manuscript totals: 73 samples, 6759 ASVs. Differences require review.
if (nsamples(ps0) != 73L || ntaxa(ps0) != 6759L) {
  warning("ITS RDS counts do not match the manuscript totals (73 samples, 6759 ASVs).")
}

sample_df <- as(sample_data(ps0), "data.frame") %>%
  rownames_to_column("SampleID") %>%
  mutate(across(c(Month, Compartment, Health_status), as.character),
         Reads = sample_sums(ps0)[SampleID])

otu_mat <- otu_taxa_rows(ps0)
tax_df <- as.data.frame(tax_table(ps0)) %>% rownames_to_column("ASV_ID") %>% as_tibble()

# ----------------------------------------------------------------------------
# Table S2 - sequencing output and taxonomic recovery
# ----------------------------------------------------------------------------

detected <- function(ids) rownames(otu_mat)[rowSums(otu_mat[, ids, drop = FALSE] > 0) > 0]
n_genera <- function(asv_ids) n_distinct(na.omit(assigned(tax_df$Genus[tax_df$ASV_ID %in% asv_ids])))

table_s2 <- sample_df %>%
  group_by(Compartment, Health_status) %>%
  summarise(
    n_samples = n(),
    total_reads = sum(Reads), mean_reads = mean(Reads), median_reads = median(Reads),
    min_reads = min(Reads), max_reads = max(Reads), sd_reads = sd(Reads),
    n_ASVs_detected = length(detected(SampleID)),
    n_genera_detected = n_genera(detected(SampleID)),
    .groups = "drop"
  )

global_totals <- tibble(
  n_samples = nsamples(ps0),
  total_reads = sum(sample_sums(ps0)),
  total_ASVs = ntaxa(ps0),
  total_genera = n_genera(tax_df$ASV_ID)
)

write_csv(table_s2, file.path(outdir, "Table_S2_ITS_sequencing_summary.csv"))
write_csv(global_totals, file.path(outdir, "Table_S2_ITS_global_totals.csv"))

# ----------------------------------------------------------------------------
# Table S3 - alpha diversity (QIIME 2 rarefaction output, depth 39,000)
# ----------------------------------------------------------------------------

metadata_raw <- read.delim(metadata_file, check.names = FALSE, stringsAsFactors = FALSE)
id_col <- intersect(c("sample-id", "SampleID", "sample_id", "#SampleID"), colnames(metadata_raw))[1]
metadata_alpha <- metadata_raw %>%
  transmute(sample_id = as.character(.data[[id_col]]), Month = as.character(Month),
            Compartment = as.character(Compartment), Health_status = as.character(Health_status))

qzv_dir <- file.path(tempdir(), "alpha_ITS")
unlink(qzv_dir, recursive = TRUE)
unzip(alpha_qzv_file, exdir = qzv_dir)
qzv_files <- list.files(qzv_dir, recursive = TRUE, full.names = TRUE)

read_alpha <- function(pattern, metric) {
  f <- grep(pattern, qzv_files, value = TRUE)
  f <- if (any(grepl("[/\\\\]data[/\\\\]", f))) grep("[/\\\\]data[/\\\\]", f, value = TRUE)[1] else f[1]
  wide <- read_csv(f, show_col_types = FALSE)
  names(wide) <- gsub("-", "_", tolower(gsub("^#+", "", names(wide))))

  long <- wide %>%
    select(sample_id, matches("^depth_[0-9]+_iter_[0-9]+$")) %>%
    mutate(sample_id = as.character(sample_id)) %>%
    pivot_longer(-sample_id, names_to = c("depth", "iter"),
                 names_pattern = "^depth_([0-9]+)_iter_([0-9]+)$", values_to = "value") %>%
    mutate(depth = as.numeric(depth))
  if (!alpha_depth %in% long$depth) stop("Depth ", alpha_depth, " not found in ", alpha_qzv_file)

  long %>%
    filter(depth == alpha_depth) %>%
    group_by(sample_id) %>%
    summarise(alpha_value = mean(as.numeric(value), na.rm = TRUE), .groups = "drop") %>%
    inner_join(metadata_alpha, by = "sample_id") %>%
    filter(tolower(Month) == "october") %>%
    mutate(Metric = metric, Health_status = factor(Health_status, levels = health_levels))
}

alpha_long <- imap_dfr(alpha_metrics, read_alpha)

alpha_summary <- alpha_long %>%
  group_by(Metric, Compartment, Health_status) %>%
  summarise(
    n = n(), mean = mean(alpha_value), sd = sd(alpha_value), se = sd / sqrt(n),
    median = median(alpha_value), q1 = quantile(alpha_value, 0.25), q3 = quantile(alpha_value, 0.75),
    min = min(alpha_value), max = max(alpha_value),
    .groups = "drop"
  )

# One-way ANOVA across all health-status groups followed by Tukey's HSD;
# the healthy-vs-diseased contrast of each compartment is reported, with a Wilcoxon test
target_contrasts <- paste0("Diseased ", plant_compartments, "-Healthy ", plant_compartments)

alpha_tukey <- alpha_long %>%
  group_by(Metric) %>%
  group_modify(~ as_tibble(TukeyHSD(aov(alpha_value ~ Health_status, data = .x))$Health_status,
                           rownames = "contrast")) %>%
  ungroup() %>%
  filter(contrast %in% target_contrasts) %>%
  transmute(Metric, Compartment = str_extract(contrast, paste(plant_compartments, collapse = "|")),
            contrast, estimate = diff, conf.low = lwr, conf.high = upr, tukey_p_adj = `p adj`)

alpha_wilcox <- alpha_long %>%
  filter(Compartment %in% plant_compartments) %>%
  group_by(Metric, Compartment) %>%
  group_modify(~ tibble(wilcox_p = wilcox.test(alpha_value ~ droplevels(Health_status),
                                               data = .x, exact = FALSE)$p.value)) %>%
  ungroup()

alpha_tests <- alpha_tukey %>%
  left_join(alpha_wilcox, by = c("Metric", "Compartment")) %>%
  mutate(result = case_when(
    tukey_p_adj < 0.05 & estimate > 0 ~ "Diseased > Healthy",
    tukey_p_adj < 0.05 & estimate < 0 ~ "Healthy > Diseased",
    TRUE ~ "No significant difference"
  ))

table_s3 <- alpha_summary %>%
  left_join(alpha_tests, by = c("Metric", "Compartment")) %>%
  mutate(Month = "October", .after = Metric)
write_csv(table_s3, file.path(outdir, "Table_S3_ITS_alpha_diversity.csv"))

# ----------------------------------------------------------------------------
# Table S4 - beta diversity
# ----------------------------------------------------------------------------

beta_global <- beta_rows(ps0, "Global", NA_character_)

beta_compartment <- map_dfr(plant_compartments, function(comp) {
  ids <- sample_df$SampleID[sample_df$Compartment == comp & grepl("^(Healthy|Diseased)", sample_df$Health_status)]
  beta_rows(subset_ps(ps0, ids), "Within-compartment", comp)
})

table_s4 <- bind_rows(beta_global, beta_compartment)
write_csv(table_s4, file.path(outdir, "Table_S4_ITS_beta_diversity.csv"))

# ----------------------------------------------------------------------------
# Table S6 - taxonomic composition (relative abundance per sample)
# ----------------------------------------------------------------------------

rel_long <- sweep(otu_mat, 2, pmax(colSums(otu_mat), 1), "/") %>%
  as.data.frame() %>%
  rownames_to_column("ASV_ID") %>%
  pivot_longer(-ASV_ID, names_to = "SampleID", values_to = "Relative_abundance") %>%
  left_join(tax_df, by = "ASV_ID") %>%
  left_join(sample_df %>% select(SampleID, Month, Compartment, Health_status), by = "SampleID")

composition <- function(rank, n_top) {
  summ <- rel_long %>%
    mutate(Taxon = replace_na(assigned(.data[[rank]]), paste0("Unclassified_", rank))) %>%
    group_by(SampleID, Month, Compartment, Health_status, Taxon) %>%
    summarise(Relative_abundance = sum(Relative_abundance), .groups = "drop") %>%
    group_by(Month, Compartment, Health_status, Taxon) %>%
    summarise(
      mean_rel_abundance_pct = 100 * mean(Relative_abundance),
      prevalence_pct = 100 * mean(Relative_abundance > 0),
      .groups = "drop"
    )
  top <- summ %>% group_by(Taxon) %>% summarise(m = mean(mean_rel_abundance_pct)) %>%
    slice_max(m, n = n_top, with_ties = FALSE) %>% pull(Taxon)
  summ %>% filter(Taxon %in% top) %>% mutate(Taxonomic_level = rank, .before = 1)
}

table_s6 <- bind_rows(
  composition("Phylum", 15), composition("Family", 20),
  composition("Genus", 20), composition("Species", 20)
) %>%
  mutate(Health_status = factor(Health_status, levels = health_levels)) %>%
  arrange(Taxonomic_level, Health_status, desc(mean_rel_abundance_pct))
write_csv(table_s6, file.path(outdir, "Table_S6_ITS_taxonomic_composition.csv"))

# ----------------------------------------------------------------------------
# Table S7 - genus-level differential abundance (DESeq2)
# ----------------------------------------------------------------------------

genus <- assigned(tax_df$Genus[match(rownames(otu_mat), tax_df$ASV_ID)])
genus_counts_all <- rowsum(otu_mat[!is.na(genus), , drop = FALSE], group = genus[!is.na(genus)], reorder = FALSE)

run_deseq_genus <- function(comp) {
  coldata <- sample_df %>%
    filter(Compartment == comp, grepl("^(Healthy|Diseased)", Health_status)) %>%
    transmute(SampleID, Condition = factor(health_binary(Health_status), levels = c("Healthy", "Diseased"))) %>%
    column_to_rownames("SampleID")

  counts <- genus_counts_all[, rownames(coldata), drop = FALSE]
  counts <- counts[rowSums(counts) >= min_total_count_genus, , drop = FALSE]

  dds <- DESeqDataSetFromMatrix(round(counts), colData = coldata, design = ~ Condition)
  dds <- DESeq(dds[rowSums(counts(dds)) > 0, ], fitType = "parametric")
  results(dds, contrast = c("Condition", "Diseased", "Healthy"), alpha = deseq_alpha) %>%
    as.data.frame() %>%
    rownames_to_column("Genus") %>%
    mutate(Month = "October", Compartment = comp, n_tested = n(), .before = 1)
}

deseq_all <- map_dfr(plant_compartments, run_deseq_genus)

genus_lookup <- tax_df %>%
  mutate(Genus = assigned(Genus)) %>%
  filter(!is.na(Genus)) %>%
  group_by(Genus) %>%
  summarise(across(c(Phylum, Family), ~ first(na.omit(assigned(.x)))), .groups = "drop")

table_s7 <- deseq_all %>%
  filter(!is.na(padj), padj < deseq_alpha) %>%
  left_join(genus_lookup, by = "Genus") %>%
  mutate(Direction = ifelse(log2FoldChange > 0, "Diseased > Healthy", "Healthy > Diseased")) %>%
  select(Month, Compartment, n_tested, Genus, Family, Phylum, baseMean, log2FoldChange, lfcSE, padj, Direction) %>%
  arrange(Compartment, padj)

write_csv(deseq_all, file.path(outdir, "DESeq2_ITS_genus_all_tested.csv"))
write_csv(table_s7, file.path(outdir, "Table_S7_ITS_DESeq2_genus.csv"))

# ----------------------------------------------------------------------------
# Tables S8 and S9 - pathogen-associated signals
# ----------------------------------------------------------------------------

# Relative abundance (% of total ITS reads per sample) of a set of ASVs
signal_samples <- function(asv_ids, compartments, signal) {
  counts <- otu_samples_rows(ps0)
  tibble(
    SampleID = rownames(counts),
    Target_reads = rowSums(counts[, intersect(asv_ids, colnames(counts)), drop = FALSE]),
    Total_reads = rowSums(counts)
  ) %>%
    left_join(sample_df %>% select(SampleID, Compartment, Health_status), by = "SampleID") %>%
    filter(Compartment %in% compartments, grepl("^(Healthy|Diseased)", Health_status)) %>%
    mutate(Signal = signal, Relative_abundance_pct = 100 * Target_reads / Total_reads)
}

# Group summaries and Mann-Whitney U tests (healthy vs diseased within each compartment)
signal_table <- function(samples) {
  samples %>%
    group_by(Signal, Compartment) %>%
    group_modify(~ {
      mw <- wilcox.test(Relative_abundance_pct ~ Health_status, data = .x, exact = FALSE)$p.value
      .x %>%
        group_by(Health_status) %>%
        summarise(
          n = n(),
          mean_rel_abundance_pct = mean(Relative_abundance_pct),
          se_rel_abundance_pct = sd(Relative_abundance_pct) / sqrt(n()),
          prevalence_pct = 100 * mean(Relative_abundance_pct > 0),
          .groups = "drop"
        ) %>%
        mutate(MannWhitney_p = mw)
    }) %>%
    ungroup()
}

# Table S8: BLAST-curated R. necatrix-compatible ASVs, root-associated compartments,
# with Benjamini-Hochberg correction across the root-associated tests
if (!all(rnecatrix_asvs %in% taxa_names(ps0))) stop("R. necatrix-compatible ASV IDs not found in ", physeq_path)
table_s8 <- signal_samples(rnecatrix_asvs, c("Rhizosphere", "Endosphere"), "R. necatrix-compatible ASVs") %>%
  signal_table()
bh <- table_s8 %>% distinct(Compartment, MannWhitney_p) %>%
  mutate(MannWhitney_p_BH = p.adjust(MannWhitney_p, method = "BH"))
table_s8 <- table_s8 %>% left_join(select(bh, Compartment, MannWhitney_p_BH), by = "Compartment")
write_csv(table_s8, file.path(outdir, "Table_S8_ITS_Rnecatrix_compatible_signals.csv"))

# Table S9: Botryosphaeriaceae (family) and Neofusicoccum (genus) in the phyllosphere
botryo_asvs <- tax_df$ASV_ID[assigned(tax_df$Family) %in% "Botryosphaeriaceae"]
neofusicoccum_asvs <- tax_df$ASV_ID[assigned(tax_df$Genus) %in% "Neofusicoccum"]
table_s9 <- bind_rows(
  signal_samples(botryo_asvs, "Phyllosphere", "Botryosphaeriaceae"),
  signal_samples(neofusicoccum_asvs, "Phyllosphere", "Neofusicoccum")
) %>%
  signal_table()
write_csv(table_s9, file.path(outdir, "Table_S9_ITS_branch_dieback_signals.csv"))

# ----------------------------------------------------------------------------
# Table S5 - tree-level robustness analysis
# ----------------------------------------------------------------------------
# Spatial subsamples are summed into one profile per Compartment x Health status x Tree.
# `Trees` holds subsample codes (e.g. A1RO-A4RO); the biological tree is the leading letter.

md_tree_in <- sample_df %>%
  filter(Compartment %in% plant_compartments, grepl("^(Healthy|Diseased)", Health_status)) %>%
  mutate(
    Health_binary = health_binary(Health_status),
    Tree = sub("^([A-Z]).*$", "\\1", Trees),
    Tree_level_ID = paste(Compartment, Health_binary, Tree, sep = "__")
  )

otu_tree <- rowsum(otu_samples_rows(ps0)[md_tree_in$SampleID, , drop = FALSE],
                   group = md_tree_in$Tree_level_ID, reorder = FALSE)
md_tree <- md_tree_in %>%
  distinct(Tree_level_ID, Compartment, Health_binary) %>%
  column_to_rownames("Tree_level_ID")
md_tree <- md_tree[rownames(otu_tree), , drop = FALSE]

table_s5 <- map_dfr(plant_compartments, function(comp) {
  ids <- rownames(md_tree)[md_tree$Compartment == comp]
  otu_c <- otu_tree[ids, colSums(otu_tree[ids, , drop = FALSE]) > 0, drop = FALSE]
  md_c <- md_tree[ids, , drop = FALSE]
  d <- vegdist(otu_c, method = "bray")

  set.seed(seed)
  per <- adonis2(d ~ Health_binary, data = md_c, permutations = n_permutations)
  bd <- permutest(betadisper(d, md_c$Health_binary), permutations = n_permutations)

  tibble(Analysis = "October tree-level", Month = "October", Compartment = comp, Term = "Health status",
         n_profiles = nrow(md_c), PERMANOVA_R2 = per$R2[1], PERMANOVA_F = per$F[1],
         PERMANOVA_p = per$`Pr(>F)`[1], betadisper_F = bd$tab$F[1], betadisper_p = bd$tab$`Pr(>F)`[1])
})

write_csv(table_s5, file.path(outdir, "Table_S5_ITS_tree_level_robustness.csv"))

# ----------------------------------------------------------------------------
# Session information
# ----------------------------------------------------------------------------

writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo_ITS.txt"))
message("ITS analysis completed. Tables written to ", outdir)
