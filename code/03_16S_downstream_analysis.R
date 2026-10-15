#!/usr/bin/env Rscript

###############################################################################
# 03_16S_downstream_analysis.R
#
# Downstream analyses of 16S rRNA gene amplicon data (July and October) for:
# "Fungal diseases reshape avocado-associated multi-kingdom microbiomes through
#  compartment-specific restructuring and network rewiring"
#
# Outputs (16S part of each supplementary table):
#   Table S2  - sequencing output and taxonomic recovery
#   Table S3  - alpha diversity
#   Table S4  - beta diversity (Bray-Curtis, Jaccard, weighted/unweighted UniFrac)
#   Table S6  - taxonomic composition (phylum and genus)
#   Table S7  - differentially abundant genera (DESeq2)
#   Table S5 - tree-level robustness analysis
#
# Inputs (working directory; see REPRODUCIBILITY_README.md):
#   Phyloseq_16S.rds
#   alpha-rarefaction_withPhy_depth_18288.qzv
###############################################################################

# ----------------------------------------------------------------------------
# Configuration
# ----------------------------------------------------------------------------

physeq_path <- "Phyloseq_16S.rds"
alpha_qzv_file <- "alpha-rarefaction_withPhy_depth_18288.qzv"
alpha_depth <- 18288
outdir <- "results_16S"

set.seed(123)
n_permutations <- 999
deseq_alpha <- 0.05
min_total_count_genus <- 10

month_levels <- c("July", "October")
compartment_levels <- c("Bulk soil", "Rhizosphere", "Endosphere", "Phyllosphere")
plant_compartments <- c("Rhizosphere", "Endosphere", "Phyllosphere")
health_levels <- c(
  "Bulk soil",
  "Healthy Rhizosphere", "Diseased Rhizosphere",
  "Healthy Endosphere", "Diseased Endosphere",
  "Healthy Phyllosphere", "Diseased Phyllosphere"
)

alpha_metrics <- c(
  "Observed features" = "observed_features\\.csv$",
  "Chao1" = "chao1\\.csv$",
  "Faith's PD" = "faith_pd\\.csv$",
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
required_inputs <- c(physeq_path, alpha_qzv_file)
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

# Missing or empty genus labels are treated as unassigned
assigned_genus <- function(x) {
  x <- as.character(x)
  x[is.na(x) | x == ""] <- NA_character_
  x
}

health_binary <- function(x) {
  case_when(
    grepl("^Healthy", as.character(x)) ~ "Healthy",
    grepl("^Diseased", as.character(x)) ~ "Diseased",
    TRUE ~ as.character(x)
  )
}

n_genera <- function(asv_ids) {
  n_distinct(na.omit(assigned_genus(tax_df$Genus[tax_df$ASV_ID %in% asv_ids])))
}

# PCoA variance, PERMANOVA, ANOSIM and betadisper for one distance matrix
beta_stats <- function(dist_obj, grouping) {
  df <- data.frame(grouping = grouping)
  per <- adonis2(dist_obj ~ grouping, data = df, permutations = n_permutations)
  ano <- anosim(dist_obj, grouping, permutations = n_permutations)
  bd <- permutest(betadisper(dist_obj, grouping), permutations = n_permutations)
  eig <- ape::pcoa(dist_obj)$values$Relative_eig
  tibble(
    PCoA_axis1_var = 100 * eig[1], PCoA_axis2_var = 100 * eig[2],
    PERMANOVA_R2 = per$R2[1], PERMANOVA_F = per$F[1], PERMANOVA_p = per$`Pr(>F)`[1],
    ANOSIM_R = as.numeric(ano$statistic), ANOSIM_p = ano$signif,
    betadisper_F = bd$tab$F[1], betadisper_p = bd$tab$`Pr(>F)`[1]
  )
}

beta_rows <- function(ps_sub, grouping, analysis, month, compartment, unifrac = FALSE) {
  if (unifrac) {
    ps_sub <- prune_taxa(intersect(taxa_names(ps_sub), phy_tree(ps_sub)$tip.label), ps_sub)
    dists <- list(
      "Weighted UniFrac" = phyloseq::distance(ps_sub, method = "unifrac", weighted = TRUE, fast = TRUE),
      "Unweighted UniFrac" = phyloseq::distance(ps_sub, method = "unifrac", weighted = FALSE, fast = TRUE)
    )
  } else {
    mat <- otu_samples_rows(ps_sub)
    dists <- list(
      "Bray-Curtis" = vegdist(mat, method = "bray"),
      "Jaccard" = vegdist((mat > 0) * 1, method = "jaccard", binary = TRUE)
    )
  }
  imap_dfr(dists, function(d, metric) {
    tibble(Analysis = analysis, Month = month, Compartment = compartment, Metric = metric,
           n_samples = nsamples(ps_sub)) %>%
      bind_cols(beta_stats(d, grouping))
  })
}

# ----------------------------------------------------------------------------
# Input data
# ----------------------------------------------------------------------------

ps0 <- readRDS(physeq_path)
ps0 <- prune_samples(sample_sums(ps0) > 0, ps0)
ps0 <- prune_taxa(taxa_sums(ps0) > 0, ps0)
message("16S input: ", nsamples(ps0), " samples; ", ntaxa(ps0), " ASVs")
# Expected manuscript totals: 152 samples, 112212 ASVs. Differences require review.
if (nsamples(ps0) != 152L || ntaxa(ps0) != 112212L) {
  warning("16S RDS counts do not match the manuscript totals (152 samples, 112212 ASVs).")
}

sample_df <- as(sample_data(ps0), "data.frame") %>%
  rownames_to_column("SampleID") %>%
  mutate(
    Month = factor(as.character(Month), levels = month_levels),
    Compartment = factor(as.character(Compartment), levels = compartment_levels),
    Health_status = factor(as.character(Health_status), levels = health_levels),
    Reads = sample_sums(ps0)[SampleID]
  )
sample_data(ps0) <- sample_data(column_to_rownames(sample_df, "SampleID"))

otu_mat <- otu_taxa_rows(ps0)
tax_df <- as.data.frame(tax_table(ps0)) %>% rownames_to_column("ASV_ID") %>% as_tibble()

# ----------------------------------------------------------------------------
# Table S2 - sequencing output and taxonomic recovery
# ----------------------------------------------------------------------------

detected <- function(ids) rownames(otu_mat)[rowSums(otu_mat[, ids, drop = FALSE] > 0) > 0]

table_s2 <- sample_df %>%
  group_by(Month, Compartment, Health_status) %>%
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

write_csv(table_s2, file.path(outdir, "Table_S2_16S_sequencing_summary.csv"))
write_csv(global_totals, file.path(outdir, "Table_S2_16S_global_totals.csv"))

# ----------------------------------------------------------------------------
# Table S3 - alpha diversity (QIIME 2 rarefaction output, depth 18,288)
# ----------------------------------------------------------------------------

qzv_dir <- file.path(tempdir(), "alpha_16S")
unlink(qzv_dir, recursive = TRUE)
unzip(alpha_qzv_file, exdir = qzv_dir)
qzv_files <- list.files(qzv_dir, recursive = TRUE, full.names = TRUE)

read_alpha <- function(pattern, metric) {
  f <- grep(pattern, qzv_files, value = TRUE)
  f <- if (any(grepl("[/\\\\]data[/\\\\]", f))) grep("[/\\\\]data[/\\\\]", f, value = TRUE)[1] else f[1]
  wide <- read_csv(f, show_col_types = FALSE)
  names(wide) <- gsub("-", "_", tolower(gsub("^#+", "", names(wide))))
  if (!"sample_id" %in% names(wide)) names(wide)[names(wide) %in% c("sampleid", "id")][1] <- "sample_id"

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
    inner_join(sample_df %>% transmute(sample_id = SampleID, Month, Compartment, Health_status),
               by = "sample_id") %>%
    mutate(Metric = metric)
}

alpha_long <- imap_dfr(alpha_metrics, read_alpha)

alpha_summary <- alpha_long %>%
  group_by(Metric, Month, Compartment, Health_status) %>%
  summarise(
    n = n(), mean = mean(alpha_value), sd = sd(alpha_value), se = sd / sqrt(n),
    median = median(alpha_value), q1 = quantile(alpha_value, 0.25), q3 = quantile(alpha_value, 0.75),
    min = min(alpha_value), max = max(alpha_value),
    .groups = "drop"
  )

# Healthy vs diseased within each compartment and month: ANOVA + Tukey's HSD, and Wilcoxon
alpha_tests <- alpha_long %>%
  filter(Compartment != "Bulk soil") %>%
  group_by(Metric, Month, Compartment) %>%
  group_modify(~ {
    tk <- TukeyHSD(aov(alpha_value ~ Health_status, data = .x))$Health_status
    wt <- wilcox.test(alpha_value ~ droplevels(Health_status), data = .x, exact = FALSE)
    tibble(contrast = rownames(tk)[1], estimate = tk[1, "diff"],
           conf.low = tk[1, "lwr"], conf.high = tk[1, "upr"],
           tukey_p_adj = tk[1, "p adj"], wilcox_p = wt$p.value)
  }) %>%
  ungroup() %>%
  mutate(result = case_when(
    tukey_p_adj < 0.05 & estimate > 0 ~ "Diseased > Healthy",
    tukey_p_adj < 0.05 & estimate < 0 ~ "Healthy > Diseased",
    TRUE ~ "No significant difference"
  ))

table_s3 <- alpha_summary %>% left_join(alpha_tests, by = c("Metric", "Month", "Compartment"))
write_csv(table_s3, file.path(outdir, "Table_S3_16S_alpha_diversity.csv"))

# ----------------------------------------------------------------------------
# Table S4 - beta diversity
# ----------------------------------------------------------------------------

beta_global <- map_dfr(month_levels, function(m) {
  ps_m <- subset_ps(ps0, sample_df$SampleID[sample_df$Month == m])
  beta_rows(ps_m, sample_data(ps_m)$Health_status, "Global", m, NA_character_)
})

compartment_subset <- function(m, comp) {
  subset_ps(ps0, sample_df$SampleID[sample_df$Month == m & sample_df$Compartment == comp])
}
binary_grouping <- function(ps) factor(health_binary(sample_data(ps)$Health_status), levels = c("Healthy", "Diseased"))

beta_compartment <- map_dfr(month_levels, function(m) map_dfr(plant_compartments, function(comp) {
  ps_c <- compartment_subset(m, comp)
  beta_rows(ps_c, binary_grouping(ps_c), "Within-compartment", m, comp)
}))

beta_unifrac <- map_dfr(month_levels, function(m) map_dfr(plant_compartments, function(comp) {
  ps_c <- compartment_subset(m, comp)
  beta_rows(ps_c, binary_grouping(ps_c), "Within-compartment", m, comp, unifrac = TRUE)
}))

table_s4 <- bind_rows(beta_global, beta_compartment, beta_unifrac)
write_csv(table_s4, file.path(outdir, "Table_S4_16S_beta_diversity.csv"))

# ----------------------------------------------------------------------------
# Table S6 - taxonomic composition (relative abundance per sample)
# ----------------------------------------------------------------------------

ps_rel <- transform_sample_counts(ps0, function(x) if (sum(x) == 0) x else x / sum(x))

composition <- function(rank, n_top) {
  summ <- psmelt(tax_glom(ps_rel, taxrank = rank, NArm = FALSE)) %>%
    mutate(Taxon = as.character(.data[[rank]]),
           Taxon = ifelse(is.na(Taxon) | Taxon %in% c("", "NA"), paste0("Unclassified_", rank), Taxon)) %>%
    group_by(Month, Compartment, Health_status, Taxon) %>%
    summarise(
      mean_rel_abundance_pct = 100 * mean(Abundance),
      prevalence_pct = 100 * mean(Abundance > 0),
      .groups = "drop"
    )
  top <- summ %>% group_by(Taxon) %>% summarise(m = mean(mean_rel_abundance_pct)) %>%
    slice_max(m, n = n_top, with_ties = FALSE) %>% pull(Taxon)
  summ %>% filter(Taxon %in% top) %>% mutate(Taxonomic_level = rank, .before = 1)
}

table_s6 <- bind_rows(composition("Phylum", 15), composition("Genus", 20)) %>%
  arrange(Taxonomic_level, Month, Compartment, Health_status, desc(mean_rel_abundance_pct))
write_csv(table_s6, file.path(outdir, "Table_S6_16S_taxonomic_composition.csv"))

# ----------------------------------------------------------------------------
# Table S7 - genus-level differential abundance (DESeq2)
# ----------------------------------------------------------------------------

run_deseq_genus <- function(m, comp) {
  ps_c <- subset_ps(ps0, sample_df$SampleID[sample_df$Month == m & sample_df$Compartment == comp &
                                              grepl("^(Healthy|Diseased)", sample_df$Health_status)])
  counts <- otu_taxa_rows(ps_c)
  genus <- assigned_genus(tax_df$Genus[match(rownames(counts), tax_df$ASV_ID)])
  genus_counts <- rowsum(counts[!is.na(genus), , drop = FALSE], group = genus[!is.na(genus)], reorder = FALSE)
  genus_counts <- genus_counts[rowSums(genus_counts) >= min_total_count_genus, , drop = FALSE]

  coldata <- data.frame(
    Condition = factor(health_binary(sample_data(ps_c)$Health_status), levels = c("Healthy", "Diseased")),
    row.names = sample_names(ps_c)
  )[colnames(genus_counts), , drop = FALSE]

  dds <- DESeqDataSetFromMatrix(round(genus_counts), colData = coldata, design = ~ Condition)
  dds <- DESeq(dds[rowSums(counts(dds)) > 0, ], fitType = "parametric")
  results(dds, contrast = c("Condition", "Diseased", "Healthy"), alpha = deseq_alpha) %>%
    as.data.frame() %>%
    rownames_to_column("Genus") %>%
    mutate(Month = m, Compartment = comp, n_tested = n(), .before = 1)
}

deseq_all <- map_dfr(month_levels, function(m) map_dfr(plant_compartments, ~ run_deseq_genus(m, .x)))

phylum_lookup <- tax_df %>%
  mutate(Genus = assigned_genus(Genus)) %>%
  filter(!is.na(Genus)) %>%
  group_by(Genus) %>%
  summarise(Phylum = paste(unique(na.omit(Phylum)), collapse = " | "), .groups = "drop")

table_s7 <- deseq_all %>%
  filter(!is.na(padj), padj < deseq_alpha) %>%
  left_join(phylum_lookup, by = "Genus") %>%
  mutate(Direction = ifelse(log2FoldChange > 0, "Diseased > Healthy", "Healthy > Diseased")) %>%
  select(Month, Compartment, n_tested, Genus, Phylum, baseMean, log2FoldChange, lfcSE, padj, Direction) %>%
  arrange(Month, Compartment, padj)

write_csv(deseq_all, file.path(outdir, "DESeq2_16S_genus_all_tested.csv"))
write_csv(table_s7, file.path(outdir, "Table_S7_16S_DESeq2_genus.csv"))

# ----------------------------------------------------------------------------
# Table S5 - tree-level robustness analysis
# ----------------------------------------------------------------------------
# The four spatial subsamples of each tree are summed into one profile per
# Month x Compartment x Health status x Tree (column `Trees` in the metadata).

md_tree_in <- sample_df %>%
  filter(Compartment %in% plant_compartments, grepl("^(Healthy|Diseased)", Health_status)) %>%
  mutate(
    Health_binary = health_binary(Health_status),
    Tree_level_ID = paste(Month, Compartment, Health_binary, Trees, sep = "__")
  )

otu_tree <- rowsum(otu_samples_rows(ps0)[md_tree_in$SampleID, , drop = FALSE],
                   group = md_tree_in$Tree_level_ID, reorder = FALSE)
md_tree <- md_tree_in %>%
  distinct(Tree_level_ID, Month, Compartment, Health_binary) %>%
  column_to_rownames("Tree_level_ID")
md_tree <- md_tree[rownames(otu_tree), , drop = FALSE]
md_tree$Health_binary <- factor(md_tree$Health_binary, levels = c("Healthy", "Diseased"))

table_s5 <- map_dfr(plant_compartments, function(comp) {
  ids <- rownames(md_tree)[md_tree$Compartment == comp]
  otu_c <- otu_tree[ids, colSums(otu_tree[ids, , drop = FALSE]) > 0, drop = FALSE]
  md_c <- md_tree[ids, , drop = FALSE]

  # Combined July + October model: health status tested after sampling month
  combined <- adonis2(vegdist(otu_c, method = "bray") ~ Month + Health_binary,
                      data = md_c, permutations = n_permutations, by = "terms") %>%
    as.data.frame() %>%
    rownames_to_column("Term") %>%
    filter(Term %in% c("Month", "Health_binary")) %>%
    transmute(Analysis = "Combined July-October", Month = NA_character_, Compartment = comp, Term,
              n_profiles = nrow(md_c), PERMANOVA_R2 = R2, PERMANOVA_F = `F`, PERMANOVA_p = `Pr(>F)`)

  # Within-month models
  by_month <- map_dfr(month_levels, function(m) {
    md_m <- md_c[md_c$Month == m, , drop = FALSE]
    d_m <- vegdist(otu_c[rownames(md_m), , drop = FALSE], method = "bray")
    per <- adonis2(d_m ~ Health_binary, data = md_m, permutations = n_permutations)
    bd <- permutest(betadisper(d_m, md_m$Health_binary), permutations = n_permutations)
    tibble(Analysis = "Within month", Month = m, Compartment = comp, Term = "Health_binary",
           n_profiles = nrow(md_m), PERMANOVA_R2 = per$R2[1], PERMANOVA_F = per$F[1],
           PERMANOVA_p = per$`Pr(>F)`[1], betadisper_F = bd$tab$F[1], betadisper_p = bd$tab$`Pr(>F)`[1])
  })

  bind_rows(combined, by_month)
}) %>%
  mutate(Term = recode(Term, Health_binary = "Health status"))

write_csv(table_s5, file.path(outdir, "Table_S5_16S_tree_level_robustness.csv"))

# ----------------------------------------------------------------------------
# Session information
# ----------------------------------------------------------------------------

writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo_16S.txt"))
message("16S analysis completed. Tables written to ", outdir)
