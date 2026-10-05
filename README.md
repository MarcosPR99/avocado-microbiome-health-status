# Avocado-associated microbiome responses to fungal disease

This repository contains analysis scripts, cleaned metadata and processed output folders associated with the manuscript:

**Fungal diseases reshape avocado-associated multi-kingdom microbiomes through compartment-specific restructuring and network rewiring**

The study analyses bacterial 16S rRNA gene and fungal ITS amplicon datasets from avocado-associated compartments, including bulk soil, rhizosphere, root endosphere and phyllosphere, comparing healthy and diseased avocado trees.

## Repository structure

```text
code/
  01a_16S_dada2_denoising.sh
  01b_16S_alignment_mafft.sh
  01c_16S_alignment_mask.sh
  01d_16S_fasttree_rooted.sh
  01e_16S_taxonomy_classification.sh
  02_ITS_import_to_phyloseq.R
  03_16S_downstream_analysis.R
  04_ITS_downstream_analysis.R

metadata/
  metadata_16S_clean.tsv
  metadata_ITS_clean.tsv

processed_tables/
  Processed output tables used for downstream analyses.

supplementary_tables/
  Supplementary tables associated with the manuscript.

docs/
  Additional documentation.
