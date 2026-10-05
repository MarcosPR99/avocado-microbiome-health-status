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

## Data availability

Raw 16S rRNA gene and ITS amplicon sequencing reads are not included in this repository.

Raw reads will be deposited in a public sequence archive and cited in the manuscript using the corresponding BioProject accession number.

**Raw sequencing reads:** accession pending.

## Code availability

This repository contains the scripts used for QIIME2-based processing and downstream analyses in R.

A permanent archived version of this repository will be generated through Zenodo after the final version of the code and tables is available.

**Archived code DOI:** pending.

## Software requirements

The analyses were performed using QIIME2 and R.

Main R packages used across the analysis scripts include:

- phyloseq
- vegan
- tidyverse
- dplyr
- tidyr
- readr
- ggplot2
- patchwork
- ape
- DESeq2
- multcompView

Some scripts were originally run in local or high-performance computing environments. File paths may need to be adapted before re-running the analyses on a different system.

## Large files

Large sequencing and analysis files are not tracked in this repository, including:

- FASTQ files
- QIIME2 `.qza` and `.qzv` files
- R serialized objects such as `.rds` or `.RDS`
- large figure files

## License

The analysis code in this repository is released under the MIT License.
