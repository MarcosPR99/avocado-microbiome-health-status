# Avocado-associated microbiome responses to fungal disease

This repository contains analysis scripts and cleaned metadata associated with the manuscript:

**Fungal diseases reshape avocado-associated multi-kingdom microbiomes through compartment-specific restructuring and network rewiring**

The study analyses bacterial 16S rRNA gene and fungal ITS amplicon datasets from avocado-associated compartments, including bulk soil, rhizosphere, root endosphere and phyllosphere, comparing healthy and diseased trees.

## Repository structure

code/
- 01a_16S_dada2_denoising.sh
- 01b_16S_alignment_mafft.sh
- 01c_16S_alignment_mask.sh
- 01d_16S_fasttree_rooted.sh
- 01e_16S_taxonomy_classification.sh
- 02_ITS_import_to_phyloseq.R
- 03_16S_downstream_analysis.R
- 04_ITS_downstream_analysis.R

metadata/
- metadata_16S_clean.tsv
- metadata_ITS_clean.tsv

processed_tables/
- Processed output tables used for downstream analyses.

supplementary_tables/
- Supplementary tables associated with the manuscript.

docs/
- Additional documentation.

## Data availability

Raw sequencing reads are not included in this repository. They will be deposited in a public sequence archive under the BioProject accession provided in the manuscript.

This repository contains cleaned metadata, processing scripts and downstream analysis code needed to reproduce the analyses reported in the manuscript.

## Notes

Some scripts were originally run in local or high-performance computing environments. File paths may need to be adapted before re-running the analyses on a different system.
