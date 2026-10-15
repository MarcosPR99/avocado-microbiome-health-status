# Avocado-associated microbiome responses to fungal disease

This repository contains processed amplicon sequencing data, sample metadata, and analysis scripts associated with the manuscript:

**Fungal diseases reshape avocado-associated multi-kingdom microbiomes through compartment-specific restructuring and network rewiring**

The study examines bacterial 16S rRNA gene and fungal ITS communities in avocado-associated compartments (bulk soil, rhizosphere, root endosphere, and phyllosphere), comparing healthy and symptomatic trees.

## Repository structure

```text
data/
  16S_ASV_table.tsv.gz
  16S_taxonomy.tsv
  ITS_ASV_table.tsv
  ITS_taxonomy.tsv

metadata/
  metadata_16S_clean.tsv
  metadata_ITS_clean.tsv

code/
  00_16S_import_quality_trimming.sh
  00_ITS_import_quality_trimming.sh
  01a_16S_dada2_denoising.sh
  01b_16S_alignment_mafft.sh
  01c_16S_alignment_mask.sh
  01d_16S_fasttree_rooted.sh
  01e_16S_taxonomy_classification.sh
  02_ITS_import_to_phyloseq.R
  03_16S_downstream_analysis.R
  04_ITS_downstream_analysis.R
  05_network_analysis.py
  README.md
```

## Processed data

- `data/16S_ASV_table.tsv.gz`: bacterial ASV count table (gzip-compressed TSV).
- `data/16S_taxonomy.tsv`: bacterial ASV taxonomic assignments.
- `data/ITS_ASV_table.tsv`: fungal ASV count table.
- `data/ITS_taxonomy.tsv`: fungal ASV taxonomic assignments.
- `metadata/metadata_16S_clean.tsv` and `metadata/metadata_ITS_clean.tsv`: sample metadata.

## Analysis code

The scripts in `code/` cover QIIME 2-based read processing, bacterial and fungal downstream analyses in R, and cross-kingdom SparCC/NetworkX network analyses in Python. See `code/README.md` for workflow details.

Analyses use QIIME 2, R (including phyloseq, vegan and DESeq2), and Python (including NetworkX). Versions and analytical parameters are documented in the manuscript and associated scripts.

Some scripts rely on intermediate QIIME 2 archives or phyloseq RDS objects not stored here and may require path adjustments or regeneration of intermediate inputs before execution. The scripts have not all been independently rerun from raw reads in this repository.

## Raw sequencing data

Raw 16S rRNA gene and ITS amplicon reads have been deposited in the NCBI Sequence Read Archive under BioProject accession **PRJNA1540021**. Raw FASTQ files are not duplicated in this repository.

## License

Analysis code is released under the MIT License.
