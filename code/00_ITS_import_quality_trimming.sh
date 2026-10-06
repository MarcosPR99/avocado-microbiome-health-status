#!/bin/bash

# QIIME2 ITS preprocessing workflow
# This script follows the same preprocessing structure used for the 16S workflow,
# adapted to the ITS dataset and file naming used in the analysis.

# Load QIIME2 environment on Picasso
module load qiime2/2023.5.1

############################################################
# 1. Import paired-end FASTQ reads
############################################################

# The ITS metadata/manifest file should include:
# - sample-id
# - forward-absolute-filepath
# - reverse-absolute-filepath
# - experimental metadata

qiime tools import \
  --input-path metadata_octubre_its_samples.tsv \
  --output-path articulo_ITS_raw.qza \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-format PairedEndFastqManifestPhred33V2

############################################################
# 2. Initial read quality summary
############################################################

qiime demux summarize \
  --p-n 50000 \
  --i-data articulo_ITS_raw.qza \
  --o-visualization articulo_ITS_raw.qzv

############################################################
# 3. Primer trimming with cutadapt
############################################################

# ITS primers/adapters were trimmed before DADA2.
# Replace the primer sequences below only if re-running the workflow from raw reads.
# The downstream analysis used the trimmed ITS object that produced:
# denoising/table_dada2_ITS_front.qza
# denoising/rep-seqs_dada2_ITS_front.qza
# taxonomy/taxonomy_dada2_ITS_front.qza

qiime cutadapt trim-paired \
  --p-cores 20 \
  --i-demultiplexed-sequences articulo_ITS_raw.qza \
  --p-front-f ITS_FORWARD_PRIMER_SEQUENCE \
  --p-front-r ITS_REVERSE_PRIMER_SEQUENCE \
  --o-trimmed-sequences articulo_ITS_trimmed.qza \
  --p-match-read-wildcards \
  --p-match-adapter-wildcards \
  --verbose

############################################################
# 4. Post-trimming quality summary
############################################################

qiime demux summarize \
  --p-n 50000 \
  --i-data articulo_ITS_trimmed.qza \
  --o-visualization articulo_ITS_trimmed.qzv
