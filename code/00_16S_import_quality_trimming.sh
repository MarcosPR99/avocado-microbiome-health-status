#!/bin/bash

# QIIME2 16S preprocessing workflow
# Working directory used in the original analysis:
# ~/ja_project/qiime2_16s_artículo_final

# Load QIIME2 environment on Picasso
module load qiime2/2023.5.1

############################################################
# 1. Import paired-end FASTQ reads
############################################################

# metadata_samples.tsv was used as a manifest and includes:
# - sample-id
# - forward-absolute-filepath
# - reverse-absolute-filepath
# - experimental metadata

qiime tools import \
  --input-path metadata_samples.tsv \
  --output-path articulo_16S_raw.qza \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-format PairedEndFastqManifestPhred33V2

############################################################
# 2. Initial read quality summary
############################################################

qiime demux summarize \
  --p-n 50000 \
  --i-data articulo_16S_raw.qza \
  --o-visualization articulo_16S_raw.qzv

############################################################
# 3. Primer trimming with cutadapt
############################################################

# 16S primers used in this analysis:
# Forward: AACMGGATTAGATACCCKG
# Reverse: ACGTCATCCCCACCTTCC

qiime cutadapt trim-paired \
  --p-cores 20 \
  --i-demultiplexed-sequences articulo_16S_raw.qza \
  --p-front-f AACMGGATTAGATACCCKG \
  --p-front-r ACGTCATCCCCACCTTCC \
  --o-trimmed-sequences articulo_16S_trimmed.qza \
  --p-match-read-wildcards \
  --p-match-adapter-wildcards \
  --verbose

############################################################
# 4. Post-trimming quality summary
############################################################

qiime demux summarize \
  --p-n 50000 \
  --i-data articulo_16S_trimmed.qza \
  --o-visualization articulo_16S_trimmed.qzv
