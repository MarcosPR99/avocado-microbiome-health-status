#!/bin/bash
#SBATCH --mem=1800gb
#SBATCH --cpus-per-task=128
#SBATCH --error=dada2.%J.err
#SBATCH --output=dada2.%J.out



qiime dada2 denoise-paired \
  --i-demultiplexed-seqs articulo_16S_trimmed.qza \
  --p-trim-left-f 0 \
  --p-trim-left-r 0 \
  --p-trunc-len-f 222 \
  --p-trunc-len-r 223 \
  --p-min-overlap 10 \
  --p-max-ee-f 2 \
  --p-max-ee-r 2 \
  --p-pooling-method 'pseudo' \
  --p-chimera-method 'pooled' \
  --p-n-threads 128 \
  --o-table table_dada2_ee2_ee2_f222r223.qza \
  --o-representative-sequences rep-seqs_dada2_ee2_ee2_f222r223.qza \
  --o-denoising-stats denoising-stats_dada2_ee2_ee2_f222r223.qza \
  --verbose
