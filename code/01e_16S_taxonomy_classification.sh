#!/bin/bash
#SBATCH --mem=1800gb
#SBATCH --cpus-per-task=128
#SBATCH --error=classify-sklearn.%J.err
#SBATCH --output=classify-sklearn.%J.out




qiime feature-classifier classify-sklearn \
  --i-classifier Taxonomy/classifier-V5-799F-V7-1193R-novogene-notrim.qza \
  --i-reads rep-seqs_dada2_ee2_ee2_f222r223.qza \
  --o-classification Taxonomy/rep-seqs_dada2_ee2_ee2_f222r223_taxonomy.qza \
  --verbose

