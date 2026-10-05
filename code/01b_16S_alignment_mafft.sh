#!/bin/bash
#SBATCH --mem=1800gb
#SBATCH --cpus-per-task=128
#SBATCH --error=alignment.%J.err
#SBATCH --output=alignment.%J.out



qiime alignment mafft \
  --i-sequences rep-seqs_dada2_ee2_ee2_f222r223.qza \
  --p-n-threads 128 \
  --o-alignment Phylogenetic_trees/ee2_ee2_f222r223_alignment.qza