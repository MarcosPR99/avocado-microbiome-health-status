#!/bin/bash
#SBATCH --mem=1800gb
#SBATCH --cpus-per-task=128
#SBATCH --error=alignment_mask.%J.err
#SBATCH --output=alignment_mask.%J.out

qiime alignment mask \
  --i-alignment Phylogenetic_trees/ee2_ee2_f222r223_alignment.qza \
  --o-masked-alignment Phylogenetic_trees/ee2_ee2_f222r223_alignment_masked.qza