#!/bin/bash
#SBATCH --mem=1800gb
#SBATCH --cpus-per-task=128
#SBATCH --error=fasttree.%J.err
#SBATCH --output=fasttree.%J.out


qiime phylogeny midpoint-root \
  --i-tree Phylogenetic_trees/ee2_ee2_f222r223_fasttree.qza \
  --o-rooted-tree Phylogenetic_trees/ee2_ee2_f222r223_fasttree_rooted.qza

