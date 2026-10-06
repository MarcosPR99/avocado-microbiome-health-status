# Analysis workflow

This folder contains the scripts used to process and analyse bacterial 16S rRNA gene and fungal ITS amplicon datasets associated with the manuscript.

Raw FASTQ files are not included in this repository. Raw reads are available through the NCBI Sequence Read Archive under BioProject accession PRJNA1540021.

## 16S rRNA gene workflow

```text
Raw paired-end FASTQ files
→ QIIME2 import using metadata_samples.tsv as manifest
→ initial demultiplexed read quality summary
→ primer trimming with cutadapt
→ post-trimming quality summary
→ DADA2 denoising
→ MAFFT alignment
→ alignment masking
→ FastTree phylogenetic tree construction
→ midpoint rooting
→ taxonomic classification
→ import into R/phyloseq
→ downstream diversity, differential abundance and network analyses
```

Scripts:

```text
00_16S_import_quality_trimming.sh
01a_16S_dada2_denoising.sh
01b_16S_alignment_mafft.sh
01c_16S_alignment_mask.sh
01d_16S_fasttree_rooted.sh
01e_16S_taxonomy_classification.sh
03_16S_downstream_analysis.R
```

## ITS workflow

The ITS workflow followed the same general QIIME2 structure, adapted to ITS file names and marker-specific processing.

```text
Raw paired-end FASTQ files
→ QIIME2 import using the ITS metadata/manifest file
→ initial demultiplexed read quality summary
→ primer/adaptor trimming with cutadapt
→ post-trimming quality summary
→ DADA2 denoising
→ taxonomic classification
→ import into R/phyloseq
→ downstream diversity, differential abundance and network analyses
```

Scripts:

```text
00_ITS_import_quality_trimming.sh
02_ITS_import_to_phyloseq.R
04_ITS_downstream_analysis.R
```

The ITS preprocessing script documents the same import, quality summary and trimming structure used for 16S. Primer sequences should be checked against the original ITS run configuration before re-running from raw FASTQ files.

## Notes

Some scripts were originally run in local or high-performance computing environments. File paths may need to be adapted before re-running the analyses on a different system.

Large intermediate files such as FASTQ, `.qza`, `.qzv`, `.rds` and `.RDS` files are intentionally not tracked in this repository.
