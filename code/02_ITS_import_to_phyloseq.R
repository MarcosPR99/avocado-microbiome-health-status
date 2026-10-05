library(dplyr)
library(phyloseq)
library(speedyseq)
library(qiime2R)
library(Biostrings)

table_path <- "denoising/table_dada2_ITS_front.qza"
tax_path   <- "taxonomy/taxonomy_dada2_ITS_front.qza"
metadata   <- "metadata_octubre_its_samples.tsv"
refseq_qza <- "denoising/rep-seqs_dada2_ITS_front.qza"

physeq_ITS <- qza_to_phyloseq(
  features = table_path,
  taxonomy = tax_path,
  metadata = metadata
)

physeq_ITS

sort(sample_sums(physeq_ITS))

physeq_ITS_clean <- subset_taxa(physeq_ITS, Kingdom != "Unassigned" | is.na(Kingdom))
saveRDS(physeq_ITS_clean, file = "Phyloseq_ITS.rds")

refseq <- read_qza(refseq_qza)
refseq <- refseq$data
writeXStringSet(refseq, "refseq_ITS_from_qza.fasta", format = "fasta")
