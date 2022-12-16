setwd("/Users/anttonalberdi/github/buffer_metagenomics_test/")

#########
# Prepare count table
#########

count_raw <- read.table("rawdata/buffer_test_final_count_table.txt",header=T,row.names=1,sep="\t")
sample_metadata <- read.table("data/sample_metadata.tsv",header=T,sep="\t")

#Remove sample suffices
colnames(count_raw) <- gsub(".Read.Count|.Covered.Fraction|.Length","",colnames(count_raw))
#Rename
colnames(count_raw) <- sample_metadata$Sample[match(colnames(count_raw), gsub("-","\\.",sample_metadata$Data))]

#Split
counts <- count_raw[ ,seq(from = 1, to = ncol(count_raw), by = 3)]  # Reads
fraction <- count_raw[ ,seq(from = 2, to = ncol(count_raw), by = 3)]  # Fraction
mag_length <- as.data.frame(t(t(count_raw[ ,seq(from = 3, to = ncol(count_raw), by = 3)][,1])))  # Length
colnames(mag_length) <- "Length"
mag_length$Genome <- rownames(fraction)

#Save
write.table(counts,"data/counts.tsv",col.names=T,row.names=T,sep="\t",quote=F)
write.table(fraction,"data/fraction.tsv",col.names=T,row.names=T,sep="\t",quote=F)

#########
# Prepare MAG attributes
#########

summary_file <- read.table("rawdata/buffer_test_combined_summary.tsv",header=T,sep="\t")

library(dplyr)
library(tidyr)
library(stringr)

# Split name column into firstname and last name
summary_file %>%
  select(user_genome, classification) %>%
  rename(Genome = user_genome) %>%
  separate(classification, c('Domain', 'Phylum','Class','Order','Family','Genus','Species'),sep=";") %>%
  mutate_at(vars(Domain,Phylum,Class,Order,Family,Genus,Species), ~ str_replace(., "[dpcofgs]__", "")) %>%
  mutate_at(vars(Genome), ~ str_replace(., ".fa", "")) %>%
  inner_join(mag_length,by="Genome") %>%
  write.table(.,"data/genome_metadata.tsv",col.names=T,row.names=F,sep="\t",quote=F)
