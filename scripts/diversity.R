setwd("/Users/anttonalberdi/github/buffer_metagenomics_test/")

#########
# Load tables
#########

counts <- read.table("data/counts.tsv",header=T,row.names=1,sep="\t")
sample_metadata <- read.table("data/sample_metadata.tsv",header=T,row.names=1,sep="\t")
genome_metadata <- read.table("data/genome_metadata.tsv",header=T,row.names=1,sep="\t")

#########
# Assess sequencing depth
#########



#########
# Run diversity analyses
#########

library(hilldiv)
library(ggplot2)

hillq1 <- hill_div(tss(counts),qvalue=1)
hillq2 <- hill_div(tss(counts),qvalue=2)

#########
# Plots
#########

#Diversity
div_test_plot(div_test(tss(counts),qvalue=1,hierarchy=sample_metadata))
div_test_plot(div_test(tss(counts),qvalue=2,hierarchy=sample_metadata))

#Depth
cbind(Depth=colSums(counts),sample_metadata) %>%
  ggplot(aes(.,x=Treatment,y=Depth)) +
  geom_boxplot() +
  theme_light()

#Depth vs. diversity
as.data.frame(cbind(Depth=colSums(counts),hillq1,sample_metadata)) %>%
  ggplot(aes(.,x=Depth,y=hillq1,color=Treatment)) +
  geom_point() +
  theme_light()

#Not enough structure for ordination
pairq1 <- pair_dis(tss(counts),qvalue=1,hierarchy=sample_metadata)
pair_dis_plot(pairq1$L1_VqN)
