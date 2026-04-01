library(readxl)
library(dplyr)
setwd("/lisc/scratch/jmf/internal/Analyses_Jay/JMF_2110_03/results/JMF-2110-03-0008_spades/final_bins/dereplicated_genomes")

table <- read_xlsx("final_QC.xlsx")

table2 <- table %>% mutate_all(~gsub(".fa", "", .))
table2 <- table2 %>% rename(Genome = genome)


coverM_mean <- read.delim("coverM_mean_output.tsv")
coverM_relat <- read.delim("coverM_relative_output.tsv")

table_full <- table2 %>% left_join(coverM_mean) %>% left_join(coverM_relat)

names(table_full) = gsub(pattern = ".interleave.fastq.gz.", replacement = " ", x = names(table_full))

writexl::write_xlsx(table_full, "JMF-21103-0008_HighQ_stats.xlsx")
