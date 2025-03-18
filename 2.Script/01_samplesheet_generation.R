# ----------------------------------------------------------------------------
# Setup | Build nf-core sample sheets for the ATAC-seq and RNA-seq runs
#
# Paper    : Regulatory T cells in the tumor microenvironment display a
#            unique chromatin accessibility profile
# Original : 2.Script/1.1.sampleSheetGener.R
# Input    : 1.Data/ATACseq_*_sample.csv
#            1.Data/RNA_sample.csv
# Output   : 3.Output/1.AlignFiles/*_mastersampleSheet.csv
# Note     : Upstream alignment and peak calling were run with nf-core/atacseq and STAR outside this repo
# ----------------------------------------------------------------------------

if(!dir.exists('3.Output/1.AlignFiles')){
  dir.create('3.Output/1.AlignFiles')
}

library(tidyverse)
library(readxl)

# ATACseq_B16MC38CD4Treg ----
df <- read.table('1.Data/ATACseq_B16MC38CD4Treg_sample.csv')
allsmp <- gsub('R1|R2', 'END', df$V1) %>% 
  unique()

dfNew <- data.frame(fastq_1 = gsub('END', 'R1', allsmp),
                    fastq_2 = gsub('END', 'R2', allsmp),
                    sampleNAme = stringr::str_extract(allsmp, 'BD[0-9]+.*S[0-9]+')) %>% 
  mutate(tissue = ifelse(grepl('LN', fastq_1), 'LN', 'TIL'),
         tumor = ifelse(grepl('B16', fastq_1), 'B16', 'MC38'),
         cell = ifelse(grepl('CD4', fastq_1), 'CD4', 'Treg'),
         time = stringr::str_extract(fastq_1, 'Day_[0-9]+'),
         time = as.numeric(gsub('Day_', '', time)),
         tumorCellTissueTime = paste(tumor, cell, tissue, time, sep = '_'),
         tumorCellTissue = paste(tumor, cell, tissue, sep = '_'),
         tumorCell = paste(tumor, cell, sep = '_')) %>% 
  group_by(tumorCellTissue) %>% 
  mutate(replicate = 1:n())

write.csv(dfNew, '3.Output/1.AlignFiles/ATACseq_B16MC38CD4Treg_mastersampleSheet.csv', quote = F, row.names = F)


dfNew %>% 
  select(tumorCellTissue, fastq_1, fastq_2, replicate) %>% 
  rename('sample' = 1) %>% 
  write.csv('3.Output/1.AlignFiles/ATACseq_B16MC38CD4Treg_input.csv', quote = F, row.names = F)

# ATACseq_skinB16 ----
df <- read.table('1.Data/ATACseq_skinB16_sample_Val.csv')
allsmp <- gsub('_val_1|_val_2', 'END2', gsub('R1|R2', 'END1', df$V1)) %>% 
  unique()

dfNew <- data.frame(fastq_1 = gsub('END1', 'R1', allsmp) %>% gsub('END2', '_val_1', .),
                    fastq_2 = gsub('END1', 'R2', allsmp) %>% gsub('END2', '_val_2', .)) %>% 
  mutate(fqName = gsub('/ix1/dvignali/cuij/Becky/test/valfq/|_R1_001_val_1.fq.gz', '', fastq_1))

meta <- read_xlsx('1.Data/Metadata.xlsx', sheet = 2)[c(32, 96:172), ] %>% 
  janitor::row_to_names(1) %>% 
  rename('fqName' = 1)

dfNew <- merge(dfNew, meta, by = 'fqName', all = T) %>% 
  separate('*title', into = c('sample', 'rep1'), sep = '_') %>% 
  group_by(sample) %>% 
  mutate(sample = gsub(' ', '.', sample)) %>% 
  arrange(rep1) %>% 
  mutate(replicate = 1:n())


write.csv(dfNew, '3.Output/1.AlignFiles/ATACseq_skinB16_mastersampleSheet.csv', quote = F, row.names = F)

dfNew[c(4, 2, 3, 16)] %>% 
  write.csv('3.Output/1.AlignFiles/ATACseq_skinB16_input.csv', quote = F, row.names = F)

# RNA ----
df <- read.table('1.Data/RNA_sample.csv') %>% 
  filter(!grepl('Unmerged', V1))

dfNew <- data.frame(fastq_1 = df$V1,
                    fastq_2 = NA) %>% 
  mutate(fqName = gsub('/ix1/dvignali/cuij/BDData/RNA/|_R1_001.fastq.gz', '', fastq_1))

meta <- read_xlsx('1.Data/Metadata.xlsx', sheet = 1)[32:96, ] %>% 
  janitor::row_to_names(1) %>% 
  rename('fqName' = 1)

dfNew <- merge(dfNew, meta, by = 'fqName') %>% 
  separate('*title', into = c('sample', 'rep1'), sep = '_') %>% 
  group_by(sample) %>% 
  mutate(sample = gsub(' ', '.', sample)) %>% 
  arrange(rep1) %>% 
  mutate(replicate = 1:n())

write.csv(dfNew, '3.Output/1.AlignFiles/RNA_mastersampleSheet.csv', quote = F, row.names = F)

dfNew[c(4, 2, 3, 19)]  %>% 
  mutate(sample = paste(sample, replicate, sep = "_"),
         fastq_2 = '',
         strandedness = 'auto') %>% 
  select(-replicate) %>% 
  write.csv('3.Output/1.AlignFiles/RNA_input.csv', quote = F, row.names = F)
