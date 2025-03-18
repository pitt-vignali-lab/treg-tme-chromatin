# ----------------------------------------------------------------------------
# Setup | Merge sample metadata across the ATAC-seq batches
#
# Paper    : Regulatory T cells in the tumor microenvironment display a
#            unique chromatin accessibility profile
# Original : 2.Script/1.2.metaDataIntgration.R
# Input    : 1.Data/Metadata.xlsx
#            1.Data/ATACseq_*_sample.csv
# Output   : 1.Data/ATAC_overall_sample_Infor.RDS
# Upstream : 01_samplesheet_generation.R
# ----------------------------------------------------------------------------

library(tidyverse)
library(DESeq2)
library(ggthemes)
library(readxl)

count <- data.table::fread('1.Data/ATAC_overall/consensus_peaks.mLb.clN.featureCounts.txt')[, -c(2:6)] %>% 
  column_to_rownames('Geneid')

names(count) <- gsub('.mLb.clN.sorted.bam', '', names(count))

mcbtregcd <- read.csv('3.Output/1.AlignFiles/ATACseq_B16MC38CD4Treg_mastersampleSheet.csv') %>% 
  mutate(samples = paste0(tumorCellTissue, '_REP', replicate),
         sampleNAme = gsub('_Day.*', '', sampleNAme)) %>% 
  select(sampleNAme, time, samples) %>% 
  dplyr::rename('wkinName' = 1)

sknpdcd <-  read.csv('3.Output/1.AlignFiles/ATACseq_skinB16_mastersampleSheet.csv') %>% 
  mutate(samples = paste0(sample, '_REP', replicate),
         fqName = gsub('_S[0-9]+|_VIG.*|Lib_', '', fqName),
         wkinName = gsub('_', '', fqName),
         time = NA) %>% 
  select(wkinName, time, samples)

smplDf <- data.frame(samples = names(count)) %>% 
  mutate(cell = ifelse(grepl('CD4', samples), 'CD4', 'Treg'),
         tissue = case_when(
           grepl('TIL', samples) ~ 'TIL',
           grepl('LN', samples) ~ 'LN',
           grepl('Spl', samples) ~ 'Spleen',
           grepl('Skin', samples) ~ 'Skin',
           grepl('rTr', samples) ~ 'Resting',
           grepl('VAT', samples) ~ 'VAT',
           grepl('aTr', samples) ~ 'Activated',
           grepl('Muscle', samples) ~ 'Muscle',
           grepl('Colon', samples) ~ 'Colon'
         ),
         project = case_when(
           grepl('B16|MC38', samples) ~ 'MC38B16TregCD4',
           grepl('Iso|Late|Early|NR', samples) ~ 'PDCD1',
           grepl('CD4.LN_REP|CD4.Skin_REP|Treg.Skin_REP|Treg.LN_REP', samples) ~ 'Skin',
           TRUE ~ 'GEO'
         ),
         replicate = as.numeric(gsub("^.*REP", "", samples))) %>% 
  merge(rbind(mcbtregcd, sknpdcd), by = 'samples', all = T)

saveRDS(smplDf, '1.Data/ATAC_overall_sample_Infor.RDS')
