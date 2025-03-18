# ----------------------------------------------------------------------------
# Suppl Table | Collect STAR mapping statistics for the RNA-seq libraries
#
# Paper    : Regulatory T cells in the tumor microenvironment display a
#            unique chromatin accessibility profile
# Original : 2.Script/3.RNAseqState.R
# Input    : 1.Data/STARlog/
# Output   : 3.Output/1.AlignFiles/RNAseqMapStat.csv
# Note     : Quality control
# ----------------------------------------------------------------------------

library(tidyverse)

fls <- list.files('1.Data/STARlog/', full.names = T)

df <- data.frame()
for (i in seq(length(fls))) {
  df <- rbind(df, readLines(fls[i])[c(6, 9, 10)] %>% 
    data.frame() %>% 
    rename('t' = 1) %>% 
    separate(t, into = c('cat', 
                          gsub("1.Data/STARlog/|Log.final.out", '', fls[i])), 
             sep = ' \\|\t') %>% 
    mutate(cat = gsub('^ +| +$', '', cat)) %>% 
    column_to_rownames('cat') %>% 
    t() %>% 
    data.frame())
}

samplesheet <- read.csv('3.Output/1.AlignFiles/RNA_input.csv')

dfNew <- merge(df, samplesheet, by.x = 0, by.y = 'sample') %>% 
  mutate(fastq_1 = gsub('/ix1/dvignali/cuij/BDData/RNA/|_R1_001.fastq.gz',
                        '',
                        fastq_1),
         Row.names = gsub('_[0-9]+', '', Row.names),
         fastq_1 = paste0('BD', gsub('.*BD', '', fastq_1)),
         fastq_1 =  gsub('_S[0-9]+', '', fastq_1),
         fastq_1 =  ifelse(grepl('[0-9]b', fastq_1),
                           gsub('b', '', gsub('_', '_0', fastq_1)),
                           fastq_1),
         fastq_1 = ifelse(grepl('BDHY', fastq_1),
                          fastq_1, gsub('_', '', fastq_1))) %>% 
  separate('Row.names', into = c('tissue', 'type', 'cell'), remove = F) 

write.csv(dfNew, '3.Output/1.AlignFiles/RNAseqMapStat.csv')
