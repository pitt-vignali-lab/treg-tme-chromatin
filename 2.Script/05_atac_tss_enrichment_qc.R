# ----------------------------------------------------------------------------
# Suppl Fig 1A-B | TSS enrichment scores and ATAC-seq mapping statistics
#
# Paper    : Regulatory T cells in the tumor microenvironment display a
#            unique chromatin accessibility profile
# Original : 2.Script/3.TSSenrch.R
# Input    : 1.Data/bams_ATAC_B16MC38CD4Treg/*.mLb.clN.sorted.bam
# Output   : 3.Output/3.PeakAnno/allATACseqTSSes.csv
#            3.Output/1.AlignFiles/ATACseqMapStat.csv
# Note     : Quality control
#            requires the aligned BAM files
# ----------------------------------------------------------------------------

library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(tidyverse)
library(GenomicAlignments)
library(GenomicRanges)
library(ATACseqQC)

txs <- transcripts(TxDb.Mmusculus.UCSC.mm10.knownGene)

bams <- list.files('1.Data/bams_ATAC_B16MC38CD4Treg/',
                   pattern = 'bam$', full.names = T)

gals <- lapply(bams, \(x) readBamFile(bamFile=x, bigFile=TRUE))
names(gals) <- gsub('1.Data/bams_ATAC_B16MC38CD4Treg/|.mLb.clN.sorted.bam', 
                    '', bams)

tsses <- sapply(gals, \(x) return(TSSEscore(x, txs)$TSSEscore))

data.frame(tsses) %>% 
  write.csv('3.Output/3.PeakAnno/allATACseqTSSes.csv', quote = F)

df <- read.csv('3.Output/3.PeakAnno/allATACseqTSSes.csv')

samplesheet <- read.csv('3.Output/1.AlignFiles/overall_input.csv')[-c(1:20),] %>% 
  mutate(X = paste0(sample, '_REP', replicate),
         fastq_1 = gsub('/ix1/dvignali/cuij/BDData/ATAC|.fastq.gz|_R1_001|_R1',
                        '', fastq_1)) %>% 
  select(-fastq_2)

sumDf <- data.table::fread('3.Output/3.PeakAnno/multiqc_picard_AlignmentSummaryMetrics.txt')[,c(1, 3 ,10)] %>% 
  mutate(ratio = PF_HQ_ALIGNED_READS*100/TOTAL_READS) %>% 
  merge(df, by.x = 'Sample', by.y = 'X') %>% 
  merge(samplesheet, by.x = 'Sample', by.y = 'X') %>% 
  separate(fastq_1, into = c('proj', 'sample_id', 'other', 'other1'), 
           remove = F, sep = '/', fill = 'left') %>% 
  mutate(other1 = gsub('_', '', 
                       gsub('_VIG.*|_S[0-9]+|_Day.*', '', 
                       paste0('BD', gsub('.*BD', '', other1)))),
         cell = ifelse(grepl('CD4', Sample), 'CD4', 'Treg'))

write.csv(sumDf, '3.Output/1.AlignFiles/ATACseqMapStat.csv')
