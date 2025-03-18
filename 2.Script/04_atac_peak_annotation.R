# ----------------------------------------------------------------------------
# Fig 2B | ChIPseeker peak annotation and genomic feature distribution
#
# Paper    : Regulatory T cells in the tumor microenvironment display a
#            unique chromatin accessibility profile
# Original : 2.Script/3.annotation.R
# Input    : 3.Output/2.DEseqRes/1.B16MC38TregCD4/de_*.RDS
# Output   : 3.Output/3.PeakAnno/peakAnnotationChipseeker.RDS
#            3.Output/3.PeakAnno/de_anno_Treg_*.csv
#            3.Output/3.PeakAnno/Pie_*.png
# Upstream : 03_atac_deseq_tumor_vs_ln.R
# Note     : Annotated peak tables feed every downstream enrichment script
# ----------------------------------------------------------------------------

if(!dir.exists('3.Output/3.PeakAnno')){
  dir.create('3.Output/3.PeakAnno')
}

library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(tidyverse)

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

# functional region annotation 
df <- read.table('1.Data/ATAC_overall/consensus_peaks.mLb.clN.featureCounts.txt', header = T)
  
peakAll <- GRanges(seqnames = df$Chr, 
                   ranges = IRanges(as.numeric(df$Start), 
                                    as.numeric(df$End)),
                   peakID = df$Geneid)

peakAnno <- annotatePeak(peakAll, 
                         tssRegion=c(-3000, 3000),
                         TxDb=txdb, 
                         annoDb="org.Mm.eg.db")


saveRDS(peakAnno, '3.Output/3.PeakAnno/peakAnnotationChipseeker.RDS')

peakAnno <- readRDS('3.Output/3.PeakAnno/peakAnnotationChipseeker.RDS')
annoDf <- data.frame(peakAnno@anno)

# subsection anno 
readRDS('3.Output/2.DEseqRes/1.B16MC38TregCD4/subsection.RDS') %>% 
  merge(annoDf, by.x = 0, by.y = 'peakID') %>%
  write.csv('3.Output/3.PeakAnno/subsectionAnnotated.csv', quote = F)

for (i in c('B16', 'MC38')) {
  de <- readRDS(paste0('3.Output/2.DEseqRes/1.B16MC38TregCD4/de_', i, '_Treg.RDS')) %>% 
    as.data.frame() %>% 
    merge(annoDf, by.x = 0, by.y = 'peakID') %>% 
    dplyr::rename('peakID' = 1) %>% 
    mutate(peakID = gsub('Interval', 'Peak', peakID))
  
   de %>% filter(log2FoldChange < 0, pvalue < 0.05) %>% 
    write.csv(paste0('3.Output/3.PeakAnno/de_anno_Treg_', i, '_(p0.05)down.csv'))
   
   de %>% filter(log2FoldChange > 0, pvalue < 0.05) %>% 
     write.csv(paste0('3.Output/3.PeakAnno/de_anno_Treg_', i, '_(p0.05)up.csv'))
   
   de %>% filter(log2FoldChange < 0, padj < 0.05) %>% 
     write.csv(paste0('3.Output/3.PeakAnno/de_anno_Treg_', i, '_(padj0.05)down.csv'))
   
   de %>% filter(log2FoldChange > 0, padj < 0.05) %>% 
     write.csv(paste0('3.Output/3.PeakAnno/de_anno_Treg_', i, '_(padj0.05)up.csv'))
   
  
}




# generate Pie plot
df <- read.table('1.Data/ATAC_overall/consensus_peaks.mRp.clN.boolean.txt', header = T) 

for (i in c('B16_Treg_TIL', 'B16_Treg_LN', 'MC38_Treg_TIL', 'MC38_Treg_LN')) {
  dfSub <- df[, c(1:6, grep(i, names(df)))] %>% 
    filter(eval(parse(text = paste0(i, '.mRp.clN.bool'))))
  
  subAnno <- annotatePeak(GRanges(seqnames = dfSub$chr, 
                                  ranges = IRanges(as.numeric(dfSub$start), 
                                                   as.numeric(dfSub$end)),
                                  peakID = dfSub$interval_id), 
                          tssRegion=c(-3000, 3000),
                          TxDb=txdb, 
                          annoDb="org.Mm.eg.db")
  
  png(paste0('3.Output/3.PeakAnno/Pie_', i, '.png'), width = 7, height = 5, units = 'in', res = 300)
  plotAnnoPie(subAnno, col = RColorBrewer::brewer.pal(12, 'Paired')[c(1:6, 8:12)])
  dev.off()
  
  pdf(paste0('3.Output/3.PeakAnno/Pie_', i, '.pdf'), width = 7, height = 5)
  plotAnnoPie(subAnno, col = RColorBrewer::brewer.pal(12, 'Paired')[c(1:6, 8:12)])
  dev.off()
}