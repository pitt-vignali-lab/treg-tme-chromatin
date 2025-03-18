# ----------------------------------------------------------------------------
# Fig 4B, 4D, 4F | Same ATAC-seq comparison after ComBat-seq batch correction
#
# Paper    : Regulatory T cells in the tumor microenvironment display a
#            unique chromatin accessibility profile
# Original : 2.Script/2.2.2.B16pdcd_DEseq_BatEffCor.R
# Input    : 1.Data/skinB16_pdcd/consensus_peaks.mLb.clN.featureCounts.txt
# Output   : 3.Output/2.DEseqRes/2.B16pdcd/pca_B16pdcd_BatEffCor.pdf
# Upstream : 18_atac_deseq_antipd1.R
# Note     : Batch-corrected counterpart of the previous script
# ----------------------------------------------------------------------------

library(readxl)
library(tidyverse)
library(DESeq2)
library(ggthemes)
library(export)
library(sva)


count <- data.table::fread('1.Data/ATAC_overall/consensus_peaks.mLb.clN.featureCounts.txt')[, -c(2:6)] %>% 
  column_to_rownames('Geneid')

names(count) <- gsub('.mLb.clN.sorted.bam', '', names(count))


sampleDf <- readRDS('1.Data/ATAC_overall_sample_Infor.RDS') %>% 
  filter(project == 'PDCD1') %>% 
  separate(samples, into = c('cell', 'tissu', 'type', 'rep'), remove = F) %>% 
  mutate(type2 = ifelse(type %in% c('Early', 'Late'), 'Res', type),
         rep = gsub('REP', '', rep),
         cllTssTp = paste0(cell, tissu, type2))


#### all groups ####
sampleDf <- sampleDf %>% 
  merge( read_xlsx('1.Data/Subinoy ATACseq analysis_CW and BD (1).xlsx', sheet = 2) %>% 
          select(sample, batch) , by.y = 'sample', by.x = 'wkinName')

cntNew <- ComBat_seq(as.matrix(count[, sampleDf$samples]), 
                     batch = sampleDf$batch,
                     group = sampleDf$cllTssTp)

dds <- DESeqDataSetFromMatrix(countData = cntNew, 
                              colData = sampleDf,
                              design = ~ cllTssTp)

dds <- DESeq(dds)

saveRDS(dds, '3.Output/2.DEseqRes/2.B16pdcd/B16pdcd_BatEffCor.RDS')

vsd <- vst(dds, blind=FALSE)

pcaData <- plotPCA(vsd, intgroup=c('cllTssTp'), returnData=TRUE)

percentVar <- round(100 * attr(pcaData, "percentVar"))

pcaData %>% 
  merge(sampleDf, by.x = 'name', by.y = 'samples') %>% 
  mutate(cllTss = paste(cell, tissu, sep = '_')) %>% 
  ggplot(aes(PC1, PC2, shape = type2, fill = cllTss)) +
  # geom_text(aes(label = name))
  geom_point(size = 5, col = 'black') +
  scale_shape_manual(values = 21:23) +
  scale_fill_brewer(palette = 'Paired') +
  scale_fill_manual(values = alpha(c("#fed589", "#FF5927", "#A6CEE3", "#013e6f"),
                                   alpha = 0.8)) +
  guides(fill = guide_legend(override.aes=list(shape=21))) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  theme(# legend.position = 'none',
    aspect.ratio=1,
    text = element_text(size = 15, face = 2),
    panel.grid = element_blank(),
    axis.text = element_text(size = 15, face = 2, color = 'black'),
    legend.key = element_blank())

ggsave('3.Output/2.DEseqRes/2.B16pdcd/pca_B16pdcd_BatEffCor.png', 
       width = 6, height = 6, units = 'in')
ggsave('3.Output/2.DEseqRes/2.B16pdcd/pca_B16pdcd_BatEffCor.pdf', 
       width = 6, height = 6)

### Treg Til ################
subSmp <- sampleDf %>% 
  filter(cell == 'Treg',
         tissu == 'TIL')

cntNew <- ComBat_seq(as.matrix(count[, subSmp$sample]), 
                     batch = subSmp$batch,
                     group = subSmp$cllTssTp)

dds <- DESeqDataSetFromMatrix(countData = cntNew, 
                              colData = subSmp,
                              design = ~ type2)

dds <- DESeq(dds)

saveRDS(dds, '3.Output/2.DEseqRes/2.B16pdcd/B16pdcd_TregTIL_BatEffCor.RDS')

vsd <- vst(dds, blind=FALSE)

pcaData <- plotPCA(vsd, intgroup=c('cllTssTp'), returnData=TRUE)

percentVar <- round(100 * attr(pcaData, "percentVar"))

pcaData %>% 
  merge(subSmp, by.x = 'name', by.y = 'samples') %>% 
  mutate(type2 = ifelse(type %in% c('Early', 'Late'), 'Resp', type),
         cllTss = paste(cell, tissu, sep = '_')) %>% 
  ggplot(aes(PC1, PC2, shape = type2)) +
  geom_point(size = 5, col = 'black', fill = alpha("#013e6f", 0.8)) +
  scale_shape_manual(values = 21:23) +
  scale_fill_brewer(palette = 'Paired') +
  # scale_fill_manual(values = alpha(c("#FDBF6F", "#FF7F00", "#A6CEE3", "#1F78B4"),
  #                                  alpha = 0.8)) +
  guides(fill = guide_legend(override.aes=list(shape=21))) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  # coord_fixed() +
  theme(# legend.position = 'none',
    aspect.ratio=1,
    text = element_text(size = 15, face = 2),
    panel.grid = element_blank(),
    axis.text = element_text(size = 15, face = 2, color = 'black'),
    legend.key = element_blank())

ggsave('3.Output/2.DEseqRes/2.B16pdcd/pca_B16pdcd_TregTIL_BatEffCor.png', 
       width = 6, height = 6, units = 'in')
ggsave('3.Output/2.DEseqRes/2.B16pdcd/pca_B16pdcd_TregTIL_BatEffCor.pdf', 
       width = 6, height = 6)

#### CD4 TIL ###########
subSmp <- sampleDf %>% 
  filter(cell == 'CD4',
         tissu == 'TIL')

cntNew <- ComBat_seq(as.matrix(count[, subSmp$sample]), 
                     batch = subSmp$batch,
                     group = subSmp$cllTssTp)

dds <- DESeqDataSetFromMatrix(countData = cntNew, 
                              colData = subSmp,
                              design = ~ type2)

dds <- DESeq(dds)

saveRDS(dds, '3.Output/2.DEseqRes/2.B16pdcd/B16pdcd_CD4TIL_BatEffCor.RDS')

vsd <- vst(dds, blind=FALSE)

pcaData <- plotPCA(vsd, intgroup=c('cllTssTp'), returnData=TRUE)

percentVar <- round(100 * attr(pcaData, "percentVar"))

pcaData %>% 
  merge(subSmp, by.x = 'name', by.y = 'samples') %>% 
  mutate(type2 = ifelse(type %in% c('Early', 'Late'), 'Resp', type),
         cllTss = paste(cell, tissu, sep = '_')) %>% 
  ggplot(aes(PC1, PC2, shape = type2)) +
  geom_point(size = 5, col = 'black', fill = alpha("#FF5927", 0.8)) +
  scale_shape_manual(values = 21:23) +
  scale_fill_brewer(palette = 'Paired') +
  # scale_fill_manual(values = alpha(c("#FDBF6F", "#FF7F00", "#A6CEE3", "#1F78B4"),
  #                                  alpha = 0.8)) +
  guides(fill = guide_legend(override.aes=list(shape=21))) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  # coord_fixed() +
  theme(# legend.position = 'none',
    aspect.ratio=1,
    text = element_text(size = 15, face = 2),
    panel.grid = element_blank(),
    axis.text = element_text(size = 15, face = 2, color = 'black'),
    legend.key = element_blank())

ggsave('3.Output/2.DEseqRes/2.B16pdcd/pca_B16pdcd_CD4TIL_BatEffCor.png', 
       width = 6, height = 6, units = 'in')
ggsave('3.Output/2.DEseqRes/2.B16pdcd/pca_B16pdcd_CD4TIL_BatEffCor.pdf', 
       width = 6, height = 6)
