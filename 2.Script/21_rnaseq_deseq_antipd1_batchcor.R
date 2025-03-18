# ----------------------------------------------------------------------------
# Fig 4C, 4E, 4G | Same RNA-seq comparison after ComBat-seq batch correction
#
# Paper    : Regulatory T cells in the tumor microenvironment display a
#            unique chromatin accessibility profile
# Original : 2.Script/2.5.2.RNAseq_DEseq_batEffCort.R
# Input    : 1.Data/RNAseq_count.RDS
#            1.Data/RNAseq_meta.RDS
# Output   : 3.Output/2.DEseqRes/5.RNAseq/RNAseq_batEffCorct.pdf
# Upstream : 20_rnaseq_deseq_antipd1.R
# Note     : Batch-corrected counterpart of the previous script
# ----------------------------------------------------------------------------

if(!dir.exists('3.Output/2.DEseqRes/5.RNAseq')){
  dir.create('3.Output/2.DEseqRes/5.RNAseq')
}

library(tidyverse)
library(DESeq2)
library(readxl)
library(sva)

count <- readRDS('1.Data/RNAseq_count.RDS')

sampleDf <- readRDS('1.Data/RNAseq_meta.RDS')

#### all groups ####

cntNew <- ComBat_seq(as.matrix(count[, sampleDf$samples]), 
                     batch = sampleDf$batch,
                     group = sampleDf$cllTssTp)

dds <- DESeqDataSetFromMatrix(countData = cntNew, 
                              colData = sampleDf,
                              design = ~ cllTssTp)

dds <- DESeq(dds)

saveRDS(dds, '3.Output/2.DEseqRes/5.RNAseq/RNAseq_batEffCorct.RDS')

vsd <- vst(dds, blind=FALSE)

pcaData <- plotPCA(vsd, intgroup=c('cllTssTp'), returnData=TRUE)

percentVar <- round(100 * attr(pcaData, "percentVar"))

pcaData %>% 
  merge(sampleDf, by.x = 'name', by.y = 'samples') %>% 
  mutate(cllTss = paste(cell, tissue, sep = '_')) %>% 
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
  labs(fill = 'Cell',
       shape = 'Type') +
  theme(# legend.position = 'none',
    aspect.ratio=1,
    text = element_text(size = 15, face = 2),
    panel.grid = element_blank(),
    axis.text = element_text(size = 15, face = 2, color = 'black'),
    legend.key = element_blank())

ggsave('3.Output/2.DEseqRes/5.RNAseq/RNAseq_batEffCorct.png', 
       width = 6, height = 6, units = 'in')
ggsave('3.Output/2.DEseqRes/5.RNAseq/RNAseq_batEffCorct.pdf', 
       width = 6, height = 6)

### Treg Til ################
subSmp <- sampleDf %>% 
  filter(cell == 'TREG',
         tissue == 'TIL') 

cntNew <- ComBat_seq(as.matrix(count[, subSmp$sample]), 
                     batch = subSmp$batch,
                     group = subSmp$cllTssTp)

dds <- DESeqDataSetFromMatrix(countData = cntNew, 
                              colData = subSmp,
                              design = ~ type2)
dds <- DESeq(dds)

saveRDS(dds, '3.Output/2.DEseqRes/5.RNAseq/RNAseq_TregTIL_batEffCorct.RDS')

vsd <- vst(dds, blind=FALSE)

pcaData <- plotPCA(vsd, intgroup=c('cllTssTp'), returnData=TRUE)

percentVar <- round(100 * attr(pcaData, "percentVar"))

pcaData %>% 
  merge(subSmp, by.x = 'name', by.y = 'samples') %>% 
  mutate(cllTss = paste(cell, tissue, sep = '_')) %>% 
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

ggsave('3.Output/2.DEseqRes/5.RNAseq/pca_RNAseq_TregTIL_batEffCorct.png', 
       width = 6, height = 6, units = 'in')
ggsave('3.Output/2.DEseqRes/5.RNAseq/pca_RNAseq_TregTIL_batEffCorct.pdf', 
       width = 6, height = 6)

#### CD4 TIL ###########
subSmp <- sampleDf %>% 
  filter(cell == 'CD4',
         tissue == 'TIL')

cntNew <- ComBat_seq(as.matrix(count[, subSmp$sample]), 
                     batch = subSmp$batch,
                     group = subSmp$cllTssTp)

dds <- DESeqDataSetFromMatrix(countData = cntNew, 
                              colData = subSmp,
                              design = ~ type2)

dds <- DESeq(dds)

saveRDS(dds, '3.Output/2.DEseqRes/5.RNAseq/RNAseq_CD4TIL_batEffCorct.RDS')

vsd <- vst(dds, blind=FALSE)

pcaData <- plotPCA(vsd, intgroup=c('cllTssTp'), returnData=TRUE)

percentVar <- round(100 * attr(pcaData, "percentVar"))

pcaData %>% 
  merge(subSmp, by.x = 'name', by.y = 'samples') %>% 
  mutate(cllTss = paste(cell, tissue, sep = '_')) %>% 
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

ggsave('3.Output/2.DEseqRes/5.RNAseq/pca_RNAseq_CD4TIL_batEffCorct.png', 
       width = 6, height = 6, units = 'in')
ggsave('3.Output/2.DEseqRes/5.RNAseq/pca_RNAseq_CD4TIL_batEffCorct.pdf', 
       width = 6, height = 6)
