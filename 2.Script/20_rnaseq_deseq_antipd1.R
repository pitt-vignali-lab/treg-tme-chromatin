# ----------------------------------------------------------------------------
# Fig 4C, 4E, 4G | RNA-seq differential expression after anti-PD-1 or isotype treatment
#
# Paper    : Regulatory T cells in the tumor microenvironment display a
#            unique chromatin accessibility profile
# Original : 2.Script/2.5.1.RNAseq_DEseq.R
# Input    : 1.Data/RNA_readTabls/*ReadsPerGene.out.tab
#            1.Data/RNA_sample.csv
# Output   : 1.Data/RNAseq_count.RDS
#            3.Output/2.DEseqRes/2.4.RNAseq/RNAseq.pdf
# Upstream : 01_samplesheet_generation.R
# ----------------------------------------------------------------------------

if(!dir.exists('3.Output/2.DEseqRes/2.4.RNAseq')){
  dir.create('3.Output/2.DEseqRes/2.4.RNAseq')
}

library(tidyverse)
library(DESeq2)

tab_files<-list.files("1.Data/RNA_readTabls/",pattern = ".tab", full.names = T)

count <- do.call(cbind,
                 lapply(tab_files, 
                        function(x) {
                          tmp <- read.table(x)[-c(1:4), c(1, 4)] %>% 
                            remove_rownames() %>% 
                            column_to_rownames('V1')
                          
                          names(tmp) <- gsub('1.Data/RNA_readTabls/|ReadsPerGene.out.tab', '', x)
                          
                          return(tmp)
                          
                        }))

saveRDS(count, '1.Data/RNAseq_count.RDS')

sampleDf <- data.frame(samples = names(count)) %>% 
  separate(samples, into = c('tissue', 'type', 'cell', 'replicate'), remove = F) %>% 
  mutate(type2 = ifelse(type %in% c('EARLY', 'LATE'), 'Res', type),
         cllTssTp = paste0(cell, tissue, type2))

saveRDS(sampleDf, '1.Data/RNAseq_meta.RDS')

#### all groups ####
dds <- DESeqDataSetFromMatrix(countData = count[,sampleDf$samples], 
                              colData = sampleDf,
                              design = ~ cllTssTp)

dds <- DESeq(dds)

saveRDS(dds, '3.Output/2.DEseqRes/2.4.RNAseq/RNAseq.RDS')

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

ggsave('3.Output/2.DEseqRes/2.4.RNAseq/RNAseq.png', 
       width = 6, height = 6, units = 'in')
ggsave('3.Output/2.DEseqRes/2.4.RNAseq/RNAseq.pdf', 
       width = 6, height = 6)

### Treg Til ################
subSmp <- sampleDf %>% 
  filter(cell == 'TREG',
         tissue == 'TIL') 

dds <- DESeqDataSetFromMatrix(countData = count[, subSmp$sample], 
                              colData = subSmp,
                              design = ~ type2)

dds <- DESeq(dds)

saveRDS(dds, '3.Output/2.DEseqRes/2.4.RNAseq/RNAseq_TregTIL.RDS')

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

ggsave('3.Output/2.DEseqRes/2.4.RNAseq/pca_RNAseq_TregTIL.png', 
       width = 6, height = 6, units = 'in')
ggsave('3.Output/2.DEseqRes/2.4.RNAseq/pca_RNAseq_TregTIL.pdf', 
       width = 6, height = 6)

#### CD4 TIL ###########
subSmp <- sampleDf %>% 
  filter(cell == 'CD4',
         tissue == 'TIL')

dds <- DESeqDataSetFromMatrix(countData = count[, subSmp$sample], 
                              colData = subSmp,
                              design = ~ type2)

dds <- DESeq(dds)

saveRDS(dds, '3.Output/2.DEseqRes/2.4.RNAseq/RNAseq_CD4TIL.RDS')

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

ggsave('3.Output/2.DEseqRes/2.4.RNAseq/pca_RNAseq_CD4TIL.png', 
       width = 6, height = 6, units = 'in')
ggsave('3.Output/2.DEseqRes/2.4.RNAseq/pca_RNAseq_CD4TIL.pdf', 
       width = 6, height = 6)
