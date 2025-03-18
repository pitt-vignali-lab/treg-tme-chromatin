# ----------------------------------------------------------------------------
# Fig 3B | Joint PCA of skin, B16 and MC38 Tregs
#
# Paper    : Regulatory T cells in the tumor microenvironment display a
#            unique chromatin accessibility profile
# Original : 2.Script/2.4.3.skinB16MC38.R
# Input    : 1.Data/ATAC_overall/
# Output   : 3.Output/2.DEseqRes/4.SkinTreg/pca_skinB16MC38.pdf
# Upstream : 10_atac_deseq_skin_treg.R
# ----------------------------------------------------------------------------

library(tidyverse)
library(DESeq2)
library(ggthemes)
library(sva)

count <- data.table::fread('1.Data/ATAC_overall/consensus_peaks.mLb.clN.featureCounts.txt')[, -c(2:6)] %>% 
  column_to_rownames('Geneid')

names(count) <- gsub('.mLb.clN.sorted.bam', '', names(count))

sampleDf <- readRDS('1.Data/ATAC_overall_sample_Infor.RDS') %>% 
  filter(project %in% c('Skin', 'MC38B16TregCD4'),
         cell == 'Treg',
         tissue %in% c('TIL', 'Skin')) %>% 
  mutate(colindi = paste0(cell, '_', tissue, 
                         case_when(
                           grepl('B16', samples) ~ '_B16',
                           grepl('MC38', samples) ~ '_MC38',
                           TRUE~''
                         )))

colSet = RColorBrewer::brewer.pal(n = 12, name = "Set3")[c(10, 11, 8)]
# colSet = c("#013e6f", "#6495ED", "#A6CEf0")
names(colSet) <- c("Treg_TIL_B16",  "Treg_TIL_MC38", "Treg_Skin")

expdata <- count[, sampleDf %>% 
                   pull(samples)]

subCondi <- sampleDf %>%
  pull(colindi) %>%
  unique()


dds <- DESeqDataSetFromMatrix(countData = expdata, 
                              colData = sampleDf, 
                              design = ~ colindi)

dds <- DESeq(dds)

# results(dds, contrast=c("tissue", 
#                         'Skin',
#                         'LN')) %>% 
#   saveRDS(paste0('3.Output/2.DEseqRes/4.SkinTreg/de_skinB16MC38.RDS'))

vsd <- vst(dds, blind=FALSE)

pcaData <- plotPCA(vsd, intgroup=c('tissue'), returnData=TRUE)

percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData %>% 
         merge(sampleDf, by.x = 'name', by.y = 'samples'), 
       aes(PC1, PC2, fill = colindi)) +
  geom_point(size = 5, shape = 21, col = 'black') +
  scale_fill_manual(breaks = subCondi,
                    values = alpha(colSet[subCondi], 0.9),
                    'Group') +
  guides(fill = guide_legend(override.aes=list(shape=21))) +
  # ggrepel::geom_text_repel(aes(label = time),
  #                          fontface = 2, size = 6, col = 'black') +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  #coord_fixed()  +
  # theme_bw() +
  theme(# legend.position = 'none',
    legend.key = element_blank(),
    aspect.ratio=1,
    text = element_text(size = 15, face = 2),
    panel.grid = element_blank(),
    axis.text = element_text(size = 15, face = 2, color = 'black'))

ggsave(paste0('3.Output/2.DEseqRes/4.SkinTreg/pca_skinB16MC38.png'),
       width = 6, height = 6, units = 'in')

ggsave(paste0('3.Output/2.DEseqRes/4.SkinTreg/pca_skinB16MC38.pdf'),
       width = 6, height = 6, units = 'in')
