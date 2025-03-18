# ----------------------------------------------------------------------------
# Fig 3A | PCA of tumour Tregs together with public non-lymphoid and resting/activated Tregs
#
# Paper    : Regulatory T cells in the tumor microenvironment display a
#            unique chromatin accessibility profile
# Original : 2.Script/2.3.exInternalCompare.R
# Input    : 1.Data/ATAC_overall/
#            public ATAC-seq from GSE112731 and GSE83315
# Output   : 3.Output/2.DEseqRes/3.GEOinduced/
# Upstream : 02_metadata_integration.R
# Note     : Public datasets must be downloaded and processed with the same pipeline first
# ----------------------------------------------------------------------------

if(!dir.exists('3.Output/2.DEseqRes/3.GEOinduced')){
  dir.create('3.Output/2.DEseqRes/3.GEOinduced')
}

library(tidyverse)
library(DESeq2)
library(ggthemes)
library(export)
library(sva)
library(factoextra)
library(ggforce)
library(cowplot)

count <- data.table::fread('1.Data/ATAC_overall/consensus_peaks.mLb.clN.featureCounts.txt')[, -c(2:6)] %>% 
  column_to_rownames('Geneid')

names(count) <- gsub('.mLb.clN.sorted.bam', '', names(count))


sampleDf <- readRDS('1.Data/ATAC_overall_sample_Infor.RDS') %>%
  filter(project %in% c('GEO', 'MC38B16TregCD4') & cell == 'Treg' |
           wkinName %in% c('BD103', 'BD105', 'BD107', 'BD109', 'BD149', 'BD145')
           ) %>% 
  mutate(TILType = ifelse(grepl("B16", samples), 'B16',
                          ifelse(grepl('MC38', samples), 'MC38', 'None')),
         tissue = ifelse(project == 'Skin' & tissue == 'LN', 'ndLN', tissue),
         tissue = ifelse(TILType %in% c('B16', 'MC38'), 
                         paste(tissue, TILType, sep = '_'),
                         tissue),
         project = ifelse(grepl('rTr|aTr', samples), 'GEO1', project))


# B16 ----
# subsmpl <- sampleDf[sampleDf$TILType != 'MC38',]

# cntNew <- ComBat_seq(as.matrix(count[, sampleDf$samples]), 
#                      batch = sampleDf$project)

dds <- DESeqDataSetFromMatrix(countData = count[, sampleDf$samples],  # count[, sampleDf$samples], cntNew,
                              colData = sampleDf,
                              design = ~ tissue)

dds <- DESeq(dds)

saveRDS(dds, '3.Output/2.DEseqRes/3.GEOinduced/exinternalCombo.RDS')

vsd <- vst(dds, blind=FALSE)

pcaData <- plotPCA(vsd, intgroup=c('tissue'), returnData=TRUE)

results <- pcaData[, 1:2]

plot_grid(fviz_nbclust(results, FUNcluster=kmeans, k.max = 10) +
            theme(aspect.ratio=1),
          fviz_nbclust(results, FUNcluster=kmeans, method="gap_stat", k.max = 10) + 
            theme_classic()+
            theme(aspect.ratio=1))

ggsave('3.Output/2.DEseqRes/3.GEOinduced/pca_clustred_thres_Kmeans.png', 
       width = 10, height = 6, units = 'in')
ggsave('3.Output/2.DEseqRes/3.GEOinduced/pca_clustred_thres_Kmeans.pdf', 
       width = 10, height = 6)

km1 <- eclust(results, "kmeans", hc_metric="eucliden",k=4)

percentVar <- round(100 * attr(pcaData, "percentVar"))

# figre with clustering 
pcaData %>% 
  merge(data.frame(clst = km1$cluster), by = 0) %>% 
  mutate(clst = paste0('C', clst)) %>% 
  ggplot(aes(PC1, PC2)) +
  geom_point(aes(fill = tissue), size = 5, col = 'black', shape = 21) +
  geom_mark_ellipse(aes(col = clst), alpha = 0.05) +
  scale_fill_brewer(palette = 'Set3') + 
  scale_color_brewer(palette = 'Dark2') +
  guides(fill = guide_legend(override.aes=list(shape=21))) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme_gray() +
  coord_fixed() +
  theme(text = element_text(size = 15, face = 2),
        axis.text = element_text(size = 15, face = 2, color = 'black'),
        legend.key = element_blank())

ggsave('3.Output/2.DEseqRes/3.GEOinduced/pca_exinternalCombo_clustred.png', 
       width = 10, height = 6, units = 'in')
ggsave('3.Output/2.DEseqRes/3.GEOinduced/pca_exinternalCombo_clustred.pdf', 
       width = 10, height = 6)

# figre with no clustering 
pcaData %>% 
  merge(data.frame(clst = km1$cluster), by = 0) %>% 
  mutate(clst = paste0('C', clst)) %>% 
  ggplot(aes(PC1, PC2)) +
  geom_point(aes(fill = tissue), size = 5, col = 'black', shape = 21) +
  # geom_mark_ellipse(aes(col = clst), alpha = 0.05) +
  scale_fill_brewer(palette = 'Set3') + 
  # scale_color_brewer(palette = 'Dark2') +
  guides(fill = guide_legend(override.aes=list(shape=21))) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme_gray() +
  coord_fixed() +
  theme(text = element_text(size = 15, face = 2),
        axis.text = element_text(size = 15, face = 2, color = 'black'),
        legend.key = element_blank())

ggsave('3.Output/2.DEseqRes/3.GEOinduced/pca_exinternalCombo_noclustred.png', 
       width = 10, height = 6, units = 'in')
ggsave('3.Output/2.DEseqRes/3.GEOinduced/pca_exinternalCombo_noclustred.pdf', 
       width = 10, height = 6)

dm<-dist(results) 
hc<-hclust(dm, method="complete") # simple dendrogram
hc$labels <- sampleDf$tissue
plot(hc, hang=-1) 
rect.hclust(hc, k=4, border="red") 
