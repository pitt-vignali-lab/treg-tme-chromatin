# ----------------------------------------------------------------------------
# Fig 3B | Differential accessibility and PCA of skin Tregs
#
# Paper    : Regulatory T cells in the tumor microenvironment display a
#            unique chromatin accessibility profile
# Original : 2.Script/2.4.1.Skin_DEseq.R
# Input    : 1.Data/ATAC_overall/
#            1.Data/ATACseq_skinB16_sample_Val.csv
# Output   : 3.Output/2.DEseqRes/4.SkinTreg/de_*.RDS
#            3.Output/2.DEseqRes/4.SkinTreg/Venn.pdf
# Upstream : 02_metadata_integration.R
# ----------------------------------------------------------------------------

if(!dir.exists('3.Output/2.DEseqRes/4.SkinTreg')){
  dir.create('3.Output/2.DEseqRes/4.SkinTreg')
}


library(tidyverse)
library(DESeq2)
library(ggthemes)

count <- data.table::fread('1.Data/ATAC_overall/consensus_peaks.mLb.clN.featureCounts.txt')[, -c(2:6)] %>% 
  column_to_rownames('Geneid')

names(count) <- gsub('.mLb.clN.sorted.bam', '', names(count))


sampleDf <- readRDS('1.Data/ATAC_overall_sample_Infor.RDS') %>% 
  filter(project == 'Skin') %>% 
  mutate(colindi = paste(cell, tissue, sep = '_'))

colSet = c("#fed589", "#FF5927", "#A6CEE3", "#013e6f")
names(colSet) <- c('CD4_LN', 'CD4_Skin', 'Treg_LN', 'Treg_Skin')

for (i in unique(sampleDf$cell)) {
  
  expdata <- count[, sampleDf %>% 
                     filter(cell == i) %>% 
                     pull(samples)]
  
  subCondi <- sampleDf %>%
    filter(cell == i) %>%
    pull(colindi) %>%
    unique()
  
  dds <- DESeqDataSetFromMatrix(countData = expdata, 
                                colData = sampleDf %>% 
                                  filter(cell == i),
                                design = ~ tissue)
  
  dds <- DESeq(dds)
  
  results(dds, contrast=c("tissue", 
                          'Skin',
                          'LN')) %>% 
    saveRDS(paste0('3.Output/2.DEseqRes/4.SkinTreg/de_', i, '.RDS'))
  
  vsd <- vst(dds, blind=FALSE)
  
  pcaData <- plotPCA(vsd, intgroup=c('tissue'), returnData=TRUE)
  
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  
  ggplot(pcaData %>% 
           merge(sampleDf, by.x = 'name', by.y = 'samples'), 
         aes(PC1, PC2, fill = colindi)) +
    geom_point(size = 5, shape = 21, col = 'black') +
    scale_fill_manual(breaks = subCondi,
                      values = alpha(colSet[subCondi], 0.8),
                      'Group') +
    guides(fill = guide_legend(override.aes=list(shape=21))) +
    ggrepel::geom_text_repel(aes(label = time),
                             fontface = 2, size = 6, col = 'black') +
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
  
  ggsave(paste0('3.Output/2.DEseqRes/4.SkinTreg/pca_', i, '.png'),
         width = 6, height = 6, units = 'in')
  
  ggsave(paste0('3.Output/2.DEseqRes/4.SkinTreg/pca_', i, '.pdf'),
         width = 6, height = 6, units = 'in')
  
}

# vennT# vpheatmap()ennT

allFile <- list.files('3.Output/2.DEseqRes/4.SkinTreg/', full.names = T)
deFile <- allFile[grepl('RDS|de', allFile)][1:2]

peaklist <- list()
for(i in deFile){
  t <- readRDS(i) %>% 
    as.data.frame() %>% 
    filter(padj < 0.05) %>% 
    rownames()
  
  peaklist[[i]] <- t
  
}

names(peaklist) <- gsub('3.Output/2.DEseqRes/4.SkinTreg/de_|.RDS', '', names(peaklist))

ggvenn::ggvenn(peaklist, 
               set_name_size = 3.5, text_size = 2.5,
               fill_color = c( "#FF5927", "#013e6f")) + coord_fixed()

ggsave('3.Output/2.DEseqRes/4.SkinTreg/Venn.png', width = 6, height = 4)
ggsave('3.Output/2.DEseqRes/4.SkinTreg/Venn.pdf', width = 6, height = 4)

ItemsList <- gplots::venn(peaklist, show.plot = FALSE)

lengths(attributes(ItemsList)$intersections)

attributes(ItemsList)$intersections$Treg

saveRDS(attributes(ItemsList)$intersections, '3.Output/2.DEseqRes/1.B16MC38TregCD4/venn_intersection.RDS')

# subsection ----
subsection <-
  readRDS('3.Output/2.DEseqRes/4.SkinTreg/de_Treg.RDS') %>% 
    as.data.frame() %>% 
    rename_with(~ paste0('Skin_', .x)) %>% 
    filter(rownames(.) %in% attributes(ItemsList)$intersections$`Treg`)

saveRDS(subsection, '3.Output/2.DEseqRes/4.SkinTreg/subsection.RDS')
