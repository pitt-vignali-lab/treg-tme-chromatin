# ----------------------------------------------------------------------------
# Fig 1B-E | Differential accessibility and PCA of Tregs, tumour vs lymph node
#
# Paper    : Regulatory T cells in the tumor microenvironment display a
#            unique chromatin accessibility profile
# Original : 2.Script/2.1.B16MC38CD4Treg_DEseq.R
# Input    : 1.Data/ATAC_overall/consensus_peaks.mLb.clN.featureCounts.txt
#            1.Data/ATAC_overall_sample_Infor.RDS
# Output   : 3.Output/2.DEseqRes/1.B16MC38TregCD4/de_*.RDS
#            3.Output/2.DEseqRes/1.B16MC38TregCD4/pca_*.pdf
# Upstream : 02_metadata_integration.R
# Note     : PCA panels are Fig 1B-C
#            replicate correlation plots are Fig 1D-E
# ----------------------------------------------------------------------------

if(!dir.exists('3.Output/2.DEseqRes')){
  dir.create('3.Output/2.DEseqRes')
}

if(!dir.exists('3.Output/2.DEseqRes/1.B16MC38TregCD4')){
  dir.create('3.Output/2.DEseqRes/1.B16MC38TregCD4')
}


library(tidyverse)
library(DESeq2)
library(ggthemes)

count <- data.table::fread('1.Data/ATAC_overall/consensus_peaks.mLb.clN.featureCounts.txt')[, -c(2:6)] %>% 
  column_to_rownames('Geneid')

names(count) <- gsub('.mLb.clN.sorted.bam', '', names(count))


sampleDf <- readRDS('1.Data/ATAC_overall_sample_Infor.RDS') %>% 
  filter(project == 'MC38B16TregCD4') %>% 
  mutate(tumor = ifelse(grepl('B16', samples), 'B16', 'MC38'),
         tumorCell = paste(tumor, cell, sep = '_'), 
        colindi = paste(cell, tissue, sep = '_'))



colSet = c("#fed589", "#FF5927", "#A6CEE3", "#013e6f")
names(colSet) <- c('CD4_LN', 'CD4_TIL', 'Treg_LN', 'Treg_TIL')

# DE and PCA
for (i in unique(sampleDf$tumorCell)) {
  
  expdata <- count[, sampleDf %>% 
                     filter(tumorCell == i) %>% 
                     pull(samples)]

  subCondi <- sampleDf %>% 
    filter(tumorCell == i) %>% 
    pull(colindi) %>% 
    unique()
  
  dds <- DESeqDataSetFromMatrix(countData = expdata, 
                                colData = sampleDf %>% 
                                  filter(tumorCell == i),
                                design = ~ tissue)
  
  dds <- DESeq(dds)
  
  results(dds, contrast=c("tissue", 
                                       'TIL',
                                       'LN')) %>% 
    saveRDS(paste0('3.Output/2.DEseqRes/1.B16MC38TregCD4/de_', i, '.RDS'))
  
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
  
  ggsave(paste0('3.Output/2.DEseqRes/1.B16MC38TregCD4/pca_', i, '.png'),
         width = 6, height = 6, units = 'in')
  
  ggsave(paste0('3.Output/2.DEseqRes/1.B16MC38TregCD4/pca_', i, '.pdf'),
         width = 6, height = 6, units = 'in')
  
}


# Correlation Coefficient graph

for (i in c('B16_Treg', 'MC38_Treg')) {
  
  expdata <- count[, sampleDf %>% 
                     filter(tumorCell == i) %>% 
                     arrange(tissue, time) %>% 
                     pull(samples)]
  
  
  
  dds <- DESeqDataSetFromMatrix(countData = expdata, 
                                colData = sampleDf %>% 
                                  filter(tumorCell == i) %>% 
                                  arrange(tissue, time),
                                design = ~ tissue)
  
  
  vsd <- vst(dds, blind=FALSE)
  
  dd <- cor(assay(vsd), method='pearson')
  
  rownames(dd) <- paste(vsd$tissue, vsd$time, sep="-")
  colnames(dd) <- paste(vsd$tissue, vsd$time, sep="-")
  
  png(paste0('3.Output/2.DEseqRes/1.B16MC38TregCD4/corrplot_', i, '.png'), 
      width = 7, height = 6, units = 'in', res = 300)
  corrplot::corrplot(dd, method="color", tl.col="black")
  dev.off()
  
  pdf(paste0('3.Output/2.DEseqRes/1.B16MC38TregCD4/corrplot_', i, '.pdf'), 
      width = 7, height = 6)
  corrplot::corrplot(dd, method="color", tl.col="black")
  dev.off()
  
}

# Pairwise comparisons of Tregs 
compDf <- data.frame(time1 = c(12, 16, 16, 20, 20, 20),
            time2 = c(8, 12, 8, 12, 8, 16))

allDe <- data.frame()

for (i in c('B16', 'MC38')) {
  for (j in 1:nrow(compDf)) {
    expdata <- count[, sampleDf %>% 
                       filter(cell == 'Treg',
                              tissue == "TIL",
                              tumor == i,
                              time %in% c(compDf$time1[j], compDf$time2[j])) %>% 
                       pull(samples)]
    
    dds <- DESeqDataSetFromMatrix(countData = expdata, 
                                  colData = sampleDf %>% 
                                    filter(cell == 'Treg',
                                           tissue == "TIL",
                                           tumor == i,
                                           time %in% c(compDf$time1[j], compDf$time2[j])) %>% 
                                    mutate(time = paste('Day ', time)),
                                  design = ~ time)
    
    dds <- DESeq(dds)
    
    de <- results(dds, contrast=c("time", 
                            paste('Day ', compDf$time1[j]),
                            paste('Day ', compDf$time2[j])))
    
    saveRDS(de, paste0('3.Output/2.DEseqRes/1.B16MC38TregCD4/Pairwise_de_', i, "_", compDf$time1[j], 'vs', compDf$time2[j], '.RDS'))
    
    allDe <- rbind(allDe,  
                   de %>%  
                           as.data.frame() %>% 
                           mutate(sample = paste0(i, '_Day', compDf$time1[j], ' vs Day', compDf$time2[j])))
    
  }
}

saveRDS(allDe, '3.Output/2.DEseqRes/1.B16MC38TregCD4/Pairwise_Treg_TIL_allpair.RDS')

for (i in c('B16', 'MC38')) {
  allDe %>%
    separate(sample, into = c('tumor', 'condi'), sep = '_') %>%
    mutate(condi = factor(condi, levels = c('Day12 vs Day8',
                                            'Day16 vs Day12',
                                            'Day16 vs Day8',
                                            'Day20 vs Day12', 
                                            'Day20 vs Day8',
                                            'Day20 vs Day16'))) %>% 
    filter(tumor == i) %>% 
    group_by(condi) %>% 
    ggplot(aes(x = log2FoldChange, y = -log10(pvalue))) +
    geom_point(aes(col = ifelse(abs(log2FoldChange) > 1 & -log10(pvalue) > 5, 'pvalue and log2FC', 
                                ifelse(abs(log2FoldChange) > 1 & -log10(pvalue) < 5, 'log2FC', 'NS'))),
               alpha = 0.5) + 
    scale_color_manual(breaks = c('NS', 'log2FC', 'pvalue and log2FC'),
                       values = c('#4d4d4d', '#1b8c42', '#ec2226')) +
    facet_wrap(~ condi, ncol = 2) +
    geom_vline(xintercept = c(-1, 1), lty = 2) + 
    geom_hline(yintercept = 5, lty = 2) +
    labs(col = '') +
    xlim(-5, 5) +
    theme_bw() +
    theme(legend.position = 'top',
          strip.background = element_blank(),
          text = element_text(face = 2, size = 15))
  
  ggsave(paste0('3.Output/2.DEseqRes/1.B16MC38TregCD4/Pairwise_de_', i, '_all.png'), 
         width = 8, height = 6, units = 'in')
  
  ggsave(paste0('3.Output/2.DEseqRes/1.B16MC38TregCD4/Pairwise_de_', i, '_all.pdf'), 
         width = 8, height = 6)
  
}


# vennT# vpheatmap()ennT
allFile <- list.files('3.Output/2.DEseqRes/1.B16MC38TregCD4/', full.names = T)
deFile <- allFile[grepl('RDS|de', allFile)][1:4]

peaklist <- list()
for(i in deFile){
  t <- readRDS(i) %>% 
    as.data.frame() %>% 
    filter(padj < 0.05) %>% 
    rownames()
  
  peaklist[[i]] <- t
  
}

names(peaklist) <- gsub('3.Output/2.DEseqRes/1.B16MC38TregCD4/de_|.RDS', '', names(peaklist))

ggvenn::ggvenn(peaklist[c(2, 4, 1, 3)], 
               set_name_size = 3.5, text_size = 2.5,
               fill_color = rep(c("#1F78B4", "#FF7F00"),
                                each = 2)) + coord_fixed()

ggsave('3.Output/2.DEseqRes/1.B16MC38TregCD4/Venn.png', width = 6, height = 4)
ggsave('3.Output/2.DEseqRes/1.B16MC38TregCD4/Venn.pdf', width = 6, height = 4)

ItemsList <- gplots::venn(peaklist, show.plot = FALSE)

lengths(attributes(ItemsList)$intersections)

attributes(ItemsList)$intersections$`B16_Treg:MC38_Treg`

saveRDS(attributes(ItemsList)$intersections, '3.Output/2.DEseqRes/1.B16MC38TregCD4/venn_intersection.RDS')

# subsection ----
subsection <- cbind(
  readRDS('3.Output/2.DEseqRes/1.B16MC38TregCD4/de_B16_Treg.RDS') %>% 
    as.data.frame() %>% 
    rename_with(~ paste0('B16_', .x)) %>% 
    filter(rownames(.) %in% attributes(ItemsList)$intersections$`B16_Treg:MC38_Treg`),
  readRDS('3.Output/2.DEseqRes/1.B16MC38TregCD4/de_MC38_Treg.RDS') %>% 
    as.data.frame() %>% 
    rename_with(~ paste0('MC38_', .x)) %>% 
    filter(rownames(.) %in% attributes(ItemsList)$intersections$`B16_Treg:MC38_Treg`))

saveRDS(subsection, '3.Output/2.DEseqRes/1.B16MC38TregCD4/subsection.RDS')

subPeak <- list(
  B16_UP = subsection %>% filter(B16_log2FoldChange > 0) %>% rownames(),
  B16_DN = subsection %>% filter(B16_log2FoldChange < 0) %>% rownames(),
  MC38_UP = subsection %>% filter(MC38_log2FoldChange > 0) %>% rownames(),
  MC38_DN = subsection %>% filter(MC38_log2FoldChange < 0) %>% rownames())

ggvenn::ggvenn(subPeak, set_name_size = 3.5, text_size = 2.5,
               fill_color = rep(c('red', "blue"), 2)) + coord_fixed()
ggsave('3.Output/2.DEseqRes/1.B16MC38TregCD4/VennSub.png', width = 6, height = 4)
ggsave('3.Output/2.DEseqRes/1.B16MC38TregCD4/VennSub.pdf', width = 6, height = 4)

# subsection 2 ----
dds <- DESeqDataSetFromMatrix(countData = count[, sampleDf$samples], 
                              colData = sampleDf,
                              design = ~ colindi)

dds <- DESeq(dds)

results(dds, contrast=c("colindi", 
                        'Treg_TIL',
                        'CD4_TIL')) %>% 
  saveRDS(paste0('3.Output/2.DEseqRes/1.B16MC38TregCD4/de_TregTILvsCD4TIL.RDS'))
