# ----------------------------------------------------------------------------
# Fig 2A, 2C-J | Peak counts, GSEA against Treg signatures and KEGG enrichment
#
# Paper    : Regulatory T cells in the tumor microenvironment display a
#            unique chromatin accessibility profile
# Original : 2.Script/4.1.FuncEnrich_B16MC38TregCD4.R
# Input    : 3.Output/3.PeakAnno/de_anno_Treg_*.csv
#            1.Data/Treg_effector_sig.xlsx
# Output   : 3.Output/4.FunctionalEnrich/1.B16MC38TregCD4/*
# Upstream : 04_atac_peak_annotation.R
# Note     : GSEA uses MSigDB C7 (GSE7852_TREG_VS_TCONV_LN_DN) and the Treg transcriptional signature
# ----------------------------------------------------------------------------

if(!dir.exists('3.Output/4.FunctionalEnrich')){
  dir.create('3.Output/4.FunctionalEnrich')
}

if(!dir.exists('3.Output/4.FunctionalEnrich/1.B16MC38TregCD4')){
  dir.create('3.Output/4.FunctionalEnrich/1.B16MC38TregCD4')
}

library(clusterProfiler)
library(msigdbr)
library(fgsea)
library(tidyverse)
library(readxl)

# ora KEGG ----
# B16
enrichKEGG(read.csv('3.Output/3.PeakAnno/de_anno_Treg_B16_(p0.05)up.csv')$geneId, 
           organism = 'mmu') %>% 
  saveRDS('3.Output/4.FunctionalEnrich/1.B16MC38TregCD4/Treg_B16_activated(p0.05).RDS')

enrichKEGG(read.csv('3.Output/3.PeakAnno/de_anno_Treg_B16_(p0.05)down.csv')$geneId, 
           organism = 'mmu') %>% 
  saveRDS('3.Output/4.FunctionalEnrich/1.B16MC38TregCD4/Treg_B16_suppressed(p0.05).RDS')

enrichKEGG(read.csv('3.Output/3.PeakAnno/de_anno_Treg_B16_(padj0.05)up.csv')$geneId, 
           organism = 'mmu') %>% 
  saveRDS('3.Output/4.FunctionalEnrich/1.B16MC38TregCD4/Treg_B16_activated(padj0.05).RDS')

enrichKEGG(read.csv('3.Output/3.PeakAnno/de_anno_Treg_B16_(padj0.05)down.csv')$geneId, 
           organism = 'mmu') %>% 
  saveRDS('3.Output/4.FunctionalEnrich/1.B16MC38TregCD4/Treg_B16_suppressed(padj0.05).RDS')

# MC38
enrichKEGG(read.csv('3.Output/3.PeakAnno/de_anno_Treg_MC38_(p0.05)up.csv')$geneId, 
           organism = 'mmu') %>% 
  saveRDS('3.Output/4.FunctionalEnrich/1.B16MC38TregCD4/Treg_MC38_activated(p0.05).RDS')

enrichKEGG(read.csv('3.Output/3.PeakAnno/de_anno_Treg_MC38_(p0.05)down.csv')$geneId, 
           organism = 'mmu') %>% 
  saveRDS('3.Output/4.FunctionalEnrich/1.B16MC38TregCD4/Treg_MC38_suppressed(p0.05).RDS')

enrichKEGG(read.csv('3.Output/3.PeakAnno/de_anno_Treg_MC38_(padj0.05)up.csv')$geneId, 
           organism = 'mmu') %>% 
  saveRDS('3.Output/4.FunctionalEnrich/1.B16MC38TregCD4/Treg_MC38_activated(padj0.05).RDS')

enrichKEGG(read.csv('3.Output/3.PeakAnno/de_anno_Treg_MC38_(padj0.05)down.csv')$geneId, 
           organism = 'mmu') %>% 
  saveRDS('3.Output/4.FunctionalEnrich/1.B16MC38TregCD4/Treg_MC38_suppressed(padj0.05).RDS')

# plot
for (i in c('B16', 'MC38')) {
  for (j in c('p0.05', 'padj0.05')) {
   
    actPath <- readRDS(paste0('3.Output/4.FunctionalEnrich/1.B16MC38TregCD4/Treg_', i, '_activated(', j, ').RDS'))@result %>% 
      mutate(Description = gsub(' - Mus musculus (house mouse)', '', Description, fixed = T)) %>% 
      filter(p.adjust < 0.05) %>% 
      mutate(pathCat = 'Activated')
    
    supPath <- readRDS(paste0('3.Output/4.FunctionalEnrich/1.B16MC38TregCD4/Treg_', i, '_suppressed(', j, ').RDS'))@result %>% 
      mutate(Description = gsub(' - Mus musculus (house mouse)', '', Description, fixed = T)) %>% 
      filter(p.adjust < 0.05) %>% 
      mutate(pathCat = 'Suppressed')
    
    int <- intersect(actPath$ID, supPath$ID)
    
    
    plotDf <- rbind(actPath %>% filter(!ID %in% int),
                    supPath %>% filter(!ID %in% int)) %>% 
      mutate(GeneRatio = sapply(GeneRatio, function(x) eval(parse(text=x))))
    
    ggplot(plotDf, 
           aes(x = GeneRatio, 
               y = reorder(Description, GeneRatio))) +
      geom_point(aes(size = Count, fill = p.adjust), shape = 21) + 
      facet_grid(~pathCat) +
      scale_fill_gradient(high = 'blue', low = 'red') +
      ggtitle("Enriched KEGG TIL Treg vs LN Treg") +
      theme_bw()
    ggsave(paste0('3.Output/4.FunctionalEnrich/1.B16MC38TregCD4/', i, '_KEGG(', j, ').png'), width = 8, height = 12, units = 'in')
    ggsave(paste0('3.Output/4.FunctionalEnrich/1.B16MC38TregCD4/', i, '_KEGG(', j, ').pdf'), width = 8, height = 12)
     
  }
}




# genes <- list(Activated = read.csv('3.Output/3.PeakAnno/de_anno_Treg_B16_(p0.05)up.csv')$geneId,
#               suppressed = read.csv('3.Output/3.PeakAnno/de_anno_Treg_B16_(p0.05)down.csv')$geneId)
# 
# 
# compKEGG <- compareCluster(geneCluster   = genes,
#                            fun           = "enrichKEGG",
#                            organism = 'mmu',
#                            pvalueCutoff  = 0.05,
#                            pAdjustMethod = "BH")
# 
# t <- compKEGG %>% 
#   group_by(ID) %>% 
#   summarise(n = n()) %>% 
#   filter(n == 1) %>% 
#   pull(ID)
# 
# dotplot(compKEGG %>% 
#           mutate(Description = gsub(' - Mus musculus (house mouse)', '', Description, fixed = T)) %>% 
#           filter(ID %in% t), 
#         showCategory = 30, 
#         title = "Enriched KEGG TIL Treg vs LN Treg")

# GSEA C7 ----

m_df <- msigdbr(species = "Mus musculus", category = "C7") 

fgsea_sets <- m_df %>% 
  split(x = .$gene_symbol, f = .$gs_name)

peakAnno <- readRDS('3.Output/3.PeakAnno/peakAnnotationChipseeker.RDS')
annoDf <- data.frame(peakAnno@anno)

# If a gene is associated with multiple peaks, the one that is closest to the gene TSS is used in these analyses. 
b16Ranks <- readRDS('3.Output/2.DEseqRes/1.B16MC38TregCD4/de_B16_Treg.RDS') %>% 
  as.data.frame() %>% 
  merge(annoDf, by.x = 0, by.y = 'peakID') %>% 
  group_by(SYMBOL) %>%
  filter(distanceToTSS == distanceToTSS[which.min(abs(distanceToTSS))]) %>% 
  group_by(SYMBOL) %>% 
  filter(padj == min(padj)) %>% 
  arrange(desc(log2FoldChange)) %>% 
  select(SYMBOL, log2FoldChange) %>% 
  na.omit() %>% 
  tibble::deframe()

fgseaRes <- fgsea(fgsea_sets, stats = b16Ranks)

saveRDS(fgseaRes, '3.Output/4.FunctionalEnrich/1.B16MC38TregCD4/gsea_Treg_B16_TILvsLN_C7.RDS')
data.table::fwrite(fgseaRes, '3.Output/4.FunctionalEnrich/1.B16MC38TregCD4/gsea_Treg_B16_TILvsLN_C7.csv')

selectedPtw <- c('GSE13522_WT_VS_IFNAR_KO_SKIN_DN',
                 'GSE14415_INDUCED_TREG_VS_FOXP3_KO_INDUCED_TREG_DN',
                 'GSE14415_FOXP3_KO_NATURAL_TREG_VS_TCONV_UP',
                 'GSE40274_CTRL_VS_FOXP3_TRANSDUCED_ACTIVATED_CD4_TCELL_DN',
                 'GSE7852_TREG_VS_TCONV_DN',
                 'GSE13522_WT_VS_IFNAR_KO_SKIN_UP',
                 'GSE7460_FOXP3_MUT_VS_WT_ACT_TCONV_DN',
                 'GSE14415_INDUCED_TREG_VS_FOXP3_KO_INDUCED_TREG_UP',
                 'GSE21678_WT_VS_FOXO1_FOXO3_KO_TREG_UP',
                 'GSE7460_TCONV_VS_TREG_LN_UP',
                 'GSE14415_FOXP3_KO_NATURAL_TREG_VS_TCONV_DN',
                 'GSE40274_LEF1_VS_FOXP3_AND_LEF1_TRANSDUCED_ACTIVATED_CD4_TCELL_DN',
                 'GSE7852_TREG_VS_TCONV_LN_DN')

fgseaRes %>%  
  filter(pathway %in% selectedPtw) %>% 
  ggplot(aes(x = NES, y = reorder(pathway, NES), fill = padj)) +
  geom_bar(stat = 'identity') +
  theme_minimal()
ggsave('3.Output/4.FunctionalEnrich/1.B16MC38TregCD4/gsea_selected.png', 
       width = 15, height = 5, units = 'in', bg = 'white')
ggsave('3.Output/4.FunctionalEnrich/1.B16MC38TregCD4/gsea_selected.pdf', 
       width = 15, height = 5, bg = 'white')


fgseaRes %>%  
  filter(padj < 0.05) %>% 
  arrange(NES) %>% 
  group_by(NES > 0) %>% 
  top_n(n = 6, wt = abs(NES)) %>% 
  ggplot(aes(x = NES, y = reorder(pathway, NES), fill = padj)) +
  geom_bar(stat = 'identity') +
  theme_minimal()
ggsave('3.Output/4.FunctionalEnrich/1.B16MC38TregCD4/gsea_unbiased.png', 
       width = 15, height = 5, units = 'in', bg = 'white')
ggsave('3.Output/4.FunctionalEnrich/1.B16MC38TregCD4/gsea_unbiased.pdf', 
       width = 15, height = 5, bg = 'white')

nes = fgseaRes[fgseaRes$pathway == "GSE7852_TREG_VS_TCONV_LN_DN", ]$NES
padj = fgseaRes[fgseaRes$pathway == "GSE7852_TREG_VS_TCONV_LN_DN", ]$padj

plotEnrichment(fgsea_sets[["GSE7852_TREG_VS_TCONV_LN_DN"]],
               b16Ranks) + 
  labs(title="GSE7852_TREG_VS_TCONV_LN_DN_B16",
       subtitle = paste("NES: ", round(nes, 3), ' p.adjust: ', format(padj, scientific = T, digit = 3)))
ggsave('3.Output/4.FunctionalEnrich/1.B16MC38TregCD4/gsea_curve_GSE7852_TREG_VS_TCONV_LN_DN_B16.png',
       width = 6, height = 4, units = 'in')
ggsave('3.Output/4.FunctionalEnrich/1.B16MC38TregCD4/gsea_curve_GSE7852_TREG_VS_TCONV_LN_DN_B16.pdf',
       width = 6, height = 4)

# mc38 small scale 
mc38Ranks <- readRDS('3.Output/2.DEseqRes/1.B16MC38TregCD4/de_MC38_Treg.RDS') %>% 
  as.data.frame() %>% 
  merge(annoDf, by.x = 0, by.y = 'peakID') %>% 
  group_by(SYMBOL) %>% 
  filter(distanceToTSS == distanceToTSS[which.min(abs(distanceToTSS))]) %>% 
  group_by(SYMBOL) %>% 
  filter(padj == min(padj)) %>% 
  arrange(desc(log2FoldChange)) %>% 
  select(SYMBOL, log2FoldChange) %>% 
  na.omit() %>% 
  tibble::deframe()

tempGsea <- fgsea(fgsea_sets["GSE7852_TREG_VS_TCONV_LN_DN"], stats = mc38Ranks)

nes = tempGsea[tempGsea$pathway == "GSE7852_TREG_VS_TCONV_LN_DN", ]$NES
padj = tempGsea[tempGsea$pathway == "GSE7852_TREG_VS_TCONV_LN_DN", ]$padj

plotEnrichment(fgsea_sets[["GSE7852_TREG_VS_TCONV_LN_DN"]],
               mc38Ranks) + 
  labs(title="GSE7852_TREG_VS_TCONV_LN_DN_MC38",
       subtitle = paste("NES: ", round(nes, 3), ' p.adjust: ', format(padj, scientific = T, digit = 3)))

ggsave('3.Output/4.FunctionalEnrich/1.B16MC38TregCD4/gsea_curve_GSE7852_TREG_VS_TCONV_LN_DN_MC38.png',
       width = 6, height = 4, units = 'in')
ggsave('3.Output/4.FunctionalEnrich/1.B16MC38TregCD4/gsea_curve_GSE7852_TREG_VS_TCONV_LN_DN_MC38.pdf',
       width = 6, height = 4)

# treg signature -----------------
tsig <- read_xlsx('1.Data/Treg_effector_sig.xlsx')[-c(1:6), 1:7] %>% 
  janitor::row_to_names(1) %>% 
  as.data.frame() %>% 
  mutate(adj.P.Val = as.numeric(adj.P.Val),
         logFC = as.numeric(logFC)) %>% 
  filter(adj.P.Val < 0.05)


for (cell in c('B16', 'MC38')) {
  for (stat in c('pos', 'neg')){
    
    ifelse(stat == 'pos', 
           tsigSub <-  tsig[tsig$logFC > 0,]$Symbol,
           tsigSub <- tsig[tsig$logFC < 0,]$Symbol)
    
    ifelse(cell == 'B16', cellRanks <-  b16Ranks, cellRanks <- mc38Ranks)
    
    tempGsea <- fgsea(list(tsig = tsigSub), stats = cellRanks)
    
    nes = tempGsea[tempGsea$pathway == "tsig", ]$NES
    padj = tempGsea[tempGsea$pathway == "tsig", ]$padj
    
    plotEnrichment(tsigSub,
                   cellRanks) + 
      labs(title="Treg transcritional signature (logFC < 0)",
           subtitle = paste(cell, " NES: ", round(nes, 3), ' p.adjust: ', format(padj, scientific = T, digit = 3)))
    
    ggsave(paste0('3.Output/4.FunctionalEnrich/1.B16MC38TregCD4/gsea_curve_Treg_transcritional_signature(logfc', stat, ')_', cell ,'.png'),
           width = 6, height = 4, units = 'in')
    ggsave(paste0('3.Output/4.FunctionalEnrich/1.B16MC38TregCD4/gsea_curve_Treg_transcritional_signature(logfc', stat, ')_', cell ,'.pdf'),
           width = 6, height = 4)
    
  }
}

# consistent tissue peak treg signature ----------------------
peakAnno <- readRDS('3.Output/3.PeakAnno/peakAnnotationChipseeker.RDS')
annoDf <- data.frame(peakAnno@anno)

# subsection anno 
subTreg <- readRDS('3.Output/2.DEseqRes/1.B16MC38TregCD4/subsection.RDS') %>% 
  merge(annoDf, by.x = 0, by.y = 'peakID') 

for (cell in c('B16', 'MC38')) {
  for (stat in c('pos', 'neg')){
    
    ifelse(stat == 'pos', 
           tsigSub <-  tsig[tsig$logFC > 0,]$Symbol,
           tsigSub <- tsig[tsig$logFC < 0,]$Symbol)
    
    newRanks <- subTreg %>% 
      rename_with(~ gsub(paste0(cell, '_'), '', .x)) %>% 
      group_by(SYMBOL) %>% 
      filter(distanceToTSS == distanceToTSS[which.min(abs(distanceToTSS))]) %>% 
      group_by(SYMBOL) %>% 
      filter(padj == min(padj)) %>% 
      select(SYMBOL, log2FoldChange) %>% 
      na.omit() %>% 
      tibble::deframe()
    
    
    tempGsea <- fgsea(list(tsig = tsigSub), stats = newRanks, nPermSimple = 50000)
    
    nes = tempGsea[tempGsea$pathway == "tsig", ]$NES
    padj = tempGsea[tempGsea$pathway == "tsig", ]$padj
    
    plotEnrichment(tsigSub,
                   newRanks) + 
      labs(title="Tissue Treg signature (logFC < 0)",
           subtitle = paste(cell, " NES: ", round(nes, 3), ' p.adjust: ', format(padj, scientific = T, digit = 3)))
    
    ggsave(paste0('3.Output/4.FunctionalEnrich/1.B16MC38TregCD4/gsea_curve_Tissue_Treg_signature(logfc', stat, ')_', cell ,'.png'),
           width = 6, height = 4, units = 'in')
    ggsave(paste0('3.Output/4.FunctionalEnrich/1.B16MC38TregCD4/gsea_curve_Tissue_Treg_signature(logfc', stat, ')_', cell ,'.pdf'),
           width = 6, height = 4)
    
  }
}

###B16
# allTregB16 <- allTreg %>%  
#   select(SYMBOL, B16_log2FoldChange) %>% 
#   na.omit() %>% 
#   tibble::deframe()
# 
# tempGsea <- fgsea(list(tsig = tsig[tsig$logFC < 0,]$Symbol), stats = allTregB16)
# 
# nes = tempGsea[tempGsea$pathway == "tsig", ]$NES
# padj = tempGsea[tempGsea$pathway == "tsig", ]$padj
# 
# plotEnrichment(tsig[tsig$logFC < 0,]$Symbol,
#                allTregB16) + 
#   labs(title="Tissue Treg signature (logFC < 0)",
#        subtitle = paste("B16 NES: ", round(nes, 3), ' p.adjust: ', format(padj, scientific = T, digit = 3)))
# 
# ggsave('3.Output/4.FunctionalEnrich/1.B16MC38TregCD4/gsea_curve_Tissue_Treg_signature(logfcNeg)_B16.png',
#        width = 6, height = 4, units = 'in')
# 
# ###MC38
# allTregMC38 <- allTreg %>%  
#   select(SYMBOL, MC38_log2FoldChange) %>% 
#   na.omit() %>% 
#   tibble::deframe()
# 
# tempGsea <- fgsea(list(tsig = tsig[tsig$logFC > 0,]$Symbol), stats = allTregMC38)
# 
# nes = tempGsea[tempGsea$pathway == "tsig", ]$NES
# padj = tempGsea[tempGsea$pathway == "tsig", ]$padj
# 
# plotEnrichment(tsig[tsig$logFC > 0,]$Symbol,
#                allTregB16) + 
#   labs(title="Tissue Treg signature (logFC > 0)",
#        subtitle = paste("MC38 NES: ", round(nes, 3), ' p.adjust: ', format(padj, scientific = T, digit = 3)))
# 
# ggsave('3.Output/4.FunctionalEnrich/1.B16MC38TregCD4/gsea_curve_Tissue_Treg_signature(logfcPos)_MC38.png',
#        width = 6, height = 4, units = 'in')
