# ----------------------------------------------------------------------------
# Fig 3D-G | GSEA of the tumour-specific peak list against the tissue Treg signature
#
# Paper    : Regulatory T cells in the tumor microenvironment display a
#            unique chromatin accessibility profile
# Original : 2.Script/4.3.FunctEnrich_Subtraction.R
# Input    : 3.Output/5.Subtraction/subsectionRes.csv
#            1.Data/Treg_effector_sig.xlsx
# Output   : 3.Output/4.FunctionalEnrich/3.Subtraction/*
# Upstream : 12_subtraction_tumor_specific.R
# Note     : Enrichment plots for genes up- and down-regulated in the tissue Treg signature
# ----------------------------------------------------------------------------

if(!dir.exists('3.Output/4.FunctionalEnrich/3.Subtraction/')){
  dir.create('3.Output/4.FunctionalEnrich/3.Subtraction/')
}

library(clusterProfiler)
library(msigdbr)
library(fgsea)
library(tidyverse)

subG <- read.csv('3.Output/5.Subtraction/subsectionRes.csv')
# subG2 <- readRDS('3.Output/2.DEseqRes/1.B16MC38TregCD4/de_TregTILvsCD4TIL.RDS') %>% 
#   as.data.frame() %>% 
#   filter(padj < 0.05)

# kegg ----
enrichKEGG(unique(subG$geneId), 
           organism = 'mmu') %>% 
  saveRDS('3.Output/4.FunctionalEnrich/3.Subtraction/subsection_ORA_kegg.RDS')

upKegg <- enrichKEGG(unique(subG[subG$B16_log2FoldChange > 0 & subG$MC38_log2FoldChange > 0,]$geneId), 
           organism = 'mmu')

dnKegg <- enrichKEGG(unique(subG[subG$B16_log2FoldChange < 0 & subG$MC38_log2FoldChange < 0,]$geneId), 
                     organism = 'mmu')

plotDf <- rbind(upKegg@result %>% filter(p.adjust < 0.05) %>% mutate(cat = 'up'),
                dnKegg@result %>% filter(p.adjust < 0.05) %>% mutate(cat = 'dn')) %>% 
  as.data.frame() %>% 
  mutate(cat = factor(cat, levels = c('up', 'dn')),
         Description = gsub(' - Mus musculus \\(house mouse\\)', '', Description))

plotDf %>%
  group_by(ID) %>%
  mutate(
    freq = n(),
    transp = -log(p.adjust),
    p.adjustSci = format(p.adjust, scientific = TRUE)
  ) %>% 
  filter(freq == 2) %>%
  ggplot(aes(x = cat, y = reorder(Description, -p.adjust))) +
  geom_point(aes(size = Count, col = p.adjust)) +
  labs(y = 'Pathway',
       x = '') +
  scale_color_continuous(labels = scales::scientific_format(),
                         trans = 'log') +
  theme_bw() +
  theme(text = element_text(face = 2, size = 12))
ggsave('3.Output/4.FunctionalEnrich/3.Subtraction/commonkegg.png', width = 8, height = 8, unit = 'in')
ggsave('3.Output/4.FunctionalEnrich/3.Subtraction/commonkegg.pdf', width = 8, height = 8)

plotDf %>% 
  group_by(ID) %>% 
  mutate(freq = n()) %>% 
  filter(freq == 1) %>% 
  mutate(ord = ifelse(cat == 'up', 1, -1)) %>% 
  ggplot(aes(x = cat, y = reorder(Description, ord))) +
  geom_point(aes(size = Count, col = p.adjust)) +
  labs(y = 'Pathway',
       x = '') +
  scale_color_continuous(labels = scales::scientific_format(),
                         trans = 'log') +
  theme_bw() +
  theme(text = element_text(face = 2, size = 12))
ggsave('3.Output/4.FunctionalEnrich/3.Subtraction/uniqkegg.png', width = 8, height = 12, unit = 'in')
ggsave('3.Output/4.FunctionalEnrich/3.Subtraction/uniqkegg.pdf', width = 8, height = 12)

# C7 ----

m_df <- msigdbr(species = "Mus musculus", category = "C7") 

fgsea_sets <- m_df %>% 
  split(x = .$gene_symbol, f = .$gs_name)

ranks_b16 <- subG %>% 
  group_by(SYMBOL) %>% 
  filter(distanceToTSS == distanceToTSS[which.min(abs(distanceToTSS))]) %>% 
  group_by(SYMBOL) %>% 
  filter(B16_padj == min(B16_padj)) %>% 
  arrange(desc(B16_log2FoldChange)) %>% 
  select(SYMBOL, B16_log2FoldChange) %>% 
  na.omit() %>% 
  tibble::deframe() 

fgseaRes_B16 <- fgsea(fgsea_sets, stats = ranks_b16)

data.table::fwrite(fgseaRes_B16, '3.Output/4.FunctionalEnrich/3.Subtraction/subtractionB16C7.csv')

ranks_mc38 <- subG %>% 
  group_by(SYMBOL) %>% 
  filter(distanceToTSS == distanceToTSS[which.min(abs(distanceToTSS))]) %>% 
  group_by(SYMBOL) %>% 
  filter(MC38_padj == min(MC38_padj)) %>% 
  arrange(desc(MC38_log2FoldChange)) %>% 
  select(SYMBOL, MC38_log2FoldChange) %>% 
  na.omit() %>% 
  tibble::deframe() 

fgseaRes_MC38 <- fgsea(fgsea_sets, stats = ranks_mc38)

data.table::fwrite(fgseaRes_MC38, '3.Output/4.FunctionalEnrich/3.Subtraction/subtractionMC38C7.csv')
