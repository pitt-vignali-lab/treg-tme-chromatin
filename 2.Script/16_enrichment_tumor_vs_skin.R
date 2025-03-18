# ----------------------------------------------------------------------------
# Fig 3J | KEGG and GO over-representation of the tumour-specific gene list
#
# Paper    : Regulatory T cells in the tumor microenvironment display a
#            unique chromatin accessibility profile
# Original : 2.Script/4.4.FuncEnrich_tumorSkin.R
# Input    : 3.Output/5.Subtraction/subsectionRes.csv
# Output   : 3.Output/4.FunctionalEnrich/4.tumorvsskin/subsection_ORA_*.csv
# Upstream : 12_subtraction_tumor_specific.R
# ----------------------------------------------------------------------------

if(!dir.exists('3.Output/4.FunctionalEnrich/4.tumorvsskin/')){
  dir.create('3.Output/4.FunctionalEnrich/4.tumorvsskin/')
}

library(clusterProfiler)
library(msigdbr)
library(fgsea)
library(tidyverse)
library(org.Mm.eg.db)

peakAnno <- readRDS('3.Output/3.PeakAnno/peakAnnotationChipseeker.RDS')
annoDf <- data.frame(peakAnno@anno)

skinDf <- readRDS('3.Output/2.DEseqRes/4.SkinTreg/subsection.RDS') %>% 
  merge(annoDf, by.x = 0, by.y = 'peakID')

tumorDf <- readRDS('3.Output/2.DEseqRes/1.B16MC38TregCD4/subsection.RDS') %>% 
  merge(annoDf, by.x = 0, by.y = 'peakID')

intesection <- intersect(skinDf$geneId,
                         tumorDf$geneId)

subG <- tumorDf %>% filter(!geneId %in% intesection) %>% 
  group_by(SYMBOL) %>% 
  mutate(Freq = n())

# kegg ----
enrichKEGG(unique(subG$geneId), 
           organism = 'mmu') %>% 
  saveRDS('3.Output/4.FunctionalEnrich/4.tumorvsskin/subsection_ORA_kegg.RDS')

kegg <- enrichKEGG(unique(subG$geneId), 
           organism = 'mmu') 
kegg@result%>% 
  write.csv('3.Output/4.FunctionalEnrich/4.tumorvsskin/subsection_ORA_kegg.csv')

ego <- enrichGO(unique(subG$geneId),
                OrgDb = 'org.Mm.eg.db',
                ont = 'BP') 

ego@result %>% 
  write.csv('3.Output/4.FunctionalEnrich/4.tumorvsskin/subsection_ORA_GO.csv')

ego3 <- gofilter(ego, level = 3)

ego3@result %>% 
  write.csv('3.Output/4.FunctionalEnrich/4.tumorvsskin/subsection_ORA_GO3.csv')

upKegg <- enrichKEGG(unique(subG[subG$B16_log2FoldChange > 0 & subG$MC38_log2FoldChange > 0,]$geneId), 
                     organism = 'mmu')

upKegg@result %>% 
  write.csv('3.Output/4.FunctionalEnrich/4.tumorvsskin/subsection_ORA_kegg_up_in_B16&MC38.csv')

dnKegg <- enrichKEGG(unique(subG[subG$B16_log2FoldChange < 0 & subG$MC38_log2FoldChange < 0,]$geneId), 
                     organism = 'mmu')
dnKegg@result %>% 
  write.csv('3.Output/4.FunctionalEnrich/4.tumorvsskin/subsection_ORA_kegg_dn_in_B16&MC38.csv')

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
# c7 gsea ----
m_df <- msigdbr(species = "Mus musculus", category = "C5", subcategory = 'BP') 

fgsea_sets <- m_df %>% 
  split(x = .$gene_symbol, f = .$gs_name)

# ranks <- subG %>% 
#   group_by(SYMBOL) %>%
#   filter(distanceToTSS == distanceToTSS[which.min(abs(distanceToTSS))]) %>% 
#   group_by(SYMBOL) %>% 
#   filter(B16_padj == min(B16_padj)) %>% 
#   arrange(desc(B16_log2FoldChange)) %>% 
#   select(SYMBOL, B16_log2FoldChange) %>% 
#   na.omit() %>% 
#   tibble::deframe()

ranks <- subG %>% 
  group_by(SYMBOL) %>%
  filter(abs(distanceToTSS) == min(abs(distanceToTSS))) %>% 
  group_by(SYMBOL) %>% 
  filter(B16_padj == min(B16_padj)) %>% 
  arrange(desc(B16_log2FoldChange)) %>% 
  select(SYMBOL, B16_log2FoldChange) %>% 
  na.omit() %>% 
  tibble::deframe()

B16Res <- fgsea(fgsea_sets, stats = ranks, minSize=10, maxSize=500, nperm=10000)
data.table::fwrite(B16Res, '3.Output/4.FunctionalEnrich/4.tumorvsskin/subsection_GSEA_C7_in_B16.csv')

ranks <- subG %>% 
  group_by(SYMBOL) %>%
  filter(distanceToTSS == distanceToTSS[which.min(abs(distanceToTSS))]) %>% 
  group_by(SYMBOL) %>% 
  filter(MC38_padj == min(MC38_padj)) %>% 
  arrange(desc(B16_log2FoldChange)) %>% 
  select(SYMBOL, B16_log2FoldChange) %>% 
  na.omit() %>% 
  tibble::deframe()

Mc38Res <- fgsea(fgsea_sets, stats = ranks)
data.table::fwrite(Mc38Res, '3.Output/4.FunctionalEnrich/4.tumorvsskin/subsection_GSEA_C7_in_MC38.csv')
