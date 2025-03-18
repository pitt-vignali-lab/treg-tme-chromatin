# ----------------------------------------------------------------------------
# Fig 3I-J | GO bar plot and fold-change visualisation of the tumour-specific list
#
# Paper    : Regulatory T cells in the tumor microenvironment display a
#            unique chromatin accessibility profile
# Original : 2.Script/4.4.1.egoVis.R
# Input    : 3.Output/4.FunctionalEnrich/4.tumorvsskin/subsection_ORA_GO.csv
# Output   : 3.Output/4.FunctionalEnrich/4.tumorvsskin/subsection_GO_Bar.pdf
#            3.Output/4.FunctionalEnrich/4.tumorvsskin/subsection_DEG.pdf
# Upstream : 16_enrichment_tumor_vs_skin.R
# ----------------------------------------------------------------------------

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

# ego <- enrichGO(unique(subG$geneId),
#                 OrgDb = 'org.Mm.eg.db',
#                 ont = 'BP', 
#                 readable = T) 

# ego@result %>% 
#   write.csv('3.Output/4.FunctionalEnrich/4.tumorvsskin/subsection_ORA_GO_readable.csv')

ego <- read.csv('3.Output/4.FunctionalEnrich/4.tumorvsskin/subsection_ORA_GO_readable.csv')

ego <- ego[c(1:12, 16, 18:21, 28, 29, 39, 45, 55), ]


ggplot(ego, aes(x = Count, y = reorder(str_wrap(str_to_sentence(Description), width = 70), -p.adjust), fill = p.adjust)) +
  geom_bar(stat = 'identity') +
  scale_fill_continuous(trans = 'reverse', labels = scales::scientific_format(digits = 4)) +
  theme_minimal() +
  labs(x = "Gene count",
    y = "Description",
    fill = "Adjusted P-value"
  ) +
  theme(text = element_text(face = 2))

ggsave('3.Output/4.FunctionalEnrich/4.tumorvsskin/subsection_GO_Bar.png', 
       width = 12, height = 5, units = 'in')
ggsave('3.Output/4.FunctionalEnrich/4.tumorvsskin/subsection_GO_Bar.pdf', 
       width = 12, height = 5, units = 'in')


volDf <- subG %>% 
  filter(abs(distanceToTSS) == min(abs(distanceToTSS)))

labDf <- unique(rbind(volDf %>%
                        filter(B16_padj < 0.05) %>% 
                        group_by(B16_log2FoldChange < 0) %>% 
                        top_n(10, -B16_padj),
                      volDf %>% filter(SYMBOL %in% strsplit(ego$geneID[1], '/')[[1]])))

b16 <- ggplot(volDf, aes(x = B16_log2FoldChange, y = -log(B16_padj))) +
  geom_point(aes(col = B16_log2FoldChange < 0,
                 alpha = ifelse(SYMBOL %in% labDf$SYMBOL, 0.8, 0.5))) +
  scale_color_manual(values = c('red', 'blue')) +
  geom_text(data =  labDf,
                           aes(label = SYMBOL),
                           fontface = 2) +
  labs(x = 'log2 fold change',
       y = '-log(pavelu)',
       title = 'B16') +
  theme_bw() +
  theme(aspect.ratio = 1,
        legend.position = 'none')


mc38 <- ggplot(volDf, aes(x = MC38_log2FoldChange, y = -log(MC38_padj))) +
  geom_point(aes(col = MC38_log2FoldChange < 0,
                 alpha = ifelse(SYMBOL %in% labDf$SYMBOL, 0.8, 0.5))) +
  scale_color_manual(values = c('red', 'blue')) +
  geom_text(data =  labDf,
            aes(label = SYMBOL),
            fontface = 2) +
  labs(x = 'log2 fold change',
       y = '-log(pavelu)',
       title = 'MC38') +
  theme_bw() +
  theme(aspect.ratio = 1,
        legend.position = 'none')

b16 + mc38

ggsave('3.Output/4.FunctionalEnrich/4.tumorvsskin/subsection_DEG.png', 
       width = 8, height = 5, units = 'in')
ggsave('3.Output/4.FunctionalEnrich/4.tumorvsskin/subsection_DEG.pdf', 
       width = 8, height = 5, units = 'in')


b16 <- ggplot(volDf, aes(x = B16_log2FoldChange, y = -log(B16_padj))) +
  geom_point(aes(col = B16_log2FoldChange < 0), alpha = 0.5) +
  scale_color_manual(values = c('red', 'blue')) +
  ggrepel::geom_text_repel(data = unique(volDf %>%
                                           filter(B16_padj < 0.05) %>% 
                                           group_by(B16_log2FoldChange < 0) %>% 
                                           top_n(10, -B16_padj)),
                           aes(label = SYMBOL)) +
  labs(x = 'log2 fold change',
       y = '-log(pavelu)',
       title = 'B16') +
  theme_bw() +
  theme(aspect.ratio = 1,
        legend.position = 'none')


mc38 <- ggplot(volDf, aes(x = MC38_log2FoldChange, y = -log(MC38_padj))) +
  geom_point(aes(col = MC38_log2FoldChange < 0), alpha = 0.5) +
  scale_color_manual(values = c('red', 'blue')) +
  ggrepel::geom_text_repel(data = unique(volDf %>%
                                           filter(MC38_padj < 0.05) %>% 
                                           group_by(MC38_log2FoldChange < 0) %>% 
                                           top_n(10, -MC38_padj)),
                           aes(label = SYMBOL)) +
  labs(x = 'log2 fold change',
       y = '-log(pavelu)',
       title = 'MC38') +
  theme_bw() +
  theme(aspect.ratio = 1,
        legend.position = 'none')

b16 + mc38
ggsave('3.Output/4.FunctionalEnrich/4.tumorvsskin/subsection_DEG_top20.png', 
       width = 8, height = 4, units = 'in')
ggsave('3.Output/4.FunctionalEnrich/4.tumorvsskin/subsection_DEG_top20.pdf', 
       width = 8, height = 4, units = 'in')




subGUni <- subG %>% 
  group_by(SYMBOL) %>% 
  filter(abs(distanceToTSS) == min(abs(distanceToTSS))) %>% 
  summarise(logFC = mean(mean(B16_log2FoldChange), mean(MC38_log2FoldChange)))

gLst <- strsplit(ego$geneID, split = '/')
names(gLst) <- ego$Description
gDf <- data.frame(unlist(gLst)) %>% 
  dplyr::rename('SYMBOL' = 1) %>% 
  rownames_to_column('ptw') %>% 
  mutate(ptw = gsub('[0-9]+$', '', ptw)) %>% 
  merge(ego, by.x = 'ptw', by.y = 'Description') %>% 
  mutate(ptw = factor(ptw, 
                      levels = ego %>% arrange(p.adjust) %>% pull(Description))) %>% 
  merge(subGUni, by = 'SYMBOL') %>% 
  group_by(ptw) %>% 
  arrange(logFC) %>% 
  mutate(order = 1:n())

ggplot(gDf %>% 
         mutate(ptw = str_wrap(ptw, width = 55)), 
       aes(x = order, y = logFC)) +
  geom_bar(stat = 'identity', width = 0.3, aes(fill = logFC < 0), alpha = 0.5) +
  geom_line(aes(col = logFC < 0)) +
  facet_grid(ptw~., switch="y") +
  theme_void() +
  theme(legend.position = 'bottom',
        strip.text.y.left = element_text(angle = 0, hjust = 1)) +
  coord_fixed(ratio = 0.5)
