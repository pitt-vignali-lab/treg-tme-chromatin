# ----------------------------------------------------------------------------
# Fig 3H | Venn comparison of skin Treg and tumour-specific peak and gene lists
#
# Paper    : Regulatory T cells in the tumor microenvironment display a
#            unique chromatin accessibility profile
# Original : 2.Script/2.4.2.compareBwSkinTumor.R
# Input    : 3.Output/5.Subtraction/subsectionRes.csv
# Output   : 3.Output/2.DEseqRes/4.SkinTreg/VennPeaks.pdf
#            3.Output/2.DEseqRes/4.SkinTreg/VennGenes.pdf
# Upstream : 12_subtraction_tumor_specific.R
# ----------------------------------------------------------------------------

library(tidyverse)
library(openxlsx)

peakAnno <- readRDS('3.Output/3.PeakAnno/peakAnnotationChipseeker.RDS')
annoDf <- data.frame(peakAnno@anno)

skinDf <- readRDS('3.Output/2.DEseqRes/4.SkinTreg/subsection.RDS') %>% 
  merge(annoDf, by.x = 0, by.y = 'peakID')

tumorDf <- readRDS('3.Output/2.DEseqRes/1.B16MC38TregCD4/subsection.RDS') %>% 
  merge(annoDf, by.x = 0, by.y = 'peakID')

intesection <- intersect(skinDf$Row.names,
                         tumorDf$Row.names)

ggvenn::ggvenn(list(tumorSet = tumorDf$Row.names,
                    skinSet = skinDf$Row.names),
               fill_color = c("#013e6f", "#A6CEE3"),
               set_name_size = 3.5, text_size = 2.5) + coord_fixed()
ggsave('3.Output/2.DEseqRes/4.SkinTreg/VennPeaks.png', width = 6, height = 4)
ggsave('3.Output/2.DEseqRes/4.SkinTreg/VennPeaks.pdf', width = 6, height = 4)

ggvenn::ggvenn(list(tumorSet = unique(tumorDf$SYMBOL),
                    skinSet = unique(skinDf$SYMBOL)),
               fill_color = c("#013e6f", "#A6CEE3"), 
               set_name_size = 3.5, text_size = 2.5) + coord_fixed()
ggsave('3.Output/2.DEseqRes/4.SkinTreg/VennGenes.png', width = 6, height = 4)
ggsave('3.Output/2.DEseqRes/4.SkinTreg/VennGenes.pdf', width = 6, height = 4)

# export the list

interDF <- merge(tumorDf %>% 
                   group_by(SYMBOL) %>% 
                   filter(abs(distanceToTSS) == min(abs(distanceToTSS))) %>% 
                   select(grep('padj|log2FoldChange', colnames(.))),
                 skinDf %>% 
                   group_by(SYMBOL) %>% 
                   filter(abs(distanceToTSS) == min(abs(distanceToTSS))) %>% 
                   select(grep('padj|log2FoldChange', colnames(.))),
                 by = 'SYMBOL')

table(interDF$SYMBOL)[table(interDF$SYMBOL) > 1]

tumorDfSub <- tumorDf %>% 
  filter(!SYMBOL %in% interDF$SYMBOL) %>% 
  group_by(SYMBOL) %>% 
  filter(abs(distanceToTSS) == min(abs(distanceToTSS))) %>% 
  select(grep('padj|log2FoldChange', colnames(.)))

table(tumorDfSub$SYMBOL)[table(tumorDfSub$SYMBOL) > 1]

skinDfSub <- skinDf %>% 
  filter(!SYMBOL %in% interDF$SYMBOL) %>% 
  group_by(SYMBOL) %>% 
  filter(abs(distanceToTSS) == min(abs(distanceToTSS))) %>% 
  select(grep('padj|log2FoldChange', colnames(.)))

table(skinDfSub$SYMBOL)[table(skinDfSub$SYMBOL) > 1]

ItemsList <- gplots::venn(list(tumorSet = unique(tumorDf$SYMBOL),
                               skinSet = unique(skinDf$SYMBOL)), 
                          show.plot = FALSE)

names(attributes(ItemsList)$intersections)


dataset_names <- list(intersection = interDF,
                      tumorSet = tumorDfSub,
                      skinSet = skinDfSub)


write.xlsx(dataset_names, 
           file = '3.Output/2.DEseqRes/4.SkinTreg/VennGenes.xlsx')
#####
cbind(skinDf[,1:7] %>% 
  filter(Row.names %in% intesection) %>% 
  column_to_rownames('Row.names'),
  tumorDf %>% 
    filter(Row.names %in% intesection) %>% 
    column_to_rownames('Row.names')) %>% 
  write.csv('3.Output/2.DEseqRes/4.SkinTreg/skinTumorJoint.csv')

skinDf %>% 
  filter(! Row.names %in% intesection) %>% 
  column_to_rownames('Row.names') %>% 
  write.csv('3.Output/2.DEseqRes/4.SkinTreg/skinOnly.csv')

tumorDf %>% 
  filter(! Row.names %in% intesection) %>% 
  column_to_rownames('Row.names') %>% 
  write.csv('3.Output/2.DEseqRes/4.SkinTreg/tumorOnly.csv')
