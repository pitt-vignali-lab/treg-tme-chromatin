# ----------------------------------------------------------------------------
# Fig 3 | KEGG over-representation for skin-only, tumour-only and shared genes
#
# Paper    : Regulatory T cells in the tumor microenvironment display a
#            unique chromatin accessibility profile
# Original : 2.Script/4.2.FuncEnrich_SkinTumorTreg.R
# Input    : 3.Output/5.Subtraction/
#            3.Output/3.PeakAnno/
# Output   : 3.Output/4.FunctionalEnrich/skinTumorTreg_ORA_*.RDS
# Upstream : 13_compare_skin_vs_tumor.R
# ----------------------------------------------------------------------------

if(!dir.exists('3.Output/4.FunctionalEnrich/2.SkinTumorTreg/')){
  dir.create('3.Output/4.FunctionalEnrich/2.SkinTumorTreg/')
}

library(clusterProfiler)
library(msigdbr)
library(fgsea)
library(tidyverse)
library(readxl)

peakAnno <- readRDS('3.Output/3.PeakAnno/peakAnnotationChipseeker.RDS')
annoDf <- data.frame(peakAnno@anno)

skinDf <- readRDS('3.Output/2.DEseqRes/4.SkinTreg/subsection.RDS') %>% 
  merge(annoDf, by.x = 0, by.y = 'peakID')

tumorDf <- readRDS('3.Output/2.DEseqRes/1.B16MC38TregCD4/subsection.RDS') %>% 
  merge(annoDf, by.x = 0, by.y = 'peakID')


interG <- intersect(tumorDf$geneId, skinDf$geneId)

enrichKEGG(interG, 
           organism = 'mmu') %>% 
  saveRDS('3.Output/4.FunctionalEnrich/skinTumorTreg_ORA_intersection.RDS')

enrichKEGG(skinDf %>% filter(! geneId %in% interG) %>% pull(geneId), 
           organism = 'mmu') %>% 
  saveRDS('3.Output/4.FunctionalEnrich/skinTumorTreg_ORA_skinOnly.RDS')

enrichKEGG(tumorDf %>% filter(! geneId %in% interG) %>% pull(geneId), 
           organism = 'mmu') %>% 
  saveRDS('3.Output/4.FunctionalEnrich/skinTumorTreg_ORA_tumorOnly.RDS')

interFuncDf <- enrichKEGG(interG, 
                      organism = 'mmu') %>% 
  as.data.frame() %>% 
  filter(p.adjust < 0.05)

skinFuncDf <- enrichKEGG(skinDf %>% filter(! geneId %in% interG) %>% pull(geneId), 
                     organism = 'mmu') %>% 
  as.data.frame() %>% 
  filter(p.adjust < 0.05)

tumorFuncDf <- enrichKEGG(tumorDf %>% filter(! geneId %in% interG) %>% pull(geneId), 
                          organism = 'mmu') %>% 
  as.data.frame() %>% 
  filter(p.adjust < 0.05)
