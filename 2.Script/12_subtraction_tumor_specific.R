# ----------------------------------------------------------------------------
# Fig 3C, 3H | Subtraction method isolating peaks unique to tumour Tregs
#
# Paper    : Regulatory T cells in the tumor microenvironment display a
#            unique chromatin accessibility profile
# Original : 2.Script/5.Subtraction.R
# Input    : 3.Output/2.DEseqRes/
# Output   : 3.Output/5.Subtraction/subsectionRes.csv
#            3.Output/5.Subtraction/Venn.png
# Upstream : 10_atac_deseq_skin_treg.R
#            09_atac_deseq_public_datasets.R
# Note     : Produces the tumour-specific peak list used by Fig 3D-J
# ----------------------------------------------------------------------------

if(!dir.exists('3.Output/5.Subtraction')){
  dir.create('3.Output/5.Subtraction')
}

library(tidyverse)


# method 1, DE ----

peakAnno <- readRDS('3.Output/3.PeakAnno/peakAnnotationChipseeker.RDS')
annoDf <- data.frame(peakAnno@anno)

allFile <- list.files('3.Output/2.DEseqRes/', full.names = T) 

peaklist <- list()
for(i in allFile[grep('jc_de', allFile)]){
  t <- readRDS(i) %>% 
    as.data.frame() %>% 
    filter(padj < 0.05) %>% 
    rownames()
  
  peaklist[[i]] <- t
  
}

names(peaklist) <- gsub('3.Output/2.DEseqRes/jc_de_|.RDS', '', names(peaklist))

ggvenn::ggvenn(peaklist[c(2, 4, 1, 3)], set_name_size = 3.5, text_size = 2.5) + coord_fixed()
ggsave('3.Output/5.Subtraction/Venn.png', width = 6, height = 4)

ItemsList <- gplots::venn(peaklist, show.plot = FALSE)

lengths(attributes(ItemsList)$intersections)

attributes(ItemsList)$intersections$`B16_Treg:MC38_Treg`

saveRDS(attributes(ItemsList)$intersections, '3.Output/5.Subtraction/venn_intersection.RDS')

b16Treg <- readRDS('3.Output/2.DEseqRes/jc_de_B16_Treg.RDS') %>% 
  as.data.frame() %>% 
  filter(rownames(.) %in% attributes(ItemsList)$intersections$`B16_Treg:MC38_Treg`)

names(b16Treg) <- paste0('B16_', names(b16Treg))
  
mc38Treg <- readRDS('3.Output/2.DEseqRes/jc_de_MC38_Treg.RDS') %>% 
  as.data.frame() %>% 
  filter(rownames(.) %in% attributes(ItemsList)$intersections$`B16_Treg:MC38_Treg`)

names(mc38Treg) <- paste0('MC38_', names(mc38Treg))

allTreg <- cbind(b16Treg, mc38Treg) %>% 
  merge(annoDf, by.x = 0, by.y = 'peakID') %>% 
  rename('peakID' = 1) %>% 
  mutate(peakID = gsub('Interval', 'Peak', peakID))

write.csv(allTreg, '3.Output/5.Subtraction/subsectionRes.csv')

ggvenn::ggvenn(list(B16UP = b16Treg %>% dplyr::filter(B16_log2FoldChange > 0) %>% rownames(),
                    MC38UP = mc38Treg %>% filter(MC38_log2FoldChange > 0) %>% rownames(),
                    B16DN = b16Treg %>% filter(B16_log2FoldChange < 0) %>% rownames(),
                    MC38DN = mc38Treg %>% filter(MC38_log2FoldChange < 0) %>% rownames()), 
               set_name_size = 3.5, text_size = 2.5) + coord_fixed()
ggsave('3.Output/5.Subtraction/Venn_consistent_check.png', width = 6, height = 4)
