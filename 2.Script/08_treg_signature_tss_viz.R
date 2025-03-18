# ----------------------------------------------------------------------------
# Fig 2 | Accessibility of Treg transcriptional signature genes around the TSS
#
# Paper    : Regulatory T cells in the tumor microenvironment display a
#            unique chromatin accessibility profile
# Original : 2.Script/6.tregTransSigViz.R
# Input    : 1.Data/Treg_effector_sig.xlsx
#            3.Output/3.PeakAnno/
# Output   : 3.Output/6.tregSignViz/*geneTSS.pdf
# Upstream : 04_atac_peak_annotation.R
# ----------------------------------------------------------------------------

if(!dir.exists('3.Output/6.tregSignViz')){
  dir.create('3.Output/6.tregSignViz')
}

library(readxl)
library(tidyverse)


peakAnno <- readRDS('3.Output/3.PeakAnno/peakAnnotationChipseeker.RDS')
annoDf <- data.frame(peakAnno@anno)

tsig <- read_xlsx('1.Data/Treg_effector_sig.xlsx')[-c(1:6), 1:7] %>% 
  janitor::row_to_names(1) %>% 
  as.data.frame() %>% 
  mutate(adj.P.Val = as.numeric(adj.P.Val),
         logFC = as.numeric(logFC)) %>% 
  filter(adj.P.Val < 0.05)

for (i in c('B16', 'MC38')) {
  dfDe <- readRDS(paste0('3.Output/2.DEseqRes/1.B16MC38TregCD4/de_', i, '_Treg.RDS')) %>% 
    as.data.frame() %>% 
    merge(annoDf, by.x = 0, by.y = 'peakID') %>% 
    group_by(SYMBOL) %>% 
    filter(distanceToTSS == min(abs(distanceToTSS)),
           padj < 0.05) %>% 
    group_by(SYMBOL) %>%
    filter(padj == min(padj)) %>% 
    select("SYMBOL", "log2FoldChange", 'padj')
  
  plotDf  <- merge(dfDe, tsig, by.x = 'SYMBOL', by.y = 'Symbol') %>% 
    mutate(cat = case_when(
      (log2FoldChange > 0 & logFC > 0) ~ 'up',
      (log2FoldChange < 0 & logFC < 0) ~ 'dn',
      T ~ 'diff'
      
    ))
  
  labDf <- rbind(plotDf %>% 
    mutate(comboValue = log2FoldChange * logFC) %>% 
    filter(cat != 'diff') %>% 
    group_by(cat) %>% 
    top_n(10, comboValue),
    plotDf %>% filter(SYMBOL %in% c('Satb1', 'Tgfbr3','Gramd3', 'Klrd1')))
  
  prop.table(table(plotDf$cat))
  
  ggplot(plotDf, aes(x = log2FoldChange, y = logFC))  +
    geom_point(aes(col = cat)) +
    ggrepel::geom_text_repel(data = labDf, aes(label = SYMBOL), fontface = 2) +
    # ggrepel::geom_text_repel(aes(label = SYMBOL), fontface = 2) +
    scale_color_manual(values = alpha(c('gray', 'blue', 'red'), 0.5)) +
    coord_fixed() +
    guides(color = 'none') +
    labs(x = paste0(i, ' gene TSS'),
         y = 'Treg signature gene') +
    theme_bw() +
    theme(# aspect.ratio=1,
          axis.text = element_text(face = 2, size = 12),
          text = element_text(face = 2, size = 12))
  
  ggsave(paste0('3.Output/6.tregSignViz/', i, 'geneTSS.png'), width = 6, height = 8, units = 'in')
  ggsave(paste0('3.Output/6.tregSignViz/', i, 'geneTSS.pdf'), width = 6, height = 8)
  
}
