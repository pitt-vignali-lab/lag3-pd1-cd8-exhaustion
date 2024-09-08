# ----------------------------------------------------------------------------
# Fig 3 | Expression of the Klr gene family, including Klrc1 (NKG2A)
#
# Paper    : LAG-3 and PD-1 synergize on CD8+ T cells to drive T cell
#            exhaustion and hinder autocrine IFN-gamma-dependent anti-tumor
#            immunity
# Original : 2.Script/rev7_klrFamil.R
# Input    : 1.Data/datForRevision.RDS
# Output   : 3.Output/CellFigure/klrFamilyGeneExp*.png
# Upstream : 04_cluster_annotation.R
# ----------------------------------------------------------------------------

library(tidyverse)
library(Seurat)
library(presto)
library(export)

dat <- readRDS('1.Data/datForRevision.RDS')

datSub <- subset(dat, genotype %in% c("E8ICre.GFP (WT)", "Pdcd1L/L Lag3L/L-yfp (WT)", "Pdcd1L/L Lag3L/L-yfp E8ICre.GFP"))

datSub$genotype <- ifelse(datSub$genotype == "Pdcd1L/L Lag3L/L-yfp E8ICre.GFP", 
                          'DKO', 'WT')

Idents(datSub) <- 'genotype'

deG <- wilcoxauc(datSub) %>% 
  filter(group == 'DKO', 
         feature %in% c('Klrc1', 'Klrc2', 'Klrc3', 'Klrd1')) 


plots <- VlnPlot(datSub, features = c('Klrc1', 'Klrc2', 'Klrc3', 'Klrd1'), 
                  assay = 'RNA', raster = F)

plotDf <- data.frame()
for (i in 1:length(plots)) {
  plotDf <- rbind(plotDf,
                  plots[[i]]$data %>% 
                    rownames_to_column('cells') %>%
                    mutate(feature = names(.)[2]) %>% 
                    rename('read' = 2))
}

plotDf %>% ggplot(aes(x = ident, y = read)) +
  geom_violin(aes(fill = ident))+
  geom_boxplot(width = 0.1) +
  geom_jitter(shape = '`') +
  geom_text(data = deG, aes(x = 1.5, y = 4.5, 
                            label = paste0('padj =',  
                                           format(padj, 
                                           scientific = TRUE,
                                           digit = 2)))) +
  facet_grid(~feature) +
  theme_bw() +
  theme(legend.position = 'none')

graph2ppt(file = '3.Output/CellFigure/klrFamilyGeneExp.pptx', upscale = T)
ggsave('3.Output/CellFigure/klrFamilyGeneExp_Labled.png', width = 6, height = 3, bg = 'white')

plotDf %>% ggplot(aes(x = ident, y = read)) +
  geom_violin(aes(fill = ident))+
  geom_jitter(shape = '`') +
  geom_boxplot(width = 0.08) +
  facet_grid(~feature) +
  theme_bw() +
  coord_fixed() +
  theme(legend.position = 'none',
        text = element_blank())
ggsave('3.Output/CellFigure/klrFamilyGeneExp_NoLabled.png', width = 6, height = 3, bg = 'white')




