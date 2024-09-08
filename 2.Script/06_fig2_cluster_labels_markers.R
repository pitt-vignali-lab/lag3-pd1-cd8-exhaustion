# ----------------------------------------------------------------------------
# Fig 2 | Final cluster labels and marker gene panels
#
# Paper    : LAG-3 and PD-1 synergize on CD8+ T cells to drive T cell
#            exhaustion and hinder autocrine IFN-gamma-dependent anti-tumor
#            immunity
# Original : 2.Script/rev1_clusterLab&mkG.R
# Input    : 1.Data/datForRevision.RDS
# Output   : 3.Output/CellFigure/mkGexpression.png
#            3.Output/CellFigure/markerGenes.pptx
# Upstream : 04_cluster_annotation.R
# ----------------------------------------------------------------------------

if(!dir.exists('3.Output/Revision')){
  dir.create('3.Output/Revision')
}

library(Seurat)
library(tidyverse)
library(presto)
library(COSG)
library(export)

dat <- readRDS("1.Data/datInt.RDS")

dat <- subset(dat, idents = 9, invert = T)

# dat <- RunUMAP(dat, dims = 1:10)
# saveRDS(dat, '1.Data/datIntRv.RDS')
# 
# DimPlot(dat) +
#   geom_hline(yintercept = -8) +
#   coord_fixed()

umapAx <- Embeddings(dat, 'umap') %>% 
  as.data.frame()

dat <- subset(dat, cells = rownames(umapAx[umapAx$UMAP_2 > -8,]))

DimPlot(dat) +
  coord_fixed()

# marker genes  ----

mkG <- cosg(dat)

hmks <- rbind(data.frame(type = 'C3: Proliferation', # 3
                         mk = c('Top2a', 'Mki67')),
              data.frame(type = 'C2: Memroy \nprecursors',
                         mk = c('Tcf7', 'Sell', 'Lef1', 'Il7r')),
              data.frame(type = 'C4: Exhuasted/ \n C6: Pre-exhuasted',
                         mk = c('Ctla4', 'Layn', 'Entpd1', 'Havcr2', 'Tox', 'Mt1')),
              data.frame(type = 'C0: Effector/effector memory-like',
                         mk = c('Ccl5', 'Cxcr6', 'Gzma', 'Gzmb', 'Prf1', 'Nkg7')),
              data.frame(type = 'C7: IFN-induced',
                         mk = c('Bst2', 'Isg15', 'Stat1', 'Irf7')),
              data.frame(type = 'C5: Co-stumulatory',
                         mk = c('Tnfsf14', 'Tnfrsf9', 'Tnfrsf4')),
              data.frame(type = 'C1: Activated',
                         mk = c('Tfdp1', 'E2f1', 'Mcm7', 'Mcm2')),
              data.frame(type = 'C8: Activated (low function)',
                         mk = c('Tfdp1', 'E2f1', 'Mcm7', 'Mcm2')))


DotPlot(dat, features = unique(hmks$mk), assay = 'RNA') +
  scale_color_gradient2(high = 'red', mid = 'white', low = 'blue', midpoint = 0)

# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5297183/

# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8220296/

# https://www.biocompare.com/Editorial-Articles/569888-A-Guide-to-T-Cell-Markers/


annoDf <- data.frame(id = 0:8,
                     anno = c('C0: Effector/effector memory-like', 
                              'C1: Activated',
                              'C2: Memory precursors', 
                              'C3: Proliferation', # effector intermediate
                              'C4: Exhuasted',
                              'C5: Co-stumulatory',
                              'C6: Pre-exhuasted',
                              'C7: IFN-induced',
                              'C8: Activated(low function)')) %>%
  mutate(anno = factor(anno, levels = c('C0: Effector/effector memory-like', 
                                        'C1: Activated',
                                        'C2: Memory precursors', 
                                        'C3: Proliferation', # effector intermediate
                                        'C4: Exhuasted',
                                        'C5: Co-stumulatory',
                                        'C6: Pre-exhuasted',
                                        'C7: IFN-induced',
                                        'C8: Activated(low function)'))) %>% 
  merge(data.frame(id = dat$seurat_clusters) %>% 
          rownames_to_column('cells'),
        by = 'id') %>% 
  dplyr::select(c(anno, cells)) %>% 
  column_to_rownames('cells')

dat <- AddMetaData(dat, annoDf)

saveRDS(dat, '1.Data/datForRevision.RDS')

DimPlot(dat, label.size = 6, pt.size = 1.5, 
        group.by = 'anno', repel = T) + 
  coord_fixed() +
  labs(title = '') +
  theme_void()

graph2ppt(file = '3.Output/CellFigure/OverallUMAP.ppt', 
          orient = ifelse(type[1] == "PPT", "landscape", "auto"),
          margins = c(top = 0, right = 0, bottom = 0, left = 0))

ggsave('3.Output/Revision/allUMAP.png', width = 10, height = 8, bg = 'white', units = 'in')


dotPDf <- DotPlot(dat, features = unique(hmks$mk), assay = 'RNA', group.by = 'anno')$data 

ggplot(dotPDf, aes(x = id, y = features.plot)) +
  geom_point(aes(size = pct.exp, fill = avg.exp.scaled), shape = 21) +
  scale_fill_gradient2(high = 'red', mid = 'white', low = 'blue', midpoint = 0) +
  coord_fixed(ratio = 0.3) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1),
        text = element_text(size = 15, face = 2, colour = 'black'))

graph2ppt(file = '3.Output/CellFigure/markerGenes.pptx',
          margins = c(top = 0, right = 0, bottom = 0, left = 0))

ggsave('3.Output/Revision/mkGExpression.png', width = 10, height = 8, bg = 'white', units = 'in')






gList <- unique(unlist(as.vector(mkG$names)))

dotPDf <- DotPlot(dat, features = gList, assay = 'RNA', group.by = 'anno')$data

ggplot(dotPDf, aes(x = id, y = features.plot)) +
  geom_point(aes(size = pct.exp, fill = avg.exp.scaled), shape = 21) +
  scale_fill_gradient2(high = 'red', mid = 'white', low = 'blue', midpoint = 0) +
  coord_fixed(ratio = 0.2) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1),
        text = element_text(face = 2, colour = 'black'))

gList <- strsplit('Cxcr6, Ccl5, S100a6, S100a4, Apex1, E2f1, Sell, Tcf7, Il7r, Mki67, Gzmb, Tpx2, Tox, Tox2, Ccr7, Cd9, Id3, Tnfsf8, Macf1, Lars2, Ifit1, Irf1, Stat1, Gzmc, Ccne1, Ccl1', ', ')[[1]]

dotPDf <- DotPlot(dat, features = gList, assay = 'RNA', group.by = 'anno')$data

ggplot(dotPDf, aes(y = id, x = features.plot)) +
  geom_point(aes(size = pct.exp, fill = avg.exp.scaled), shape = 21) +
  scale_fill_gradient2(high = 'red', mid = 'white', low = 'blue', midpoint = 0) +
  scale_y_discrete(limits=rev) +
  coord_fixed(ratio = 1) +
  labs(y = '', x = '') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1),
        text = element_text(size = 20, face = 2, colour = 'black'),
        legend.position = 'bottom')
ggsave('3.Output/CellFigure/mkGexpression.png', width = 12, height = 7, bg = 'white', units = 'in')
ggsave('3.Output/Revision/mkGExpression.png', width = 10, height = 8, bg = 'white', units = 'in')
