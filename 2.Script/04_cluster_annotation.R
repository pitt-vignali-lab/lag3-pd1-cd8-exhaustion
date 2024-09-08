# ----------------------------------------------------------------------------
# Setup | Annotate clusters into functional states
#
# Paper    : LAG-3 and PD-1 synergize on CD8+ T cells to drive T cell
#            exhaustion and hinder autocrine IFN-gamma-dependent anti-tumor
#            immunity
# Original : 2.Script/annotation.R
# Input    : 1.Data/datInt.RDS
# Output   : 1.Data/datForRevision.RDS
# Upstream : 02_cluster_labeling.R
# ----------------------------------------------------------------------------

library(Seurat)
library(tidyverse)
library(presto)

if(!dir.exists("3.Output/2.ClstLab")){
  dir.create("3.Output/2.ClstLab")
}

dat <- readRDS("1.Data/datInt.RDS")

mkGs <- c('Ptprc', 'Cd44', 'Ccr7', 'Sell', 'Il7r', 'Il2rb', 'Cd28',
          'Cd27', 'Spn', 'Fas', 'Klrg1', 'Prf1', 'Gzmb')

DotPlot(dat, features = mkGs, assay = 'RNA') +
  coord_flip(clip = 'off') +
  scale_x_discrete(limits=rev)
  

DotPlot(dat, features = c('Il7'), assay = 'RNA') +
  coord_flip(clip = 'off') +
  scale_x_discrete(limits=rev)

dat <- FindNeighbors(dat, dims = 1:15)

dat <- FindClusters(dat, resolution = 0.3,
                    method = 4)

DimPlot(dat, label = T) +
  geom_point(shape = '`')
  coord_fixed()

DimPlot(dat, split.by = 'seurat_clusters') +
  coord_fixed()

DotPlot(dat, 
            features = c('Cd69', 'Itgae', 'Itga1', 'Cxcr6', 'Sell'),
        cols = c('blue', 'red'),
        assay = "RNA")

DotPlot(dat, 
        features = c('Cxcr5'),
        cols = c('blue', 'red'),
        assay = "RNA")

DotPlot(dat, 
        features = c('Il4', 'Il5', 'Il13'))


DotPlot(dat, 
        features = c('Bst2', 'Irf1', 'Irf2', 'Irf7', 'Stat1', "Stat2", 'Mx1'))
