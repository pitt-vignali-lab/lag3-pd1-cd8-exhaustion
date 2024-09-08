# ----------------------------------------------------------------------------
# Fig 2e-g | Prepare loom input and export embeddings for RNA velocity
#
# Paper    : LAG-3 and PD-1 synergize on CD8+ T cells to drive T cell
#            exhaustion and hinder autocrine IFN-gamma-dependent anti-tumor
#            immunity
# Original : 2.Script/new_scVelo/1.init.R
# Input    : 1.Data/datForRevision.RDS
#            1.Data/Velo/*.loom
# Output   : coord.csv and col.csv exports read by the velocity notebooks
# Upstream : 04_cluster_annotation.R
# Note     : Spliced / unspliced counts produced with Velocyto 0.17
# ----------------------------------------------------------------------------

library(Seurat)
library(scales)
library(tidyverse)

dat <- readRDS("1.Data/datForRevision.RDS")


# dat <- readRDS('1.Data/datIntRv.RDS')
dat$anno <- Idents(dat)

table(dat$sample)

DimPlot(dat) + coord_fixed()

for (sub in c('overall', 'Normal CD8', 'PML CD8')) {
  
  if (sub == 'overall') {
    datSub <- dat} else {
      datSub <- subset(dat, sample == sub)}
  
  # DimPlot(dat) +
  #   coord_fixed()
  
  # cell barcodes 
  write.csv(Cells(datSub), 
            file = paste(sub, "cellID.csv", sep = '_'), 
            row.names = FALSE)
  
  # umap 
  Embeddings(datSub, reduction = "umap") %>% 
    write.csv(file = paste(sub, "umap.csv", sep = '_'))
  
  # cluster
  datSub$anno %>% 
    write.csv(file = paste(sub, "clusters.csv", sep = '_'))
  
  # color collection 
  data.frame(id = levels(datSub$anno),
             col = hue_pal()(length(table(Idents(datSub)))))%>% 
    arrange(id) %>% 
    write.csv(file = paste(sub, "col.csv", sep = '_'))
  
  # figure size
  plt <- DimPlot(datSub) +
    coord_fixed()
  
  # x-range
  x = round(diff(layer_scales(plt)$x$range$range), 1)
  
  # y-range
  y = round(diff(layer_scales(plt)$y$range$range), 1)
  
  write.csv(data.frame(x, y), 
            file =  paste(sub, "coord.csv", sep = '_'))
}


# old ----
# newClst <- read.csv("1.Data/TCR_defined_clusters.csv", row.names = 1)
# names(newClst) <- "newClst"

# dat <- AddMetaData(dat, newClst)

# dat <- subset(dat, sample == "PML CD8")
# Idents(dat) <- "newClst"

DimPlot(dat)

write.csv(Cells(dat), file = "cellID.csv", row.names = FALSE)

embd <- Embeddings(dat, reduction = "umap")
write.csv(embd, file = "umap.csv")

clst <- Idents(dat) # dat@meta.data["newClst"]
write.csv(clst, file = "clusters.csv")

# Idents(dat) <- factor(Idents(dat), levels = 0:10)

hue_pal()(length(table(Idents(dat))))

plt <- DimPlot(dat) +
  coord_fixed()

plt <- DimPlot(subset(dat, idents = c(0, 7))) +
  coord_fixed()

# x-range
round(diff(layer_scales(plt)$x$range$range), 1)

# y-range
round(diff(layer_scales(plt)$y$range$range), 1)


head(Cells(dat))

plt <- DimPlot(subset(dat, idents = c(0, 7))) +
  coord_fixed()
