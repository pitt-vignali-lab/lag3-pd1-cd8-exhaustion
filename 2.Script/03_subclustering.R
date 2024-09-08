# ----------------------------------------------------------------------------
# Fig 2 | Re-cluster the clusters of interest for finer structure
#
# Paper    : LAG-3 and PD-1 synergize on CD8+ T cells to drive T cell
#            exhaustion and hinder autocrine IFN-gamma-dependent anti-tumor
#            immunity
# Original : 2.Script/3.sbuClst.R
# Input    : 1.Data/datInt.RDS
# Output   : 1.Data/datIntNoSmallClst.RDS
#            3.Output/3.SubClst/
# Upstream : 02_cluster_labeling.R
# ----------------------------------------------------------------------------

library(Seurat)
library(tidyverse)
library(clustree)
library(presto)

if(!dir.exists("3.Output/3.SubClst")){
  dir.create("3.Output/3.SubClst")
}

dat <- readRDS("1.Data/datInt.RDS")

# sub clustering multiple clusters ----
datSub <- subset(dat, idents = c(0, 7))

datSub=CreateSeuratObject(
  counts = datSub@assays$RNA@counts,
  meta.data = datSub@meta.data)

datList <- SplitObject(datSub, split.by = "sample")
datList <- lapply(datList, function(x){
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, verbose = FALSE)
})
features <- SelectIntegrationFeatures(object.list = datList)
datList <- lapply(X = datList, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})
anchors <- FindIntegrationAnchors(object.list = datList, 
                                  reference = c(1, 2), 
                                  reduction = "rpca",
                                  dims = 1:50)
datInt <- IntegrateData(anchorset = anchors, dims = 1:50)
datInt <- ScaleData(datInt, verbose = FALSE)
datInt <- RunPCA(datInt, verbose = FALSE)


pct <- datInt[["pca"]]@stdev/sum(datInt[["pca"]]@stdev) * 100
cumu <- cumsum(pct)
co1 <- which(cumu > 90 & pct < 5)[1]
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), 
            decreasing = T)[1] + 1
pcs <- min(co1, co2)
plot_df <- data.frame(pct = pct, 
                      cumu = cumu, 
                      rank = 1:length(pct))
ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + 
  geom_text() + 
  geom_vline(xintercept = 90, color = "grey") + 
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()


datInt <- RunUMAP(datInt, dims = 1:pcs)
datInt <- FindNeighbors(datInt, dims = 1:pcs)
datInt <- FindClusters(datInt, resolution = seq(0, 0.5, 0.05))

clustree(datInt, prefix = "integrated_snn_res.")

Idents(datInt) <- 'integrated_snn_res.0.3'

table(Idents(datInt))

# figure 
datSub <- subset(dat, idents = c(0, 7))
datSub <- AddMetaData(datSub, data.frame(Idents(datInt)))
# DimPlot(datSub, group.by = "Idents.datInt.")

df <- data.frame(x = datSub[["umap"]]@cell.embeddings[,1],
                 y = datSub[["umap"]]@cell.embeddings[,2],
                 sample = datSub$sample,
                 tissue = datSub$tissue,
                 genotype = datSub$genotype,
                 oldLab = Idents(datSub),
                 newLab = datSub$Idents.datInt.)

df$genotype <- factor(df$genotype, levels = c("E8ICre.GFP (WT)",
                                              "Pdcd1L/L Lag3L/L-yfp (WT)",
                                              "Pdcd1L/L E8ICre.GFP",
                                              "Lag3L/L-yfp E8ICre.GFP",
                                              "Pdcd1L/L Lag3L/L-yfp E8ICre.GFP"))


dfLab <- df %>% group_by(newLab) %>% 
  summarise(x = mean(x), y = mean(y))

ggplot() +
  geom_point(data = df, aes(x , y, col = newLab), 
             shape=16, alpha = 1, size = 2, stroke = 0) +
  geom_text(data = dfLab, aes(x , y, label = newLab), 
            size=8, col = "black") +
  xlab("UMAP_1") +
  ylab("UMAP_2") +
  theme_bw() +
  theme(text = element_text(face = 2),    
        legend.position="none")

ggplot() +
  geom_point(data = df, aes(x , y, col = newLab), 
             shape=16, alpha = 1, size = 2, stroke = 0) +
  facet_grid(sample + tissue ~ genotype, switch = "y") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        # panel.grid.minor = element_blank(),
        # panel.grid.major = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 15, 
                                  face = 2))


write.csv(wilcoxauc(datInt), "3.Output/3.SubClst/markerGsforsubClst_0_7.csv")
write.csv(data.frame(Idents(datInt)),
          "3.Output/3.SubClst/clstrsforsubClst_0_7.csv")


clst <- c(1)

datSub <- subset(dat, idents = clst)

datSub=CreateSeuratObject(
  counts = datSub@assays$RNA@counts,
  meta.data = datSub@meta.data)

datList <- SplitObject(datSub, split.by = "sample")
datList <- lapply(datList, function(x){
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, verbose = FALSE)
})
features <- SelectIntegrationFeatures(object.list = datList)
datList <- lapply(X = datList, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})
anchors <- FindIntegrationAnchors(object.list = datList, 
                                  reference = c(1, 2), 
                                  reduction = "rpca",
                                  dims = 1:50)
datInt <- IntegrateData(anchorset = anchors, dims = 1:50)
datInt <- ScaleData(datInt, verbose = FALSE)
datInt <- RunPCA(datInt, verbose = FALSE)


pct <- datInt[["pca"]]@stdev/sum(datInt[["pca"]]@stdev) * 100
cumu <- cumsum(pct)
co1 <- which(cumu > 90 & pct < 5)[1]
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), 
            decreasing = T)[1] + 1
pcs <- min(co1, co2)
plot_df <- data.frame(pct = pct, 
                      cumu = cumu, 
                      rank = 1:length(pct))
ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + 
  geom_text() + 
  geom_vline(xintercept = 90, color = "grey") + 
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()


datInt <- RunUMAP(datInt, dims = 1:pcs)
datInt <- FindNeighbors(datInt, dims = 1:pcs)
datInt <- FindClusters(datInt, resolution = seq(0, 0.5, 0.05))

clustree(datInt, prefix = "integrated_snn_res.")
Idents(datInt) <- 'integrated_snn_res.0.15'
table(Idents(datInt))
DimPlot(datInt)

# figure 
datSub <- subset(dat, idents = clst)
datSub <- AddMetaData(datSub, data.frame(Idents(datInt)))
DimPlot(datSub, group.by = "Idents.datInt.")

df <- data.frame(x = datSub[["umap"]]@cell.embeddings[,1],
                 y = datSub[["umap"]]@cell.embeddings[,2],
                 sample = datSub$sample,
                 tissue = datSub$tissue,
                 genotype = datSub$genotype,
                 oldLab = Idents(datSub),
                 newLab = datSub$Idents.datInt.)

df$genotype <- factor(df$genotype, levels = c("E8ICre.GFP (WT)",
                                              "Pdcd1L/L Lag3L/L-yfp (WT)",
                                              "Pdcd1L/L E8ICre.GFP",
                                              "Lag3L/L-yfp E8ICre.GFP",
                                              "Pdcd1L/L Lag3L/L-yfp E8ICre.GFP"))


dfLab <- df %>% group_by(newLab) %>% 
  summarise(x = mean(x), y = mean(y))

ggplot() +
  geom_point(data = df, aes(x , y, col = newLab), 
             shape=16, alpha = 1, size = 2, stroke = 0) +
  geom_text(data = dfLab, aes(x , y, label = newLab), 
            size=8, col = "black") +
  xlab("UMAP_1") +
  ylab("UMAP_2") +
  theme_bw() +
  theme(text = element_text(face = 2),    
        legend.position="none")

ggplot() +
  geom_point(data = df, aes(x , y, col = newLab), 
             shape=16, alpha = 1, size = 2, stroke = 0) +
  facet_grid(sample + tissue ~ genotype, switch = "y") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 15, 
                                  face = 2))

# write file
write.csv(wilcoxauc(datInt), 
          paste0("3.Output/3.SubClst/markerGsforsubClst_", 
                 paste0(clst, collapse = "_"), ".csv"))
write.csv(data.frame(Idents(datInt)),
          paste0("3.Output/3.SubClst/clstrsforsubClst_", 
                 paste0(clst, collapse = "_"), ".csv"))

# single sample ----
datSub <- subset(dat, idents = c(7), sample == "Normal CD8")

datSub=CreateSeuratObject(
  counts = datSub@assays$RNA@counts,
  meta.data = datSub@meta.data)

datSub <- NormalizeData(datSub) %>% 
  FindVariableFeatures() %>%
  ScaleData() %>% 
  RunPCA()

datSub <- RunUMAP(datSub)
datSub <- FindNeighbors(datSub)
datSub <- FindClusters(datSub, resolution = seq(0, 0.5, 0.05))

clustree(datSub, prefix = "RNA_snn_res.")

Idents(datSub) <- "RNA_snn_res.0.15"

DotPlot(datSub, 
        features = c(strsplit("Bst2, Irf1, Irf2, Irf7, Stat1, Stat2, Mx1", ", ")[[1]]#, # high
                     # strsplit("Tcf7, Havcr2", ", ")[[1]] # low
        ), 
        col.min = 0, cols = c("white", "blue"), assay = "RNA") 

DimPlot(datSub, split.by = "genotype")