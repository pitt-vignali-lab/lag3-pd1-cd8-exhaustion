# ----------------------------------------------------------------------------
# Setup | Integrate the polyclonal and pMEL scRNA-seq datasets
#
# Paper    : LAG-3 and PD-1 synergize on CD8+ T cells to drive T cell
#            exhaustion and hinder autocrine IFN-gamma-dependent anti-tumor
#            immunity
# Original : 2.Script/1. datIntegration.R
# Input    : 1.Data/202104.RDS
#            1.Data/202201.RDS
# Output   : 1.Data/datInt.RDS
#            3.Output/1.DatInit/
# Note     : Reciprocal PCA integration
#            PC10 cutoff, resolution 0.35, 10 clusters. Samples demultiplexed from TotalSeq-C hashes with HTODemux (positive.quantile 0.999)
# ----------------------------------------------------------------------------

library(Seurat)
library(tidyverse)

if(!dir.exists("3.Output/1.DatInit")){
  dir.create("3.Output/1.DatInit")
}

datNor <- readRDS("1.Data/202104.RDS")

#TCR 
tcr <- read.csv("1.Data/filtered_contig_annotations.csv")
tcr$barcode <- gsub("-1", "", tcr$barcode)
tcr <- tcr[!duplicated(tcr$barcode), ]
rownames(tcr) <- tcr[,1]
datNor <- AddMetaData(object = datNor, metadata =tcr)


datNor <- subset(datNor, tissues == "Tumor")

datNor$genotype <- datNor$types
datNor$tissue <- datNor$tissues
datNor$sample <- "Normal CD8"

datPml <- readRDS("1.Data/202201.RDS")

# add genotype column
datPml <- AddMetaData(datPml, merge(datPml@meta.data["sample"] %>% rownames_to_column(var = "cell"),
                            data.frame(sample = unique(datPml@meta.data$sample),
                                       genotype = c("Pdcd1L/L Lag3L/L-yfp E8ICre.GFP",
                                                    "Pdcd1L/L E8ICre.GFP",
                                                    "Pdcd1L/L Lag3L/L-yfp (WT)",
                                                    "Lag3L/L-yfp E8ICre.GFP")),
                            by = "sample") %>% select(-sample) %>% column_to_rownames(var = "cell"))

datPml$tissue <- "Tumor"
datPml$sample <- "PML CD8"

dat <- merge(datNor, datPml, add.cell.ids = c("nor", "pml"))

# QC ----
dat$log10GenesPerUMI <- log10(dat$nFeature_RNA) / log10(dat$nCount_RNA)
dat$mitoRatio <- PercentageFeatureSet(object = dat, pattern = "^mt-")/100

dat@meta.data %>% 
  ggplot(aes(x=sample, fill=sample)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")


dat@meta.data %>% 
  ggplot(aes(color=sample, x=nCount_RNA, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  expand_limits(x = 0) +
  geom_vline(xintercept = 300) + 
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("nCount_RNA")


dat@meta.data %>% 
  ggplot(aes(color=sample, x=nFeature_RNA, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 200) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("nFeature_RNA")


dat@meta.data %>% 
  ggplot(aes(y=log10(nFeature_RNA), x= sample, fill= sample)) + 
  geom_boxplot(alpha = 0.2) + 
  # scale_x_log10() + 
  theme_classic() +
  ylab("lgo10(nGene)") +
  geom_vline(xintercept = 300) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("nFeature_RNA")

dat@meta.data %>% 
  ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250) +
  facet_wrap(~sample) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("nCount_RNA vs. nFeature_RNA")

dat@meta.data %>% 
  ggplot(aes(color=sample, x=mitoRatio, fill=sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 0.5) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Mitochondrial counts ratio")

dat@meta.data %>%
  ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.80) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Complexity")

# batch effects ----
datList <- SplitObject(dat, split.by = "sample")

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


pct <- datInt[["pca"]]@stdev / sum(datInt[["pca"]]@stdev) * 100

cumu <- cumsum(pct)

co1 <- which(cumu > 90 & pct < 5)[1]

co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1

pcs <- min(co1, co2) # 10

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

datInt <- FindClusters(datInt, resolution = 0.35)

# saveRDS(datInt, "1.Data/datInt.RDS")

# clustree(datInt, prefix = "integrated_snn_res.")
# 
# Idents(datInt) <- 'integrated_snn_res.0.35'

table(Idents(datInt))

DimPlot(datInt, label = T)

DimPlot(datInt, group.by = "sample")

# figure generation ----

df <- data.frame(x = datInt[["umap"]]@cell.embeddings[,1],
                 y = datInt[["umap"]]@cell.embeddings[,2],
                 sample = datInt$sample,
                 tissue = datInt$tissue,
                 genotype = datInt$genotype,
                 lab = Idents(datInt))

df$genotype <- factor(df$genotype, levels = c("E8ICre.GFP (WT)",
                                              "Pdcd1L/L Lag3L/L-yfp (WT)",
                                              "Pdcd1L/L E8ICre.GFP",
                                              "Lag3L/L-yfp E8ICre.GFP",
                                              "Pdcd1L/L Lag3L/L-yfp E8ICre.GFP"))

dfLab <- df %>% group_by(lab) %>% 
  summarise(x = mean(x), y = mean(y))

ggplot() +
  geom_point(data = df, aes(x , y, col = lab, fill = lab), 
             shape=16, alpha = 0.5, size = 3, stroke = 0) +
  geom_text(data = dfLab, aes(x , y, label = lab), 
            size=8, col = "black") +
  xlab("UMAP_1") +
  ylab("UMAP_2") +
  theme_void() +
  theme(text = element_text(face = 2),    
        legend.position="none")

ggsave("3.Output/1.DatInit/indexUmap.png", width = 8, height = 8, 
       bg = "white", dpi = 300)


ggplot() +
  geom_point(data = df, aes(x , y, col = lab, fill = lab), 
             shape=16, alpha = 0.5, size = 5, stroke = 0) +
  facet_grid(sample + tissue ~ genotype, switch = "y") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 15, 
                                  face = 2))

ggsave("3.Output/1.DatInit/splitUmap.png", width = 20, height = 8, 
       bg = "white", dpi = 300)

# remove mini tiny cluster ----

datInt <- subset(datInt, idents = 9, invert = T)

datInt <- RunUMAP(datInt, dims = 1:20,
                  umap.method = 'uwot-learn')

DimPlot(datInt, label = T) +
  coord_fixed()

# saveRDS(datInt, "1.Data/datIntNoSmallClst.RDS")

# bar ----

datSub <- subset(datInt, sample == "Normal CD8" & tissue == "Tumor")

t1 <- data.frame(prop.table(table(Idents(datSub), 
                                 datSub$genotype)))

datSub <- subset(datInt, sample == "PML CD8" & tissue == "Tumor")

t2 <- data.frame(prop.table(table(Idents(datSub), 
                                  datSub$genotype)))


t1$sample <- "Normal CD8"
t2$sample <- "PML CD8"

t <- rbind(t1, t2)

names(t) <- c("Cluster", "Genotype", "Freq", "Sample")

t$Genotype <- factor(t$Genotype, levels = rev(c("E8ICre.GFP (WT)",
                                              "Pdcd1L/L Lag3L/L-yfp (WT)",
                                              "Pdcd1L/L E8ICre.GFP",
                                              "Lag3L/L-yfp E8ICre.GFP",
                                              "Pdcd1L/L Lag3L/L-yfp E8ICre.GFP")))

ggplot(t, aes(x = Freq, y = Genotype, fill = Cluster)) +
  geom_bar(col = "black", 
           stat = "identity",
           position=position_fill(reverse=T)) +
  scale_x_continuous(expand = c(0.015, 0)) +
  facet_wrap(~Sample) +
  theme(text = element_text(size = 18, face = 2),
        legend.position = "bottom", 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.ticks.y = element_blank(),
        strip.background = element_blank())+
  guides(fill = guide_legend(nrow = 2))

ggsave("3.Output/1.DatInit/celldistBarbyGeno.png", width = 30, height = 5, 
       bg = "white", dpi = 300)

t$Genotype <- factor(t$Genotype, levels = c("E8ICre.GFP (WT)",
                                                "Pdcd1L/L Lag3L/L-yfp (WT)",
                                                "Pdcd1L/L E8ICre.GFP",
                                                "Lag3L/L-yfp E8ICre.GFP",
                                                "Pdcd1L/L Lag3L/L-yfp E8ICre.GFP"))
t$Cluster <- factor(t$Cluster, levels = 9:0)

ggplot(t, aes(x = Freq, y = Cluster, fill = Genotype)) +
  geom_bar(col = "black", 
           stat = "identity",
           position=position_fill(reverse=T)) +
  scale_x_continuous(expand = c(0.018, 0)) +
  facet_wrap(~Sample) +
  theme(text = element_text(size = 18, face = 2),
        legend.position = "bottom", 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.ticks.y = element_blank(),
        strip.background = element_blank())+
  guides(fill = guide_legend(nrow = 2))

ggsave("3.Output/1.DatInit/celldistBarbyClst.png", width = 20, height = 10, 
       bg = "white", dpi = 300)

DimPlot(datInt)

DimPlot(datInt, split.by = "genotype")
