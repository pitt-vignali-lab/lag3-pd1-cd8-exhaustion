# ----------------------------------------------------------------------------
# Fig 3 | Differential expression within the exhausted cluster by genotype
#
# Paper    : LAG-3 and PD-1 synergize on CD8+ T cells to drive T cell
#            exhaustion and hinder autocrine IFN-gamma-dependent anti-tumor
#            immunity
# Original : 2.Script/rev3_clust4.R
# Input    : 1.Data/datForRevision.RDS
# Output   : 3.Output/Revision/clst4DeGsInEachGeno.csv
#            3.Output/Revision/clst4DeGsInEachGenoHtmp.png
# Upstream : 04_cluster_annotation.R
# ----------------------------------------------------------------------------

library(Seurat)
library(presto)
library(ComplexHeatmap)
library(miloR)

dat <- readRDS('1.Data/datForRevision.RDS')

datPml <- subset(dat, sample == 'PML CD8')

# cell number ration ----
t = cbind(data.frame(genotype = datPml$genotype,
           id = datPml$seurat_clusters) %>% 
  group_by(genotype) %>% 
  summarize(freq = n()) ,


data.frame(genotype = datPml$genotype,
           id = datPml$seurat_clusters) %>% 
  filter(id == 4) %>% 
  group_by(genotype) %>% 
  summarize(freq4 = n()))

 
t$ratio = t$freq4/t$freq

ggplot(t[,-3], aes(x = genotype,  y = ratio * 100)) +
  geom_bar(stat = "identity")

# de genes in each genotype ----
datPml4 <- subset(dat, idents = 4, sample == 'PML CD8')

Idents(datPml4) <- 'genotype'

deGSub <- wilcoxauc(datPml4) %>% 
  filter(padj < 0.05)

feats <- unique(deGSub$feature)[-grep('^mt-', unique(deGSub$feature))]

deGSub %>% filter(feature %in% feats) %>% 
  write.csv('3.Output/Revision/clst4DeGsInEachGeno.csv')

length(unique(deGSub$feature))



deMat <- AverageExpression(datPml4, features = feats, 
                           assays = 'RNA')[[1]]

htmp <- BBmisc::normalize(deMat, range = c(-2, 2))

colnames(htmp) <- colnames(deMat)

png('3.Output/Revision/clst4DeGsInEachGenoHtmp.png', width = 5, height = 12,
    units = 'in', res = 300)
Heatmap(htmp,
        column_names_rot = 45,
        column_names_gp = gpar(fontsize = 8),
        row_names_gp = gpar(fontsize = 5))

dev.off()


dotPDf <- DotPlot(datPml4, features = feats, assay = 'RNA')$data

ggplot(dotPDf, aes(x = id, y = features.plot)) +
  geom_point(aes(fill = avg.exp.scaled, size = pct.exp),
             shape = 21) +
  scale_fill_gradient2(high = 'red', mid = 'white', low = 'blue', midpoint = 0)

# clust4 vs clust6 ----
datSub64 <- subset(dat, idents = c(4, 6))

deG64 <- wilcoxauc(datSub64) %>% 
  filter(padj < 0.05,
         group == 4)

datSub64$sampGeno <- paste(datSub64$sample, datSub64$genotype)
datSub64$sampGenoAnno <- paste(datSub64$sample, datSub64$genotype, datSub64$anno)

datSce <- as.SingleCellExperiment(datSub64)

milo <- Milo(datSce)

# reducedDim(milo, "UMAP") <- reducedDim(datSce, "UMAP")

milo

milo <- buildGraph(milo, k = 50, d = 30)

milo <- makeNhoods(milo, prop = 0.2, 
                   k = 50, 
                   d = 30, 
                   refined = TRUE)

plotNhoodSizeHist(milo)

milo <- countCells(milo, 
                   meta.data = data.frame(colData(milo)), 
                   sample="sampGenoAnno")

head(nhoodCounts(milo))

designDf <- data.frame(colData(milo))[,c('sample', 'genotype', 
                                         'anno', 'sampGeno', 'sampGenoAnno')] %>% 
  distinct() %>% 
  remove_rownames() %>% 
  # mutate(tissueType = paste(tissue, anno, sep = '_')) %>% 
  column_to_rownames(var = "sampGenoAnno") 

# designDf <- designDf[colnames(nhoodCounts(milo)), , drop=FALSE]

# designDf$sampGenoAnno <- rownames(designDf)

designDf

milo <- calcNhoodDistance(milo, d = 30)

# DA test 
da_results <- testNhoods(milo, 
                         design = ~ anno, 
                         design.df = designDf)
head(da_results)

da_results %>%
  arrange(SpatialFDR) %>%
  head() 

ggplot(da_results, aes(PValue)) + geom_histogram(bins=50)

ggplot(da_results, aes(logFC, -log10(SpatialFDR))) + 
  geom_point() +
  geom_hline(yintercept = 1) ## Mark significance threshold (10% FDR)

milo <- buildNhoodGraph(milo)

plotReducedDim(milo, dimred = "UMAP", colour_by="anno", 
               text_size = 3, point_size=0.5) +
  guides(fill="none") +
  coord_fixed() +

plotNhoodGraphDA(milo, 
                 da_results, 
                 layout="UMAP",
                 alpha=0.1) +
  scale_fill_gradient2(high = 'red', mid = 'white', low = 'blue', midpoint = 0) +
  labs(fill = 'logFC') +
  coord_fixed() 

da_results <- annotateNhoods(milo, 
                             da_results, 
                             coldata_col = "anno")


da_results$NhoodGroup <- da_results$anno

da_nhood_markers <- findNhoodGroupMarkers(milo, 
                                          da_results, 
                                          subset.row = rownames(milo))

mk <- da_nhood_markers %>% 
  filter(`adj.P.Val_Tumor primed-TCM` < 0.05)

inteG <- intersect(mk$GeneID, deG64$feature)

inteG[-grep('^Hist', inteG)]

write.csv(deG64[deG64$feature %in% inteG[-grep('^Hist', inteG)],],
          '3.Output/Revision/milo&DeOfClst64.csv')
