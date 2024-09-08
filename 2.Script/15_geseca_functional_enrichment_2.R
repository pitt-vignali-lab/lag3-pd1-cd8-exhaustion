# ----------------------------------------------------------------------------
# Fig 2 | GESECA follow-up and simplified GO view
#
# Paper    : LAG-3 and PD-1 synergize on CD8+ T cells to drive T cell
#            exhaustion and hinder autocrine IFN-gamma-dependent anti-tumor
#            immunity
# Original : 2.Script/rev4_funcEnrich1.R
# Input    : 3.Output/Revision/funcEnrich.RDS
# Output   : 3.Output/Revision/goSimplifyed.pdf
#            3.Output/CellFigure/funcEnrichDataset.csv
# Upstream : 14_geseca_functional_enrichment.R
# ----------------------------------------------------------------------------

library(msigdbr)
library(Seurat)
library(presto)
library(tidyverse)
library(fgsea)
library(fastmatch)
library(ComplexHeatmap)
library(circlize)
library(BBmisc)
library(extrafont)

obj <- readRDS("1.Data/datForRevision.RDS")

es <- AverageExpression(obj, assays = 'RNA')[[1]]
colnames(es) <- paste0('C', colnames(es))

obj <- SCTransform(obj, verbose = FALSE, variable.features.n = 10000)

length(VariableFeatures(obj)) # make sure it's a full gene universe of 10000 genes
#> [1] 10000
obj <- RunPCA(obj, assay = "SCT", verbose = FALSE,
              rev.pca = TRUE, reduction.name = "pca.rev",
              reduction.key="PCR_", npcs = 50)

E <- obj@reductions$pca@feature.loadings



# KEGG 

m_df <- msigdbr(species = "Mus musculus", 
                category = "C2",
                subcategory = 'KEGG')

pathways <- split(m_df$gene_symbol, m_df$gs_name)

gesecaRes <- geseca(pathways, E, minSize = 5, maxSize = 500, center = FALSE, eps=1e-100)


# for heatmap
topPathways <- gesecaRes %>% filter(padj < 0.05) %>% pull(pathway)


# for umap 
topPathways <- c("KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION",
                 "KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION",
                 "KEGG_MAPK_SIGNALING_PATHWAY",
                 "KEGG_CELL_CYCLE",
                 "KEGG_T_CELL_RECEPTOR_SIGNALING_PATHWAY",
                 "KEGG_JAK_STAT_SIGNALING_PATHWAY",
                 "KEGG_RIG_I_LIKE_RECEPTOR_SIGNALING_PATHWAY",
                 "KEGG_NOD_LIKE_RECEPTOR_SIGNALING_PATHWAY")
# plotGesecaTable(gesecaRes %>% filter(padj < 0.05), 
#                 pathways, E, titles = E$title)

titles <- gsub("_", ' ', sub("KEGG_", "", topPathways))

trace(plotCoregulationProfileReduction, edit=TRUE)

ps <- plotCoregulationProfileReduction(pathways[topPathways], obj,
                                       title=titles,
                                       pt.size = 1, label = T,
                                       reduction="umap", colors = c("blue", "transparent", "red")) 

png('3.Output/CellFigure/keggHeatmapUMAP.png', width = 7, height = 12, units = 'in', res = 300)

cowplot::plot_grid(plotlist=ps[1:8], ncol=2) 

dev.off()

# plotting 
topPathwaySub <- topPathways[c(1:3, 5, 8, 9, 23, 28, 32, 58)]
pathwaysP <- pathways[topPathwaySub]


E <- t(scale(t(es), scale = FALSE))

pathwaysP <- lapply(pathwaysP, function(p) {
  unname(as.vector(na.omit(fmatch(p, rownames(E)))))
})

prjs <- t(do.call(cbind, lapply(pathwaysP, function(p) {
  colSums(E[p, ])/(length(p))
})))

ctMat1 <- normalize(as.matrix(prjs), range = c(-3, 3))
colnames(ctMat1) <- colnames(prjs)
rownames(ctMat1) <- gsub("KEGG", '', rownames(ctMat1))
rownames(ctMat1) <- gsub("_", ' ', rownames(ctMat1))

col_fun = circlize::colorRamp2(c(-3, 0, 3), c("blue", "white", "red"))

png('3.Output/CellFigure/keggHeatmapSub.png', width = 8, height = 3, units = 'in', res = 300) # 12.5
Heatmap(ctMat1,
        col = col_fun,
        cluster_rows = F,
        cluster_columns = F,
        row_split = 1:nrow(ctMat1), 
        row_title = rownames(ctMat1),
        row_title_rot = 0,
        row_title_gp = grid::gpar(fontfamily ='Arial', fontface = 2, fontsize = 12),
        # row_names_side = "left",
        row_labels = format(gesecaRes[gesecaRes$pathway %in% topPathwaySub,]$padj, 
                            digits = 3,
                            scientific = TRUE),
        row_names_gp = grid::gpar(fontfamily ='Arial', fontface = 2, fontsize = 8),
        # row_labels_side = 'right',
        rect_gp = gpar(col = "black", lwd = 1),
        width = ncol(ctMat1)*unit(6, "mm"), 
        height = nrow(ctMat1)*unit(6, "mm"),
        #column_names_rot = 0,
        column_names_gp = grid::gpar(fontfamily ='Arial', fontface = 2, fontsize = 10), 
        
        show_heatmap_legend = FALSE)
dev.off()

# draw legend 
col_fun = colorRamp2(c(-3, 0, 3), c("blue", "white", "red"))
lgd = Legend(col_fun = col_fun, title = "Z-score", direction = "horizontal", legend_width = unit(6, "cm"))

png('3.Output/CellFigure/keggHeatmapLegend.png', width = 5, height = 1, units = 'in', res = 300)

plot.new()
draw(lgd, x = unit(0.5, "npc"), y = unit(0.5, "npc"))

dev.off()


keggDf <- data.frame()
for (i in names(ps)) {
  df <- ps[[i]]$data %>% 
    add_column(pathwayId = i)
  
  keggDf <- rbind(keggDf, df)
}

plotGesecaTable(gesecaRes |> head(10), pathways, E=exprs(es), titles = es$title)

sumDf <- keggDf %>% 
  mutate(posV = pathway > 0) %>% 
  group_by(ident, pathwayId) %>% 
  mutate(ration = sum(posV)/n()) %>% 
  select(ident, pathwayId, ration) %>% 
  unique() %>% 
  filter(ration > 0.5)

sumDf %>% ggplot(aes(x = ident, y = pathwayId)) +
  geom_point(aes(size = ration))


cowplot::plot_grid(plotlist=ps[1:10], ncol=2)

# 

library(GEOquery)
library(limma)


gse200250 <- getGEO("GSE200250", AnnotGPL = TRUE)[[1]]

es <- gse200250
es <- es[, grep("Th2_", es$title)]
es$time <- as.numeric(gsub(" hours", "", es$`time point:ch1`))
es <- es[, order(es$time)]

exprs(es) <- normalizeBetweenArrays(log2(exprs(es)), method="quantile")

es <- es[order(rowMeans(exprs(es)), decreasing=TRUE), ]
es <- es[!duplicated(fData(es)$`Gene ID`), ]
rownames(es) <- fData(es)$`Gene ID`
es <- es[!grepl("///", rownames(es)), ]
es <- es[rownames(es) != "", ]

fData(es) <- fData(es)[, c("ID", "Gene ID", "Gene symbol")]

es <- es[head(order(rowMeans(exprs(es)), decreasing=TRUE), 12000), ]
head(exprs(es))

es$`time point:ch1`
es$title

library(msigdbr)
pathwaysDF <- msigdbr("mouse", category="H")
pathways <- split(as.character(pathwaysDF$entrez_gene), pathwaysDF$gs_name)

es <- AverageExpression(obj, assays = 'RNA')[[1]]
colnames(es) <- paste0('C', colnames(es))

gesecaRes <- geseca(pathways, es, minSize = 15, maxSize = 500)

plotGesecaTable(gesecaRes %>% filter(padj < 0.05), pathways, E=es)
