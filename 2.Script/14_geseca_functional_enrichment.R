# ----------------------------------------------------------------------------
# Fig 2 | Gene set co-regulation analysis (GESECA) across clusters
#
# Paper    : LAG-3 and PD-1 synergize on CD8+ T cells to drive T cell
#            exhaustion and hinder autocrine IFN-gamma-dependent anti-tumor
#            immunity
# Original : 2.Script/rev4_funcEnrich.R
# Input    : 1.Data/datForRevision.RDS
# Output   : 3.Output/Revision/funcEnrich.RDS
#            3.Output/CellFigure/keggHeatmap*.png
# Upstream : 04_cluster_annotation.R
# Note     : fgsea 1.28.0 GESECA
#            KEGG gene sets from msigdbr
# ----------------------------------------------------------------------------

library(Seurat)
library(presto)
library(fgsea)
library(msigdbr)
library(tidyverse)
library(simplifyEnrichment)
library(clusterProfiler)
library(org.Mm.eg.db)
library(export)

dat <- readRDS("1.Data/datForRevision.RDS")

Idents(dat) <- 'anno'

deG <- wilcoxauc(dat)

mkG <- cosg(dat)

goList <- list()

for (i in 1:length(unique(deG$group))) {
  gList <- deG %>% 
    filter(group == unique(deG$group)[i]) %>% 
    arrange(desc(logFC)) %>% 
    # filter(logFC != 0) %>% 
    dplyr::select(feature, logFC) %>% 
    deframe()
  
  ego <- gseGO(geneList     = gList,
               OrgDb        = org.Mm.eg.db,
               ont          = "BP",
               keyType      = "SYMBOL",
               minGSSize    = 3,
               maxGSSize    = 500,
               pvalueCutoff = 0.05,
               verbose      = FALSE)
  
  # ego <- enrichplot::pairwise_termsim(ego)
  
  goList[[i]] = ego@result
  
  names(goList)[i] = unique(deG$group)[i]
  
}

saveRDS(goList, '3.Output/Revision/funcEnrich.RDS')

t <- do.call(rbind, goList)
t$cluster <- rep(names(goList), sapply(goList, nrow))
write.csv(t, '3.Output/CellFigure/funcEnrichDataset.csv')

trace(simplifyGOFromMultipleLists, edit=TRUE)

# lt3 <- lapply(goList, function(x) x$ID[x$p.adjust < 0.001])
# pdf('3.Output/Revision/goSimplifyed.pdf', width = 12, height = 6)

lt3 = lapply(goList, function(x) x$ID[x$p.adjust < 0.05])
simplifyGOFromMultipleLists(lt3, 
                            padj_cutoff = 0.05,
                            ont = 'BP',
                            db = 'org.Mm.eg.db')

graph2ppt(file="3.Output/Revision/funcEnrich.pptx", width=15, height=8)

dev.off()

# newFunc(goList, 
#                             padj_cutoff = 0.05,
#                             ont = 'BP',
#                             db = 'org.Mm.eg.db')
# 
# 
# as.list(body(simplifyGOFromMultipleLists))
