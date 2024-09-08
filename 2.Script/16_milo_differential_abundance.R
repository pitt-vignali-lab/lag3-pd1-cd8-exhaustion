# ----------------------------------------------------------------------------
# Fig 2 | Differential abundance testing on the k-NN graph
#
# Paper    : LAG-3 and PD-1 synergize on CD8+ T cells to drive T cell
#            exhaustion and hinder autocrine IFN-gamma-dependent anti-tumor
#            immunity
# Original : 2.Script/rev2_milo.R
# Input    : 1.Data/datForRevision.RDS
# Output   : 1.Data/milorData.RData
#            3.Output/Revision/milo*.png
#            3.Output/CellFigure/Milo_*.png
# Upstream : 04_cluster_annotation.R
# Note     : miloR 1.10.0 with SingleCellExperiment, scater and scran
# ----------------------------------------------------------------------------

library(Seurat)
library(tidyverse)
library(miloR)
library(SingleCellExperiment)
library(scater)
library(scran)
library(dplyr)
library(patchwork)
library(igraph)
library(ggraph)
library(rts)
library(export)

dat <- readRDS('1.Data/datForRevision.RDS')

DimPlot(dat, group.by = 'anno', label = T) +
  coord_fixed()

dat$GenoAnno <- paste(dat$genotype, dat$anno, sep = '||')

head(dat$GenoAnno)
# overall milor ----
datSce <- as.SingleCellExperiment(dat)

milo <- Milo(datSce)

# reducedDim(milo, "UMAP") <- reducedDim(datSce, "UMAP")

# milo

milo <- buildGraph(milo, k = 10, d = 30)

milo <- makeNhoods(milo, prop = 0.1, 
                   k = 10, 
                   d = 30, 
                   refined = TRUE)

plotNhoodSizeHist(milo)

milo <- countCells(milo, 
                   meta.data = data.frame(colData(milo)), 
                   sample="GenoAnno")

head(nhoodCounts(milo))

designDf <- data.frame(colData(milo))[,c('genotype', 'anno', "GenoAnno")] %>% 
  distinct() %>% 
  remove_rownames() %>% 
  # mutate(tissueType = paste(tissue, anno, sep = '_')) %>% 
  column_to_rownames(var = "GenoAnno") 

# designDf <- designDf[colnames(nhoodCounts(milo)), , drop=FALSE]

# designDf$sampGenoAnno <- rownames(designDf)

designDf

milo <- calcNhoodDistance(milo, d = 30)

# DA test 
da_results <- testNhoods(milo, 
                         design = ~ genotype, 
                         design.df = designDf)
head(da_results)

da_results %>%
  arrange(SpatialFDR) %>%
  head() 

ggplot(da_results, aes(PValue)) + geom_histogram(bins=50)

ggplot(da_results, aes(logFC, -log10(SpatialFDR))) + 
  geom_point() +
  geom_hline(yintercept = -log(0.05)) ## Mark significance threshold (10% FDR)

milo <- buildNhoodGraph(milo)

load('1.Data/milorData.RData')

# overall graph 
plotNhoodGraphDA(milo, 
                 da_results, 
                 layout="UMAP",
                 alpha=0.1) +
  scale_fill_gradient2(high = 'red', mid = 'white', low = 'blue', midpoint = 0) +
  labs(fill = 'logFC') +
  coord_fixed() +
  labs(x = '', 
       y = '') +
  guides(edge_width = 'none', size = 'none')+ 
  ggforce::theme_no_axes() +
  theme(text = element_text(face = 2, size = 15), legend.position = 'bottom')

ggsave('3.Output/CellFigure/Milo_Overall.png', width = 5, height = 5.5, units = 'in')
# ggsave('3.Output/Revision/milo.png', width = 6, height = 6, units = 'in')

# seperate graph 

da_results <- annotateNhoods(milo, 
                             da_results, 
                             coldata_col = "genotype")
nh_graph <- nhoodGraph(milo)

eDf <- data.frame()
vDf <- data.frame()

for (i in unique(da_results$genotype)) {
  signif_res <- da_results %>% 
    filter(genotype == i,
           SpatialFDR < 0.1)
  if(nrow(signif_res) > 0){
    signif_res$name <- unlist(nhoodIndex(milo)[signif_res$Nhood])
    
    nh_graphSub <- igraph::induced_subgraph(nh_graph,
                                            vids = which(as.numeric(V(nh_graph)$name) %in% signif_res$name))
    
    eDfTemp <- as_data_frame(nh_graphSub) %>% 
      filter(from != to) %>% 
      mutate(cat = i)
    
    vDfTemp <- data.frame(name = V(nh_graphSub)$name,
                          size = V(nh_graphSub)$size) %>% 
      mutate(vcat = i)
    
    vDfTemp <- merge(signif_res[c('name', 'logFC')], vDfTemp, by = 'name')
    
    vDf = rbind(vDf, vDfTemp)
    eDf = rbind(eDf, eDfTemp)
  }
  
}

sub_graph <- graph_from_data_frame(eDf, vertices = vDf)

sub_graph <- permute(sub_graph, order(vertex_attr(sub_graph)$size,
                                      decreasing = TRUE))

layout <- reducedDim(milo, 'UMAP')[as.numeric(vertex_attr(sub_graph)$name), ]

bg <- reducedDim(milo, 'UMAP')[as.numeric(vertex_attr(nh_graph)$name), ] %>% 
  as.data.frame() %>% 
  rownames_to_column('cells') %>% 
  merge(data.frame(genotype = dat$genotype,
                   anno = dat$anno) %>% 
          # filter(time != 'baseline') %>% 
          rownames_to_column('cells'),
        by = 'cells') 

V(sub_graph)$size = V(sub_graph)$size/4

ggraph(sub_graph, 
       layout = layout) + 
  # geom_point(data = bg,
  #            aes(UMAP_1, UMAP_2, col = anno),
  #            alpha = 0.5,
  #            size = 2,
  #            show.legend = F) +
  facet_grid(.~genotype) +
  geom_edge_link0(aes(width = weight), 
                  edge_colour = "grey50", edge_alpha = 0.2) + 
  geom_node_point(aes(fill = logFC, 
                      size = size), 
                  shape = 21, 
                  stroke = 0.3) + 
  scale_size(range = c(1, 5), name = "Nhood size") + 
  scale_edge_width(range = c(1, 5), name = "overlap size") +
  facet_nodes(factor(vcat, 
                     levels = c("Pdcd1L/L Lag3L/L-yfp (WT)",
                                "E8ICre.GFP (WT)",
                                "Lag3L/L-yfp E8ICre.GFP",
                                "Pdcd1L/L E8ICre.GFP",
                                "Pdcd1L/L Lag3L/L-yfp E8ICre.GFP"))~., 
              nrow = 1, strip.position = "top") + #, margins = TRUE
  scale_fill_gradient2(low = 'blue', high = 'red', mid = 'white', midpoint = 0) +
  coord_fixed() +
  theme_bw() +
  labs(x = NULL, y = NULL) + 
  guides(edge_width = 'none', size = 'none')+ 
  ggforce::theme_no_axes() +
  theme(text = element_text(size = 12, face = 2),
        strip.background = element_blank(),
        strip.text = element_blank(), # element_text(size = 15, face = 2),
        legend.position = 'bottom')

ggsave('3.Output/CellFigure/Milo_split_NoLabel.png', width = 20, height = 5, units = 'in')
# ggsave('3.Output/Revision/miloByGenotypeNoColor.png', width = 16, height = 6, units = 'in')


# double deficient only with highlight 
da_results <- annotateNhoods(milo, 
                             da_results, 
                             coldata_col = "genotype")
nh_graph <- nhoodGraph(milo)

eDf <- data.frame()
vDf <- data.frame()

i <- "Pdcd1L/L Lag3L/L-yfp E8ICre.GFP"
signif_res <- da_results %>% 
  filter(genotype == i,
         SpatialFDR < 0.1)
if(nrow(signif_res) > 0){
  signif_res$name <- unlist(nhoodIndex(milo)[signif_res$Nhood])
  
  nh_graphSub <- igraph::induced_subgraph(nh_graph,
                                          vids = which(as.numeric(V(nh_graph)$name) %in% signif_res$name))
  
  eDfTemp <- as_data_frame(nh_graphSub) %>% 
    filter(from != to) %>% 
    mutate(cat = i)
  
  vDfTemp <- data.frame(name = V(nh_graphSub)$name,
                        size = V(nh_graphSub)$size) %>% 
    mutate(vcat = i)
  
  vDfTemp <- merge(signif_res[c('name', 'logFC')], vDfTemp, by = 'name')
  
  vDf = rbind(vDf, vDfTemp)
  eDf = rbind(eDf, eDfTemp)
}


dk_graph <- graph_from_data_frame(eDf, vertices = vDf)


dk_graph <- permute(dk_graph, order(vertex_attr(dk_graph)$size,
                                    decreasing = TRUE))

layout <- reducedDim(milo, 'UMAP')[as.numeric(vertex_attr(dk_graph)$name), ]


density_quantiles <- function(x, y, quantiles) {
  dens <- MASS::kde2d(x, y, n = 500)
  df   <- cbind(expand.grid(x = dens$x, y = dens$y), z = c(dens$z))
  r    <- terra::rast(df)
  ind  <- sapply(seq_along(x), function(i) cellFromXY(r, cbind(x[i], y[i])))
  ind  <- ind[order(-r[ind][[1]])]
  vals <- r[ind][[1]]
  ret  <- approx(seq_along(ind)/length(ind), vals, xout = quantiles)$y
  replace(ret, is.na(ret), max(r[]))
}


ggraph(dk_graph, 
       layout = layout) +
  geom_edge_link0(aes(width = weight), 
                  edge_colour = "grey50", 
                  edge_alpha = 0.2) + 
  geom_node_point(aes(fill = logFC, 
                      size = size/15), 
                  shape = 21, 
                  stroke = 0.3) + 
  scale_size(range = c(1, 5), name = "Nhood size") + 
  scale_edge_width(range = c(1, 5), name = "overlap size") + #, margins = TRUE
  scale_fill_gradient2(low = 'blue', high = 'red', mid = 'white', midpoint = 0) +
  geom_density_2d(data = bg %>% 
                          filter(anno %in% c('C7: IFN-induced', 'C0: Effector/effector memory-like')),
                        aes(UMAP_1, UMAP_2, col = 'Effector State'),
                        breaks = density_quantiles(bg %>% 
                                                     filter(anno %in% 
                                                              c('C7: IFN-induced', 'C0: Effector/effector memory-like')) %>% 
                                                     pull(UMAP_1), 
                                                   bg %>% 
                                                     filter(anno %in% 
                                                              c('C7: IFN-induced', 'C0: Effector/effector memory-like')) %>% 
                                                     pull(UMAP_2), 
                                                   c(0, 0.25, 0.75, 1)),
                  linewidth = 1, lty = 5) +
  geom_density_2d(data = bg %>% 
                    filter(anno %in% c('C1: Activated', 'C3: Proliferation')),
                  aes(UMAP_1, UMAP_2, col = 'Activated State'),
                  breaks = density_quantiles(bg %>% 
                                               filter(anno %in% 
                                                        c('C1: Activated', 'C3: Proliferation')) %>% 
                                               pull(UMAP_1), 
                                             bg %>% 
                                               filter(anno %in% 
                                                        c('C1: Activated', 'C3: Proliferation')) %>% 
                                               pull(UMAP_2), 
                                             c(0, 0.4, 0.69)),
                  linewidth = 1, lty = 5)+
  geom_density_2d(data = bg %>% 
                    filter(anno %in% c('C4: Exhuasted')),
                  aes(UMAP_1, UMAP_2, col = 'Exhausted State'),
                  breaks = density_quantiles(bg %>% 
                                               filter(anno %in% 
                                                        c('C4: Exhuasted')) %>% 
                                               pull(UMAP_1), 
                                             bg %>% 
                                               filter(anno %in% 
                                                        c('C4: Exhuasted')) %>% 
                                               pull(UMAP_2), 
                                             c(0, 0.40, 0.9)),
                  linewidth = 1, lty = 5) +
  scale_color_manual(breaks = c('Activated State', 'Effector State', 'Exhausted State'),
                     values = rev(viridis::viridis(3, option = "D")),
                     name = 'Cell state') +
  coord_fixed() +
  labs(x = '',
       y = '') + 
  guides(edge_width = 'none', size = 'none')+ 
  ggforce::theme_no_axes() +
  guides(fill = 'none') +
  theme(text = element_text(face = 2, size = 15), 
        legend.position = 'bottom',
        legend.direction = "vertical")

ggsave('3.Output/CellFigure/Milo_DKonly.png', width = 5, height = 6, units = 'in')
# ggsave('3.Output/Revision/miloDKonly.png', width = 7, height = 7, units = 'in')

## Add log normalized count to Milo object
da_results$NhoodGroup <- da_results$genotype 

da_results <- da_results %>% 
  filter(SpatialFDR < 0.05)

da_nhood_markers <- findNhoodGroupMarkers(milo, 
                                          da_results, 
                                          subset.row = rownames(milo))

write.csv(da_nhood_markers, '3.Output/Revision/mkGsOfSigMiloGroup.csv')

  
save.image(file = '1.Data/milorData.RData')




  #  pMEL data milor ----
datSub <- subset(dat, sample == 'PML CD8') # PML CD8

# dat$sampGeno <- paste(dat$sample, dat$genotype)
datSub$GenoAnno <- paste(datSub$genotype, datSub$anno, sep = '||')

# dat <- subset(dat, sample == 'Normal CD8')

datSce <- as.SingleCellExperiment(datSub)

miloPmel <- Milo(datSce)

# reducedDim(milo, "UMAP") <- reducedDim(datSce, "UMAP")

# milo

miloPmel <- buildGraph(miloPmel, k = 10, d = 30)

miloPmel <- makeNhoods(miloPmel, prop = 0.1, 
                   k = 10, 
                   d = 30, 
                   refined = TRUE)

plotNhoodSizeHist(miloPmel)

miloPmel <- countCells(miloPmel, 
                   meta.data = data.frame(colData(miloPmel)), 
                   sample="GenoAnno")

head(nhoodCounts(miloPmel))

designDf <- data.frame(colData(miloPmel))[,c('genotype', 'anno', "GenoAnno")] %>% 
  distinct() %>% 
  remove_rownames() %>% 
  # mutate(tissueType = paste(tissue, anno, sep = '_')) %>% 
  column_to_rownames(var = "GenoAnno") 

# designDf <- designDf[colnames(nhoodCounts(milo)), , drop=FALSE]

# designDf$sampGenoAnno <- rownames(designDf)

designDf

miloPmel <- calcNhoodDistance(miloPmel, d = 30)

# DA test 
da_results <- testNhoods(miloPmel, 
                         design = ~ genotype, 
                         design.df = designDf)
head(da_results)

da_results %>%
  arrange(SpatialFDR) %>%
  head() 

ggplot(da_results, aes(PValue)) + geom_histogram(bins=50)

ggplot(da_results, aes(logFC, -log10(SpatialFDR))) + 
  geom_point() +
  geom_hline(yintercept = 1) ## Mark significance threshold (10% FDR)

miloPmel <- buildNhoodGraph(miloPmel)

# overall graph 
plotNhoodGraphDA(miloPmel, 
                 da_results, 
                 layout="UMAP",
                 alpha=0.05) +
  scale_fill_gradient2(high = 'red', mid = 'white', low = 'blue', midpoint = 0) +
  labs(fill = 'logFC') +
  coord_fixed() +
  theme_bw()

ggsave('3.Output/Revision/miloPmel.png', width = 6, height = 6, units = 'in')

# seperate graph 

da_results <- annotateNhoods(miloPmel, 
                             da_results, 
                             coldata_col = "genotype")
nh_graph <- nhoodGraph(miloPmel)

eDf <- data.frame()
vDf <- data.frame()

for (i in unique(da_results$genotype)) {
  signif_res <- da_results %>% 
    filter(genotype == i,
           SpatialFDR < 0.05)
  if(nrow(signif_res) > 0){
    signif_res$name <- unlist(nhoodIndex(miloPmel)[signif_res$Nhood])
    
    nh_graphSub <- igraph::induced_subgraph(nh_graph,
                                            vids = which(as.numeric(V(nh_graph)$name) %in% signif_res$name))
    
    eDfTemp <- as_data_frame(nh_graphSub) %>% 
      filter(from != to) %>% 
      mutate(cat = i)
    
    vDfTemp <- data.frame(name = V(nh_graphSub)$name,
                          size = V(nh_graphSub)$size) %>% 
      mutate(vcat = i)
    
    vDfTemp <- merge(signif_res[c('name', 'logFC')], vDfTemp, by = 'name')
    
    vDf = rbind(vDf, vDfTemp)
    eDf = rbind(eDf, eDfTemp)
  }
  
}

sub_graph <- graph_from_data_frame(eDf, vertices = vDf)

sub_graph <- permute(sub_graph, order(vertex_attr(sub_graph)$size,
                                      decreasing = TRUE))

layout <- reducedDim(miloPmel, 'UMAP')[as.numeric(vertex_attr(sub_graph)$name), ]

bg <- reducedDim(miloPmel, 'UMAP')[as.numeric(vertex_attr(nh_graph)$name), ] %>% 
  as.data.frame() %>% 
  rownames_to_column('cells') %>% 
  merge(data.frame(genotype = dat$genotype,
                   anno = dat$anno) %>% 
          # filter(time != 'baseline') %>% 
          rownames_to_column('cells'),
        by = 'cells') 

ggraph(sub_graph, 
       layout = layout) + 
  # geom_point(data = bg,
  #            aes(UMAP_1, UMAP_2, col = anno),
  #            alpha = 0.5,
  #            size = 2,
  #            show.legend = F) +
  facet_grid(.~genotype) +
  geom_edge_link0(aes(width = weight), 
                  edge_colour = "grey50", edge_alpha = 0.2) + 
  geom_node_point(aes(fill = logFC, 
                      size = size/15), 
                  shape = 21, 
                  stroke = 0.3) + 
  scale_size(range = c(1, 5), name = "Nhood size") + 
  scale_edge_width(range = c(1, 5), name = "overlap size") +
  facet_nodes(factor(vcat, 
                     levels = c("Pdcd1L/L Lag3L/L-yfp (WT)",
                                "Lag3L/L-yfp E8ICre.GFP",
                                "Pdcd1L/L E8ICre.GFP",
                                "Pdcd1L/L Lag3L/L-yfp E8ICre.GFP"))~., 
              ncol = 1, strip.position = "right") + #, margins = TRUE
  scale_fill_gradient2(low = 'blue', high = 'red', mid = 'white', midpoint = 0) +
  coord_fixed() +
  theme_bw() +
  labs(x = NULL, y = NULL) +
  theme(text = element_text(size = 12, face = 2),
        strip.background = element_blank(),
        strip.text = element_text(size = 5, face = 2),
        legend.position = 'right')

ggsave('3.Output/Revision/miloPmelByGenotypeNoColor.png', width = 6, height = 16, units = 'in')


##############################




plotNhoodGraphDA(miloPmel, 
                 da_results, 
                 layout="UMAP",
                 alpha=0.05) +
  scale_fill_gradient2(high = 'red', mid = 'white', low = 'blue', midpoint = 0) +
  labs(fill = 'logFC')  +
  geom_point(data = data.frame(Embeddings(dat, reduction = 'umap'),
                               anno = dat$anno),
             aes(UMAP_1, UMAP_2, col = anno), alpha = 0.1) +
  coord_fixed()



t <- plotNhoodGraphDA(miloPmel, 
                      da_results, 
                      layout="UMAP",
                      alpha=0.05) +
  scale_fill_gradient2(high = 'red', mid = 'white', low = 'blue', midpoint = 0) +
  labs(fill = 'logFC') +
  coord_fixed() 

t$data %>% head

library(ggisoband)

t +
  # geom_point(data = data.frame(Embeddings(dat, reduction = 'umap'),
  #                              anno = dat$anno),
  #            aes(UMAP_1, UMAP_2, col = anno), alpha = 0.5) +
  geom_density_2d(data = t$data %>% filter(colour_by != 0), 
                  aes(x, y),
                  binwidth = 0.009) +
  coord_fixed() 

ggsave('3.Output/Revision/miloPmelWBins.png', width = 6, height = 6, units = 'in')

densPlot <- geom_density_2d(data = t$data %>% filter(colour_by != 0), 
                            aes(x, y),
                            binwidth = 0.008)

ggplot() + 
  geom_point(data = data.frame(Embeddings(dat, reduction = 'umap'),
                               anno = dat$anno),
             aes(UMAP_1, UMAP_2, col = anno), alpha = 0.5) +
  geom_density_2d(data = t$data %>% filter(colour_by != 0), 
                  aes(x, y),
                  binwidth = 0.009,
                  col = 'black',
                  lty = 2)

# milo <- logNormCounts(milo)

# da_results <- annotateNhoods(milo, 
#                              da_results, 
#                              coldata_col = "GenoAnno")
# 
# da_results %>% 
#   separate(GenoAnno, into = c('genotype', 'celltype'), sep = '\\|\\|') %>% 
#   # mutate(is_signif = ifelse(SpatialFDR < 0.05, 
#   #                           1, 0)) %>% 
#   mutate(logFC_shape = ifelse(SpatialFDR < 0.05,
#                               'sig', 'non-sig')) %>% 
#   # arrange(GenoAnno) %>% 
#   mutate(Nhood = factor(Nhood, 
#                         levels = unique(Nhood))) %>% 
#   ggplot(aes(logFC, genotype, shape = logFC_shape, col = logFC)) + 
#   scale_color_gradient2(high = 'red', mid = 'white', low = 'blue', midpoint = 0) +
#   facet_grid(celltype~.) +
#   guides(color = "none") + 
#   xlab('GenoAnno') + 
#   ylab("Log Fold Change") + 
#   geom_quasirandom(alpha = 1) + 
#   scale_shape_manual(values = c(1, 16)) + 
#   theme(strip.text.y = element_text(angle = 0))
# 
# da_results %>% 
#   separate(GenoAnno, into = c('genotype', 'celltype'), sep = '\\|\\|') %>%  
#   filter(logFC > 0) %>% 
#   group_by(genotype, celltype) %>% 
#   summarise(sigRatio = sum(FDR < 0.05)/n()) %>% 
#   ggplot(aes(y = genotype, x = sigRatio)) +
#   geom_bar(stat = 'identity') +
#   facet_grid(celltype~.)

# da_results$NhoodGroup <- da_results$sample
# 
# da_nhood_markers <- findNhoodGroupMarkers(milo, 
#                                           da_results, 
#                                           subset.row = rownames(milo))
# 
# mk <- da_nhood_markers %>% 
#   filter(`adj.P.Val_PML CD8` < 0.05)

# specific network by condition


## Add log normalized count to Milo object
da_results$NhoodGroup <- da_results$genotype 
da_results <- da_results %>% 
  filter(SpatialFDR < 0.05)

da_nhood_markers <- findNhoodGroupMarkers(milo, 
                                          da_results, 
                                          subset.row = rownames(milo))

write.csv(da_nhood_markers, '3.Output/Revision/mkGsOfSigMiloGroup.csv')
