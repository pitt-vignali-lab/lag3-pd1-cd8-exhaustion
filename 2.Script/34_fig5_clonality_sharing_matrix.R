# ----------------------------------------------------------------------------
# Fig 5 | Clonotype sharing between clusters as a graph and heat map
#
# Paper    : LAG-3 and PD-1 synergize on CD8+ T cells to drive T cell
#            exhaustion and hinder autocrine IFN-gamma-dependent anti-tumor
#            immunity
# Original : 2.Script/4.7.cellFateClonalitySharMatrx.R
# Input    : 1.Data/datForRevision.RDS
#            1.Data/tcrCombinedUpdated.csv
# Output   : 3.Output/4.TCR/tcrFate_graph_*
#            3.Output/CellFigure/TCR_cmmn_clonality*.png
# Upstream : 29_fig4_tcr_combined.R
# Note     : Network layout with igraph 2.0.1 and ggraph 2.1.0
# ----------------------------------------------------------------------------

library(Seurat)
library(tidyverse)
library(export)

if(!dir.exists("3.Output/4.TCR")){
  dir.create("3.Output/4.TCR")
}

dat <- readRDS("1.Data/datForRevision.RDS")

table(dat$sample)

DimPlot(dat, label = T) +
  coord_fixed()

datNor <- subset(dat, sample == "Normal CD8")

tcrDf <- read.csv("1.Data/tcrCombinedUpdated.csv")[-1] %>%
  # mutate(barcode = paste0("nor_", barcode)) %>% 
  filter(barcode %in% rownames(datNor@meta.data)) 

allAnno <- unique(dat$anno)

allplt <- data.frame()
for (gnty in c('all', unique(dat$genotype))) {
  if(gnty == 'all'){
    datSub = datNor
  } else {
    datSub = subset(datNor, 
                    cells = tcrDf$barcode,
                    genotype == gnty)
  }
  
  tcrDf1 <- data.frame(type = datSub$anno) %>% 
    rownames_to_column('barcode') %>% 
    merge(tcrDf, by = 'barcode') %>% 
    select(-barcode)
  
  plotDf <- data.frame()
  for (anno1 in unique(tcrDf1$type)) {
    for (anno2 in unique(tcrDf1$type)[unique(tcrDf1$type) != anno1]) {
      temDf = data.frame(anno1 = anno1,
                         anno2 = anno2,
                         overlap = length(intersect(tcrDf1[tcrDf1$type == anno1,]$cdr3,
                                                    tcrDf1[tcrDf1$type == anno2,]$cdr3)))
      plotDf <- rbind(plotDf, temDf)
    }
  }
  
  plotDf[c(1, 2)] <- t(apply(plotDf[c(1, 2)], 1, 
                             function(x) sort(t(x))))
  
  plotDf <- unique(plotDf) %>% 
    mutate(type = gnty)
  
  allplt <- rbind(allplt, plotDf)
}

rbind(data.frame(anno1 = 'C4: Exhuasted', 
           anno2 = allAnno[c(2, 3, 7, 8)], overlap = 0, type = "Pdcd1L/L E8ICre.GFP"),
      data.frame(anno2 = 'C4: Exhuasted', 
                 anno1 = allAnno[c(5, 1, 4, 6)], overlap = 0, type = "Pdcd1L/L E8ICre.GFP"),
allplt) %>%
  mutate(type = factor(type, levels = c('all', "E8ICre.GFP (WT)", 
                                        "Pdcd1L/L Lag3L/L-yfp (WT)", 
                                        "Lag3L/L-yfp E8ICre.GFP", 
                                        "Pdcd1L/L E8ICre.GFP", 
                                        "Pdcd1L/L Lag3L/L-yfp E8ICre.GFP"))) %>% 
  ggplot(aes(anno2, anno1, fill = overlap)) +
  geom_tile(col = 'black') +
  scale_fill_gradient(low = 'white', high = 'orange') +
  geom_text(aes(label = overlap), size = 4.5, fontface = 2) +
  facet_grid(~type) +
  coord_fixed() +
  scale_x_discrete(limits=rev) +
  labs(x = '', y = '',
       fill =  'No. of common \nTCR clonality') +
  theme(panel.background = element_blank(),
        # strip.background = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(face = 2))

graph2ppt(file = '3.Output/CellFigure/TCR_cmmn_clonality.pptx', scaling = 200)
ggsave('3.Output/CellFigure/TCR_cmmn_clonality_Labeled.png', width = 15, height = 4, units = 'in')

rbind(data.frame(anno1 = 'C4: Exhuasted', 
                 anno2 = allAnno[c(2, 3, 7, 8)], overlap = 0, type = "Pdcd1L/L E8ICre.GFP"),
      data.frame(anno2 = 'C4: Exhuasted', 
                 anno1 = allAnno[c(5, 1, 4, 6)], overlap = 0, type = "Pdcd1L/L E8ICre.GFP"),
      allplt) %>%
  mutate(type = factor(type, levels = c('all', "E8ICre.GFP (WT)", 
                                        "Pdcd1L/L Lag3L/L-yfp (WT)", 
                                        "Lag3L/L-yfp E8ICre.GFP", 
                                        "Pdcd1L/L E8ICre.GFP", 
                                        "Pdcd1L/L Lag3L/L-yfp E8ICre.GFP"))) %>% 
  ggplot(aes(anno2, anno1, fill = overlap)) +
  geom_tile(col = 'black') +
  scale_fill_gradient(low = 'white', high = 'orange') +
  geom_text(aes(label = overlap), size = 4.5, fontface = 2) +
  facet_grid(~type) +
  coord_fixed() +
  scale_x_discrete(limits=rev) +
  labs(x = '', y = '',
       fill =  'No. of common \nTCR clonality') +
  theme(panel.background = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        text = element_text(face = 2))

ggsave('3.Output/CellFigure/TCR_cmmn_clonality_NoLabeled.png', width = 15, height = 3, units = 'in')

ggsave(paste0('3.Output/4.TCR/', gsub('/', '_', gnty), 'cmmn_clonitiy.png'),
       width = 6, height = 5)





# version 1 ----
for (type in c('all', unique(datNor$types))) {
  if (type == 'all'){
    datSub <- subset(datNor, cells = tcrDf$barcode)
  } else {
    datSub <- subset(datNor, cells = tcrDf$barcode,
                     types == type)
  }
  
  
  
  # grid heatmap ----
  clstDf <- data.frame(Embeddings(datSub, 'umap'),
                       Idents(datSub)) %>% 
    merge(tcrDf, by.x = 0, by.y = 'barcode') %>% 
    rename(clust = 4) 
  
  # DimPlot(datSub, label = T)
  # 
  # ggplot(clstDf, aes(UMAP_1, UMAP_2, col = paste0('c', clust))) +
  #   geom_point()
  # umapDf <- merge(Embeddings(datSub, 'umap'), clstDf,
  #                 by.x = 0, by.y = 'Row.names')
  
  
  cmnTCR <- data.frame()
  
  for (i in 0:7) {
    for (j in (i+1):8) {
      sharedTCR <- intersect(clstDf %>% 
                               filter(clust == i) %>% 
                               pull(cdr3),
                             clstDf %>% 
                               filter(clust == j) %>% 
                               pull(cdr3))
      
      temDf <- data.frame(clst1 = i, clst2 = j, 
                          shared = length(sharedTCR))
      
      cmnTCR <- rbind(cmnTCR, temDf)
      
    }
  }
  
  # sig test ----
  
  for (sigi in unique(cmnTCR$shared)) {
    smp1 <- cmnTCR %>% 
      filter(shared < sigi) %>% 
      pull(shared)
    
    smp2 <- cmnTCR %>% 
      filter(shared >= sigi) %>% 
      pull(shared)
    
    if (wilcox.test(smp1, smp2, exact=FALSE)$p.value < 0.05) break
    cat(sigi)
  }
  
  boxplot(cmnTCR$shared)
  
  # plotting ----
  cmnTCR <- cmnTCR %>% 
    mutate(sig = ifelse(shared >= sigi, 'Sig-Hi', 'Low'))
  
  ggplot(cmnTCR, aes(x = paste0('C', clst1), 
                     y = paste0('C', clst2))) +
    geom_raster(aes(fill = shared)) + 
    geom_tile(fill = '#00000000', col = 'gray')+
    geom_tile(fill = '#00000000', 
              aes(col = sig, size = sig)) +
    scale_fill_gradient(low = "white",  
                        high = "orange",
                        breaks = seq(min(cmnTCR$shared), 
                                     max(cmnTCR$shared), 20))+ 
    scale_color_manual(values = c('#00000000', 'black')) +
    scale_size_manual("sig", values = c(1, 1)) +
    scale_y_discrete(limits=rev) +
    geom_text(aes(label = shared), size = 3) +
    coord_fixed() +
    labs(x = '', y = '',
         fill =  'No. of common \nTCR clonalities') +
    theme(panel.background = element_blank(),
          axis.ticks = element_blank())
  
  ggsave(paste0('3.Output/4.TCR/tcrFate_htmp_', gsub("/", "_", type), '.png'), 
         width = 5, height = 5, units = 'in')
  ggsave(paste0('3.Output/4.TCR/tcrFate_htmp_', gsub("/", "_", type), '.pdf'), 
         width = 5, height = 5)
  
  # graph ----

  
  DimPlot(datSub, label = T)
  
  # ggplot(umapDf, aes(UMAP_1, UMAP_2, col = paste('c', clust))) +
  #   geom_point()
  
  nodeDf <- clstDf %>% 
    group_by(clust) %>% 
    summarise(umap1 = median(UMAP_1),
              umap2 = median(UMAP_2)) %>% 
    merge(data.frame(table(clstDf$clust)), 
          by.x = 'clust', by.y = 'Var1') 
  
  edgeDf <- cmnTCR %>% 
    filter(shared != 0) %>% 
    merge(nodeDf[-4], by.x = 'clst1', by.y = 'clust') %>% 
    merge(nodeDf[-4], by.x = 'clst2', by.y = 'clust')
  
  
  edgeDfSub <- edgeDf # %>% filter(clst1 == 6 | clst2 == 6)
  
  ggplot() +
    # geom_point(data = umapDf, aes(UMAP_1, UMAP_2, fill = as.character(clust)),
    #            alpha = 0.1, shape = 21, col = '#00000000', size = 10) +
    geom_curve(data = edgeDfSub, 
               aes(x = umap1.x, xend = umap1.y, 
                   y = umap2.x, yend = umap2.y,
                   col =  sig),
               linewidth = ifelse(edgeDfSub$sig == 'Sig-Hi',
                                  edgeDfSub$shared/20, 0.5),
               linetype = ifelse(edgeDfSub$sig == 'Sig-Hi', 
                                 'solid', 'dashed'),
               curvature = 0) +
    scale_color_manual(values = c('gray', 'darkred')) +
    geom_point(data = nodeDf, 
               aes(umap1, umap2,
                   fill = as.character(clust)),
               shape = 21, col = '#00000000',
               size = log10(nodeDf$Freq) * 8,
               show.legend = F) +
    geom_text(data = nodeDf, aes(umap1, umap2,
                                 label = paste0("C", clust)),
              size = 4) +
    xlim(-8, 8) +
    ylim(-8, 8) +
    coord_fixed() +
    theme_void()
  
  ggsave(paste0('3.Output/4.TCR/tcrFate_graph_', gsub("/", "_", type), '.png'), 
         width = 5, height = 5, units = 'in', bg = 'white')
  ggsave(paste0('3.Output/4.TCR/tcrFate_graph_', gsub("/", "_", type), '.pdf'), 
         width = 5, height = 5)
  
  cat(type, "is done! \n")
}

