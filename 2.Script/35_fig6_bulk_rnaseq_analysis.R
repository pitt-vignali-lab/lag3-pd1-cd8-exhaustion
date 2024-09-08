# ----------------------------------------------------------------------------
# Fig 6 | Clustering and density views of the companion CD8 TIL dataset
#
# Paper    : LAG-3 and PD-1 synergize on CD8+ T cells to drive T cell
#            exhaustion and hinder autocrine IFN-gamma-dependent anti-tumor
#            immunity
# Original : 2.Script/rev8_tonyData.R
# Input    : 1.Data/18071_cd8_seurat_object_231116.rds
#            1.Data/18071_CD8_TIL_DEGs_by_treatment_220425.xlsx
# Output   : 3.Output/CellFigure/tony*.png
# ----------------------------------------------------------------------------

library(Seurat)
library(tidyverse)
library(COSG)
library(ggpubr)

dat <- readRDS('1.Data/18071_cd8_seurat_object_231116.rds')

dat <- subset(dat, treatment != 'off_study' & 
                timepoint %in% c('Baseline', 'Week4') & 
                sample_type == 'TIL')

Idents(dat) <- 'leiden_0.5'

DimPlot(dat, group.by = 'leiden_0.5', label = T) +
  coord_fixed()

# Pull the gene expression data ----
datSub <- subset(dat, idents =  0, sample_type == 'TIL' )

plotDf <- AggregateExpression(datSub, group.by = c('patient_id', 'timepoint', 'treatment'), 
                              features = c('IFNG', 'IFNGR1', 'CD38', 'CCL5', 'HLA-DRA', 'HLA-DRB1'),
                              assays = 'RNA')[[1]] %>% 
  as.matrix() %>% 
  t %>% 
  as.data.frame()

plotDf %>% 
  rownames_to_column('pid_tim') %>% 
  separate(pid_tim, into = c('pid', 'time', 'treatment'), sep = '_') %>% 
  filter(time %in% c('Baseline', 'Week4')) %>% 
  write.csv('3.Output/CellFigure/HumanGeneExpPanel.csv')
# ggplot(aes(x = time, y=V1, group = pid)) +
#   geom_boxplot(aes(group = time)) +
#   geom_point() +
#   geom_line()


# marker genes ----

annot <- data.frame(id = unique(Idents(dat)), 
                    anno = factor(c('Clinical Cluster 4', 
                             'Clinical Cluster 5', 
                             'Clinical Cluster 1', 
                             'Clinical Cluster 2', 
                             'Clinical Cluster 3', 
                             'Clinical Cluster 9', 
                             'Clinical Cluster 7', 
                             'Clinical Cluster 8', 
                             'Clinical Cluster 3', 
                             'Clinical Cluster 6'),
                             levels = paste('Clinical Cluster', 1:9, sep = ' ')))


datNew <- AddMetaData(dat, 
                      data.frame(id = Idents(dat)) %>% 
                        rownames_to_column('cells') %>% 
                        merge(annot, by = 'id') %>%
                        select(-id) %>% 
                        column_to_rownames('cells'))

Idents(datNew) <- 'anno'

DimPlot(datNew, label = T) +
  coord_fixed() +
  NoLegend()

datNew$timeTreat <- paste(datNew$timepoint, datNew$treatment, sep = '_')

condDf <- data.frame(treat = datNew$treatment,
           time = datNew$timepoint)

plotDf <- data.frame()
for (i in c('KLRC1', 'KLRD1', 'KLRC2', 'KLRC3')) {
  tempDf <- FeaturePlot(datNew, features = i)$data %>% 
    merge(condDf, by = 0) %>% 
    rename('exp' = 5) %>% 
    mutate(feat = i)
  
  plotDf = rbind(plotDf, tempDf)
  
}

labDf <- plotDf %>% 
  group_by(ident) %>% 
  summarise(UMAP_1 = median(UMAP_1),
            UMAP_2 = median(UMAP_2)) %>% 
  mutate(identShort = gsub('Clinical Cluster ', '', ident))

p = ggplot()
for (j in sort(unique(plotDf$exp))) {
  p = p + 
    geom_point(data = plotDf %>% filter(exp == j),
               aes(UMAP_1, UMAP_2, col = exp),
               size = 0.5)
}

p +
  geom_text(data = labDf, aes(UMAP_1, UMAP_2, label = identShort)) +
  facet_grid(treat ~ feat+time) +
  coord_fixed() +
  scale_color_gradient(low = 'gray', high = 'blue') +
  theme_bw()

ggplot(plotDf, aes(UMAP_1, UMAP_2, col = exp)) +
  geom_point(size = 0.5) +
  facet_grid(treat ~ feat+time) +
  coord_fixed() +
  scale_color_gradient(low = 'white', high = 'blue')

mkG <- cosg(datNew)

fts <- c('LAG3', 'PDCD1', 'TIGIT', 'IFNG', 'TOX', 'STMN1', 'MKI67', 'CCNA2')

lgd <- get_legend(
  DotPlot(datNew, features = fts)$data %>%
    ggplot(aes(x = features.plot, y = id, 
               size =  pct.exp,
               col = avg.exp.scaled)) +
    geom_point() +
    scale_color_gradient2(high = 'red', low = 'blue', mid = 'white', midpoint = 0) +
    labs(col = 'Average \nExpression', size = 'Percent \nExpressed') +
    theme(legend.position = 'bottom',
          text = element_text(size = 15, face = 2),
          legend.key=element_blank())
)

mkP <- DotPlot(datNew, features = fts)$data %>%
  ggplot(aes(x = features.plot, y = id, 
             size =  pct.exp,
             fill = ifelse(avg.exp.scaled > 0, 'UP', 'Dn'))) +
    geom_point(aes(alpha = abs(avg.exp.scaled)/max(abs(avg.exp.scaled))), shape = 21) +
  scale_y_discrete(limits=rev) +
  scale_fill_manual(values = c('blue', 'red')) +
  labs(x = 'Feature', y = 'Identity') +
  theme_bw() +
  theme(legend.position = 'none',
        text = element_text(size = 15, face = 2))

cowplot::plot_grid(mkP, lgd, 
                   ncol = 1, rel_heights = c(6, 1))
ggsave('3.Output/CellFigure/tonyClustMarker.png', width = 6, height = 5, units = 'in')

# figure generation -----
datSub <- subset(dat, idents = 0)

plotDf <- data.frame(Embeddings(datSub, reduction = 'umap'),
                     treat = datSub$treatment,
                     sample = datSub$sample_type,
                     pid = datSub$patient_id) %>% 
  filter(treat != 'off_study')

# hex fiugre 
ggplot(plotDf, aes(x = UMAP_1, y = UMAP_2))  +
  geom_hex(bins = 30) +
  scale_fill_continuous(type = "viridis") + 
  coord_fixed() +
  ggforce::theme_no_axes()+
  theme(strip.background = element_blank(),
        strip.text = element_text(face = 2, size = 15))

ggsave('3.Output/CellFigure/tonyHexUMAPOverall.png', width = 5, height = 6, 
       units = 'in')

ggplot(plotDf, aes(x = UMAP_1, y = UMAP_2))  +
  geom_hex(bins = 30) +
  facet_grid(~sample) +
  scale_fill_continuous(type = "viridis") + 
  coord_fixed() +
  ggforce::theme_no_axes()+
  theme(strip.background = element_blank(),
        strip.text = element_text(face = 2, size = 15))

ggsave('3.Output/CellFigure/tonyHexUMAPBySample.png', width = 10, height = 6, 
       units = 'in')

ggplot(plotDf, aes(x = UMAP_1, y = UMAP_2))  +
  geom_hex(bins = 30) +
  facet_grid(treat~sample) +
  scale_fill_continuous(type = "viridis") + 
  coord_fixed() +
  ggforce::theme_no_axes()+
  theme(strip.background = element_blank(),
        strip.text = element_text(face = 2, size = 15))

ggsave('3.Output/CellFigure/tonyHexUMAPTreatBySample.png', width = 10, height = 18, 
       units = 'in')

# maybe useful ----
# normal point figure
ggplot(plotDf, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(shape = 21, stroke = 1, fill = '#ffd322',  col = 'black', size = 2.2) +
  geom_point(shape = 16, col = '#ffd322', size = 2) +
  coord_fixed() +
  ggforce::theme_no_axes()

ggsave('3.Output/CellFigure/tonyUMAPOverall.png', width = 5, height = 6, 
       units = 'in')


ggplot(plotDf, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(shape = 21, stroke = 1, fill = '#ffd322',  col = 'black', size = 2.2) +
  geom_point(shape = 16, col = '#ffd322', size = 2) +
  facet_grid(~sample) +
  coord_fixed() +
  ggforce::theme_no_axes() +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = 2, size = 15))

ggsave('3.Output/CellFigure/tonyUMAPBySample.png', width = 10, height = 6, 
       units = 'in')

ggplot(plotDf, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(shape = 21, stroke = 1, fill = '#ffd322',  col = 'black', size = 2.2) +
  geom_point(shape = 16, col = '#ffd322', size = 2) +
  facet_grid(treat~sample) +
  coord_fixed() +
  ggforce::theme_no_axes() +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = 2, size = 15))

ggsave('3.Output/CellFigure/tonyUMAPTreatBySample.png', width = 10, height = 18, 
       units = 'in')


# shadow fiugre 
ggplot(plotDf, aes(x = UMAP_1, y = UMAP_2))  +
  stat_density_2d(aes(fill = ..level..), geom = "polygon") +
  scale_fill_distiller(palette=4, direction= 1) +
  coord_fixed() +
  ggforce::theme_no_axes()+
  theme(strip.background = element_blank(),
        strip.text = element_text(face = 2, size = 15))

ggsave('3.Output/CellFigure/tonyDensityUMAPOverall.png', width = 5, height = 6, 
       units = 'in')



ggplot(plotDf, aes(x = UMAP_1, y = UMAP_2)) +
  stat_density_2d(aes(fill = ..level..), geom = "polygon") +
  scale_fill_distiller(palette=4, direction= 1) +
  facet_grid(~sample) +
  coord_fixed() +
  ggforce::theme_no_axes() +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = 2, size = 15))

ggsave('3.Output/CellFigure/tonyDensityUMAPBySample.png', width = 10, height = 6, 
       units = 'in')

ggplot(plotDf, aes(x = UMAP_1, y = UMAP_2)) +
  stat_density_2d(aes(fill = ..level..), geom = "polygon") +
  scale_fill_distiller(palette=4, direction= 1) +
  facet_grid(treat~sample) +
  coord_fixed() +
  ggforce::theme_no_axes() +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = 2, size = 15))

ggsave('3.Output/CellFigure/tonyDensityUMAPTreatBySample.png', width = 10, height = 18, 
       units = 'in')

datSub <- subset(dat, sample_type == 'TIL')

mkG <- cosg(dat)

mkG$names[, '0']


DotPlot(subset(dat, idents = 0), split.by = 'sample_type', features = mkG$names[, '0']) +
  theme(axis.text.x = element_text(angle = 90))

DimPlot(dat, label = T) +
  coord_fixed()

dat$hlgt <-  ifelse(dat$leiden_0.5 == 0, 'aim', 'other')


DimPlot(subset(dat, idents = 0), label = T) +
  coord_fixed()

plotDf <- data.frame(Embeddings(dat, reduction = 'umap'),
                     clst = dat$leiden_0.5) %>% 
  mutate(hlgt = ifelse(clst == 0, 'aim', 'other')) %>%
  group_by(hlgt) %>% 
  sample_frac(0.3)

ggplot(plotDf, aes(x = UMAP_1, y = UMAP_2, col = hlgt)) +
  geom_point()

p <- ggplot()

for (i in rev(unique(plotDf$hlgt))) {
  p = p +
    geom_point(data = plotDf %>% filter(hlgt == i), 
               aes(UMAP_1, UMAP_2, fill = hlgt), shape = 21, stroke = 1) +
    geom_point(data = plotDf %>% filter(hlgt == i), 
               aes(UMAP_1, UMAP_2, col = hlgt), shape = 16)
}

# wwith label 
p +
  scale_color_manual(values = c('#ffd322', 'gray50')) +
  scale_fill_manual(values = c('#ffd322', 'gray50')) +
  coord_fixed() +
  theme_void() +
  NoLegend()

# cell population stats ----
head(plotDf)

prop.table(table(plotDf$pid, plotDf$sample, plotDf$treat), margin = 2) %>% 
  data.frame() %>% 
  ggplot(aes(x = Var2, y = Freq)) +
  geom_boxplot()+
  stat_compare_means(method = 'wilcox.test')


prop.table(table(plotDf$pid, plotDf$sample, plotDf$treat), margin = 1) %>% 
  data.frame() %>% # filter(Var1 == 'UoP002')
  ggplot(aes(x = Var2, y = Freq)) +
  geom_boxplot()+
  facet_grid(~Var3) +
  stat_compare_means(method = 'wilcox.test')


prop.table(table(c(1, 1,  2, 2, 3)))
