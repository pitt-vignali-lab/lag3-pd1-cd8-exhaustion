# ----------------------------------------------------------------------------
# Fig 2a-d | Labelled and unlabelled versions of the UMAP and distribution bars
#
# Paper    : LAG-3 and PD-1 synergize on CD8+ T cells to drive T cell
#            exhaustion and hinder autocrine IFN-gamma-dependent anti-tumor
#            immunity
# Original : 2.Script/rev1_clusterLab&mkGFormalFigure.R
# Input    : 1.Data/datForRevision.RDS
# Output   : 3.Output/CellFigure/overallUMAP*.png
#            3.Output/CellFigure/statBar*_poly.png
#            3.Output/CellFigure/statBar*_pMEL.png
# Upstream : 04_cluster_annotation.R
# Note     : Uses SCpubr for the final panel styling
# ----------------------------------------------------------------------------

library(Seurat)
library(tidyverse)
library(SCpubr)
library(ggrepel)
library(grid)

dat <- readRDS('1.Data/datForRevision.RDS')


do_DimPlot(dat, group.by = 'anno') +
  coord_fixed()

plotDf <- data.frame(Embeddings(dat, reduction = 'umap'),
                     clst = dat$anno,
                     gntp = dat$genotype,
                     smpl = dat$sample) %>% 
  mutate(gntp = factor(gntp, levels = c("E8ICre.GFP (WT)", 
                                            "Pdcd1L/L Lag3L/L-yfp (WT)",
                                            "Lag3L/L-yfp E8ICre.GFP",
                                            "Pdcd1L/L E8ICre.GFP",
                                            "Pdcd1L/L Lag3L/L-yfp E8ICre.GFP")))
labDf <- plotDf %>% 
  group_by(clst) %>% 
  summarise(UMAP_1 = median(UMAP_1),
            UMAP_2 = median(UMAP_2))

ggplot(plotDf, aes(UMAP_1, UMAP_2, col = clst)) +
  geom_point()

# Overall UMAP ----
p <- ggplot()

for (i in rev(unique(plotDf$clst))) {
 p = p +
   geom_point(data = plotDf %>% filter(clst == i), 
              aes(UMAP_1, UMAP_2, fill = clst), shape = 21, stroke = 1) +
   geom_point(data = plotDf %>% filter(clst == i), 
              aes(UMAP_1, UMAP_2, col = clst), shape = 16)
}

# wwith label 
p +
  geom_label_repel(data = labDf,
                   aes(UMAP_1, UMAP_2, label = clst, col = clst),
                   size = 5, fontface = 2) +
  coord_fixed() +
  theme_void() +
  NoLegend()

ggsave('3.Output/CellFigure/overallUMAPLabeled.png', width = 6, height = 6, units = 'in', bg = 'white')

p +
  geom_label_repel(data = labDf,
                   aes(UMAP_1, UMAP_2, label = gsub(':.*', '', clst), col = clst),
                   size = 6, fontface = 2) +
  coord_fixed() +
  theme_void() +
  NoLegend()

ggsave('3.Output/CellFigure/UMAPOverallLabeledSimple.png', width = 6, height = 6, units = 'in', bg = 'white')

# no label
p +
  coord_fixed() +
  theme_void() +
  NoLegend()

ggsave('3.Output/CellFigure/overallUMAPNoLabel.png', width = 5, height = 5, units = 'in', bg = 'white')

# split UMAP ----
# with label
subP <- p +
  facet_grid(smpl ~ gntp) +
  coord_fixed()+
  theme_bw() +
  theme(legend.position = "none",
        axis.ticks = element_blank(), 
        axis.text = element_blank(),
        axis.title = element_blank(),
        strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

g <- ggplotGrob(subP)
rm_grobs <- g$layout$name %in% c("panel-2-1")
g$grobs[rm_grobs] <- NULL
g$layout <- g$layout[!rm_grobs, ]
# g$layout[g$layout$name == "axis-b-3-3", c("t", "b")] = c(14.5, 14.5)
grid.newpage()

png('3.Output/CellFigure/splitUMAPLabeled.png', 
    width = 10, height = 5, units = "in", res = 300)
grid.draw(g)
dev.off()
# no label 
subP <- p +
  facet_grid(smpl ~ gntp) +
  coord_fixed()+
  theme_bw() +
  theme(legend.position = "none",
        axis.ticks = element_blank(), 
        axis.text = element_blank(),
        axis.title = element_blank(),
        strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_blank())
g <- ggplotGrob(subP)
rm_grobs <- g$layout$name %in% c("panel-2-1")
g$grobs[rm_grobs] <- NULL
g$layout <- g$layout[!rm_grobs, ]
# g$layout[g$layout$name == "axis-b-3-3", c("t", "b")] = c(14.5, 14.5)
grid.newpage()

png('3.Output/CellFigure/splitUMAPNoLabel.png', 
    width = 10, height = 5, units = "in", res = 300)
grid.draw(g)
dev.off()
# statistic bars ----
# polyclonal 
t1 <- data.frame(prop.table(table(plotDf[plotDf$smpl == "Normal CD8",]$clst, 
                                  plotDf[plotDf$smpl == "Normal CD8",]$gntp)))

names(t1) <- c("Cluster", "Genotype", "Freq")

ggplot(t1, aes(x = Freq, y = Genotype, fill = Cluster)) +
  geom_bar(col = "black", 
           stat = "identity",
           position=position_fill(reverse=T)) +
  scale_x_continuous(expand = c(0.02, 0),
                     labels = scales::percent) +
  scale_y_discrete(expand = expansion(0, 1/4), 
                   limits=rev,
                   labels = c(bquote(bold(atop(italic(Pdcd1)^{L/L}~italic(Lag3)^{L/L-YFP}, E8I^Cre.GFP))), 
                              bquote(bold(italic(Pdcd1)^{L/L-YFP}~E8I^Cre.GFY)), 
                              bquote(bold(italic(Lag3)^{L/L-YFP}~E8I^Cre.GFY)), 
                              bquote(bold(italic(Pdcd1)^{L/L}~italic(Lag3)^{L/L-YFP})), 
                              bquote(bold(E8I^Cre.GFP)))) +
  xlab("Percentage of sample") +
  ylab("") +
  ggtitle(bquote(bold(Polyclonal~CD8^'+'~TIL))) +
  coord_fixed(ratio = 0.05) +
  theme(text = element_text(family = "Arial",
                            size = 20, face = 2),
        axis.text.y = element_text(color = c("purple", "red", "blue", "black", "gray50")),
        legend.position = "none",
        # legend.direction = 'horizontal',
        # legend.text = element_text(size = 12), 
        plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.ticks.y = element_blank(),
        strip.background = element_blank())+
  guides(fill = guide_legend(nrow = 2))

graph2ppt(file = '3.Output/CellFigure/statBarLabeled_poly.pptx')
ggsave('3.Output/CellFigure/statBarLabeled_poly.png', width = 7, height = 3, units = 'in', bg = 'white')

ggplot(t1, aes(x = Freq, y = Genotype, fill = Cluster)) +
  geom_bar(col = "black", 
           stat = "identity",
           position=position_fill(reverse=T)) +
  scale_x_continuous(expand = c(0.02, 0),
                     labels = scales::percent) +
  scale_y_discrete(expand = expansion(0, 1/4), 
                   limits=rev,
                   labels = c(bquote(bold(atop(italic(Pdcd1)^{L/L}~italic(Lag3)^{L/L-YFP}, E8I^Cre.GFP))), 
                              bquote(bold(italic(Pdcd1)^{L/L-YFP}~E8I^Cre.GFY)), 
                              bquote(bold(italic(Lag3)^{L/L-YFP}~E8I^Cre.GFY)), 
                              bquote(bold(italic(Pdcd1)^{L/L}~italic(Lag3)^{L/L-YFP})), 
                              bquote(bold(E8I^Cre.GFP)))) +
  xlab("") +
  ylab("") +
  # ggtitle(bquote(bold(Polyclonal~CD8^'+'~TIL))) +
  coord_fixed(ratio = 0.05) +
  theme(text = element_blank(),
        # axis.text.y = element_text(color = c("purple", "red", "blue", "black", "gray50")),
        legend.position = "none",
        # legend.direction = 'horizontal',
        # legend.text = element_text(size = 12), 
        # plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.ticks.y = element_blank(),
        strip.background = element_blank())+
  guides(fill = guide_legend(nrow = 2))
graph2ppt(file = '3.Output/CellFigure/statBarNoLabel_poly.pptx')
ggsave('3.Output/CellFigure/statBarNoLabel_poly.png', width = 7, height = 3, units = 'in', bg = 'white')

# pMEL
t2 <- data.frame(prop.table(table(plotDf[plotDf$smpl == "PML CD8",]$clst, 
                                  plotDf[plotDf$smpl == "PML CD8",]$gntp)))

names(t2) <- c("Cluster", "Genotype", "Freq")

t2 <- t2 %>% filter(Genotype != "E8ICre.GFP (WT)")

t2$Genotype <- factor(t2$Genotype, levels = c("NA", "Pdcd1L/L Lag3L/L-yfp (WT)", 
                                              "Lag3L/L-yfp E8ICre.GFP",
                                              "Pdcd1L/L E8ICre.GFP",
                                              "Pdcd1L/L Lag3L/L-yfp E8ICre.GFP"))

ggplot(t2, aes(x = Freq, y = Genotype, fill = Cluster)) +
  geom_bar(col = "black", 
           stat = "identity",
           position=position_fill(reverse=T)) +
  scale_x_continuous(expand = c(0.02, 0),
                     labels = scales::percent) +
  scale_y_discrete(expand = expansion(0, 1/5),
                   limits=rev,
                   labels = c(bquote(bold(atop(italic(Pdcd1)^{L/L}~italic(Lag3)^{L/L-YFP}, E8I^Cre.GFP))),
                              bquote(bold(italic(Pdcd1)^{L/L-YFP}~E8I^Cre.GFY)),
                              bquote(bold(italic(Lag3)^{L/L-YFP}~E8I^Cre.GFY)),
                              bquote(bold(italic(Pdcd1)^{L/L}~italic(Lag3)^{L/L-YFP}))
                              # bquote(bold(E8I^Cre.GFP))
                   )) +
  xlab("Percentage of sample") +
  ylab("") +
  coord_fixed(ratio = 0.05) +
  ggtitle(bquote(bold(pMEL~CD8^'+'~TIL))) +
  theme(text = element_text(family = "Arial",
                            size = 20, face = 2),
        axis.text.y = element_text(color = c("purple", "red", "blue", "black")),
        legend.position = "none",
        # legend.direction = 'horizontal',
        # legend.text = element_text(size = 12), 
        plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.ticks.y = element_blank(),
        strip.background = element_blank())+
  guides(fill = guide_legend(nrow = 2))

graph2ppt(file = '3.Output/CellFigure/statBarLabeled_pMEL.pptx')
ggsave('3.Output/CellFigure/statBarLabeled_pMEL.png', width = 7, height = 3, units = 'in', bg = 'white')

ggplot(t2, aes(x = Freq, y = Genotype, fill = Cluster)) +
  geom_bar(col = "black", 
           stat = "identity",
           position=position_fill(reverse=T)) +
  scale_x_continuous(expand = c(0.02, 0),
                     labels = scales::percent) +
  scale_y_discrete(expand = expansion(0, 1/5),
                   limits=rev,
                   labels = c(bquote(bold(atop(italic(Pdcd1)^{L/L}~italic(Lag3)^{L/L-YFP}, E8I^Cre.GFP))),
                              bquote(bold(italic(Pdcd1)^{L/L-YFP}~E8I^Cre.GFY)),
                              bquote(bold(italic(Lag3)^{L/L-YFP}~E8I^Cre.GFY)),
                              bquote(bold(italic(Pdcd1)^{L/L}~italic(Lag3)^{L/L-YFP}))
                              # bquote(bold(E8I^Cre.GFP))
                   )) +
  xlab("") +
  ylab("") +
  coord_fixed(ratio = 0.05) +
  # ggtitle(bquote(bold(pMEL~CD8^'+'~TIL))) +
  theme(text = element_blank(),
        # axis.text.y = element_text(color = c("purple", "red", "blue", "black")),
        legend.position = "none",
        # legend.direction = 'horizontal',
        # legend.text = element_text(size = 12), 
        plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.ticks.y = element_blank(),
        strip.background = element_blank())+
  guides(fill = guide_legend(nrow = 2))

graph2ppt(file = '3.Output/CellFigure/statBarNoLabel_pMEL.pptx')
ggsave('3.Output/CellFigure/statBarNoLabel_pMEL.png', width = 7, height = 3, units = 'in', bg = 'white')
