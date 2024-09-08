# ----------------------------------------------------------------------------
# Fig 2a-f | Render the published Figure 2 UMAP, bar and marker panels
#
# Paper    : LAG-3 and PD-1 synergize on CD8+ T cells to drive T cell
#            exhaustion and hinder autocrine IFN-gamma-dependent anti-tumor
#            immunity
# Original : 2.Script/PictureFroManual/fig2.R
# Input    : 1.Data/datInt.RDS
# Output   : 3.Output/FormalFigures/2A_Overall_umap.png
#            3.Output/FormalFigures/2B_split_umap.png
#            3.Output/FormalFigures/2C_polycolonal_bar.png
#            3.Output/FormalFigures/2D_pMEL_bar.png
# Upstream : 02_cluster_labeling.R
# ----------------------------------------------------------------------------

library(Seurat)
library(tidyverse)
library(extrafont)
library(scales)
library(ggthemes)

if(!dir.exists("3.Output/FormalFigures")){
  dir.create("3.Output/FormalFigures")
}

loadfonts(device = "win", quiet = TRUE)

dat <- readRDS("1.Data/datInt.RDS")

plt <- DimPlot(dat)

# x-range
round(diff(layer_scales(plt)$x$range$range)/2, 1)
# y-range
round(diff(layer_scales(plt)$y$range$range)/2, 1)

plotDf <- data.frame(Embeddings(dat, reduction = "umap"),
                     Idents(dat),
                     dat$genotype,
                     dat$sample)

names(plotDf)[3:5] <- c("ident", "genotype", "sample")

plotDf$ident <- paste0("C", plotDf$ident)

plotDf$genotype <- factor(plotDf$genotype, levels = c("E8ICre.GFP (WT)", 
                                                      "Pdcd1L/L Lag3L/L-yfp (WT)",
                                                      "Lag3L/L-yfp E8ICre.GFP",
                                                      "Pdcd1L/L E8ICre.GFP",
                                                      "Pdcd1L/L Lag3L/L-yfp E8ICre.GFP"))

labDf <- plotDf %>% group_by(ident) %>% 
  summarise(x = mean(UMAP_1), y = mean(UMAP_2))

# 2A ----
ggplot(plotDf, aes(UMAP_1, UMAP_2, col = ident)) +
  geom_point(size = 1.8, alpha = 0.8) +
  coord_fixed() +
  theme_void() +
  guides(color = guide_legend(title = "Clusters", 
                              nrow=2,
                              byrow=TRUE,
                              override.aes = list(size=3))) +
  theme(legend.position = c(0.6, 0.1),
        legend.direction = "vertical",
        legend.title.align = 0.5,
        legend.key=element_rect(fill = "gray90", color = NA),
        # legend.box.background = element_rect(colour = "black"),
        legend.text = element_text(family = "Arial",
                                   size = 11.5, 
                                   face = 2),
        legend.title = element_text(family = "Arial",
                                   size = 12.5, 
                                   face = 2))
ggsave("3.Output/FormalFigures/2A_Overall_umap.png", 
       width = 6.6, height = 8.6, units = "in", dpi = 300, bg = "white")

# 2B ----
library(grid)

p <- ggplot(plotDf, aes(UMAP_1, UMAP_2, col = ident)) +
  geom_point(size = 2.5, alpha = 0.8) +
  coord_fixed() +
  facet_grid(sample ~ genotype) +
  theme_bw() +
  guides(color = guide_legend(title = "Clusters", 
                              nrow=2,
                              byrow=TRUE,
                              override.aes = list(size=3))) +
  theme(legend.position = "none",
        axis.ticks = element_blank(), 
        axis.text = element_blank(),
        axis.title = element_blank(),
        strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

p

g <- ggplotGrob(p)
rm_grobs <- g$layout$name %in% c("panel-2-1")
g$grobs[rm_grobs] <- NULL
g$layout <- g$layout[!rm_grobs, ]
# g$layout[g$layout$name == "axis-b-3-3", c("t", "b")] = c(14.5, 14.5)
grid.newpage()

png("3.Output/FormalFigures/2B_split_umap.png", 
    width = 15, height = 8, units = "in", res = 300)
grid.draw(g)
dev.off()


# 2C ----
t1 <- data.frame(prop.table(table(plotDf[plotDf$sample == "Normal CD8",]$ident, 
                                  plotDf[plotDf$sample == "Normal CD8",]$genotype)))

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

ggsave("3.Output/FormalFigures/2C_polycolonal_bar.png", 
       width = 15, height = 5, units = "in", dpi = 300, bg = "white")

# 2D ----
t2 <- data.frame(prop.table(table(plotDf[plotDf$sample == "PML CD8",]$ident, 
                                  plotDf[plotDf$sample == "PML CD8",]$genotype)))

names(t2) <- c("Cluster", "Genotype", "Freq")

t2 <- t2 %>% filter(Genotype != "E8ICre.GFP (WT)")

t2$Genotype <- factor(t2$Genotype, levels = c("Pdcd1L/L Lag3L/L-yfp (WT)", 
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

ggsave("3.Output/FormalFigures/2D_pMEL_bar.png", 
       width = 15, height = 4, units = "in", dpi = 300, bg = "white")

# legend 
library(grid)
library(gridExtra)
library(cowplot)

gplot <- ggplot(t1, aes(x = Freq, y = Genotype, fill = Cluster)) +
  geom_bar(col = "black", 
           stat = "identity") +
  theme(legend.direction = 'horizontal',
        legend.text = element_text(size = 11.5, face = 2),
        legend.title = element_text(size = 12.5, face = 2))
legend <- get_legend(gplot)
grid.newpage()
png("3.Output/FormalFigures/2C_2D_legend.png", 
    height = 1, width = 4, units = "in", res = 300)
grid.draw(legend) 
dev.off()

# 2E ----
DotPlot(dat, 
        assay = "RNA",
        features = c(strsplit("Klrg1, Cx3cr1, S1pr5, Tnf, Ifng, Gzmb, Prf1", ", ")[[1]], # high
                     strsplit("Sell, Tcf7, Lef1, Ccr7, Il7r", ", ")[[1]] # low
        ), 
        col.min = 0, cols = c("gray90", "blue"))+
  ggtitle(bquote(bold("Effector and/or effector memory-like CD8"^'+'~"cells marker genes (cluster 0)"))) +
  scale_y_discrete(labels = paste0("C", 0:9))+
  geom_rect(data = data.frame(start = c(0, 7.5), end = c(7.5, 12.5), 
                              Status = c("Shold be high", "Should be low")),
            aes(ymin=0, ymax=10.5, xmin=start,
                xmax=end, fill = Status), alpha =0.2, inherit.aes = FALSE ) +
  theme_excel_new() +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(family = "Arial",
                            size = 20, face = 2),
        panel.grid.major = element_blank(),
        legend.position = "none")

ggsave("3.Output/FormalFigures/2E_mkgene_cluster0.png", 
       width = 9, height = 5, units = "in", dpi = 300, bg = "white")

# 2F ----
DotPlot(dat, 
        assay = "RNA",
        features = c(strsplit("Bst2, Irf1, Irf2, Irf7, Stat1, Stat2, Mx1", ", ")[[1]]#, # high
                          # strsplit("Tcf7, Havcr2", ", ")[[1]] # low
                     ), col.min = 0, cols = c("gray90", "blue")) +
  ggtitle(bquote(bold("IFN-stimulated CD8"^'+'~"cells marker genes (cluster 7)"))) +
  scale_y_discrete(labels = paste0("C", 0:9))+
  geom_rect(data = data.frame(start = c(0), end = c(7.5), 
                              Status = c("Shold be high")),
            aes(ymin=0, ymax=10.5, xmin=start,
                xmax=end, fill = Status), alpha =0.2, inherit.aes = FALSE ) +
  theme_excel_new() +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(family = "Arial",
                            size = 20, face = 2),
        panel.grid.major = element_blank(),
        legend.position = "none")

ggsave("3.Output/FormalFigures/2F_mkgene_cluster7.png", 
       width = 5.25, height = 5, units = "in", dpi = 300, bg = "white")

# legend 

gplot <- DotPlot(dat, 
                 assay = "RNA",
                 features = c(strsplit("Klrg1, Cx3cr1, S1pr5, Tnf, Ifng, Gzmb, Prf1", ", ")[[1]], # high
                              strsplit("Sell, Tcf7, Lef1, Ccr7, Il7r", ", ")[[1]] # low
                 ), 
                 col.min = 0, cols = c("gray90", "blue"))+
  ggtitle(bquote(bold("Effector and/or effector memory-like CD8"^'+'~"cells marker genes (cluster 0)"))) +
  scale_y_discrete(labels = paste0("C", 0:9))+
  geom_rect(data = data.frame(start = c(0, 7.5), end = c(7.5, 12.5), 
                              Status = c("Shold be high", "Should be low")),
            aes(ymin=0, ymax=10.5, xmin=start,
                xmax=end, fill = Status), alpha =0.2, inherit.aes = FALSE ) + 
  # theme_excel_new() +
  # guides(col= guide_colorbar(barwidth = unit(2, "in"))) +
  theme(legend.text = element_text(size = 11.5, face = 2))

legend <- get_legend(gplot)
grid.newpage()
grid.draw(legend) 

png("3.Output/FormalFigures/2E_2F_legend.png", 
    height = 4, width = 2, units = "in", res = 300)
grid.draw(legend) 
dev.off()
