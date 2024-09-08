# ----------------------------------------------------------------------------
# Fig 2 | Volcano plots of differentially expressed genes per cluster
#
# Paper    : LAG-3 and PD-1 synergize on CD8+ T cells to drive T cell
#            exhaustion and hinder autocrine IFN-gamma-dependent anti-tumor
#            immunity
# Original : 2.Script/rev1_deGVolcanoPlot.R
# Input    : 1.Data/datForRevision.RDS
# Output   : 3.Output/Revision/volcanoClst*.png
# Upstream : 04_cluster_annotation.R
# Note     : presto wilcoxauc
# ----------------------------------------------------------------------------

library(tidyverse)
library(presto)
library(ggrepel)
library(Seurat)
library(extrafont)

dat <- readRDS('1.Data/datForRevision.RDS')

Idents(dat) <- 'anno'

DimPlot(dat, label = T) +
  coord_fixed()

deG <- wilcoxauc(dat)

# volcano plot with highlighted genes (cluster 0 7 )----
## cluster 0
deGSub <- deG %>% 
  filter(group == "C0: Effector/effector memory-like", 
         padj < 0.05) %>% 
  mutate(sig = -log(padj)) %>% 
  mutate(sig = ifelse(is.infinite(sig),
                      max(sig[!is.infinite(sig)]) +  runif(length(sig[is.infinite(sig)])) * 550,
                      sig))

hg = c('Ifngr1', 'Gzmb', 'Nkg7', 'Cxcr6', 'Ccl5')

deGSub %>% 
  ggplot(aes(x = logFC, 
             y = sig, 
             col = ifelse(logFC > 0, 'UP', 'DN'))) +
  geom_point(alpha = 0.3, size = 3) +
  geom_label_repel(data = deGSub %>% filter(feature %in% c(hg,
                                                           deGSub %>% filter(!feature %in% hg) %>% 
                                                             group_by(logFC > 0) %>% 
                                                             top_n(5, abs(logFC)) %>% 
                                                             pull(feature))),
                   aes(label = feature, 
                       alpha = ifelse(feature %in% hg, 'hl', 'nm')), 
                   size = 5, fontface = 2,
                   show.legend = F) +
  scale_y_log10() +
  xlim(-max(abs(deGSub$logFC)), max(abs(deGSub$logFC))) +
  scale_color_manual(breaks = c('DN', 'UP'),
                     values = c('blue', 'red'),
                     labels = c('Down-regulated genes',
                                'Up-regulated genes')) +
  scale_alpha_manual(breaks = c('hl', 'nm'),
                     values = c(1, 0.5)) +
  labs(x = bquote(bold(log^{fold~change})),
       y = expression(bold(atop(Significance, (-log^{adjusted~italic(p)~value})))),
       col = 'Regulation \nstatus', max.overlaps = Inf) +
  guides(color = guide_legend(override.aes = list(size=3)))+
  # coord_fixed(ratio = ratio.values) +
  theme_bw() +
  theme(legend.position = "bottom",
        text = element_text(family = "Arial", size = 18, face = 2), 
        aspect.ratio=1)

ggsave('3.Output/CellFigure/volcanoClst0Eff_labeled.png', 
       width = 8, height = 5, units = "in", dpi = 300, bg = "white")

deGSub %>% 
  ggplot(aes(x = logFC, 
             y = sig, 
             col = ifelse(logFC > 0, 'UP', 'DN'))) +
  geom_point(alpha = 0.3, size = 3) +
  geom_label_repel(data = deGSub %>% filter(feature %in% c(hg,
                                                           deGSub %>% filter(!feature %in% hg) %>% 
                                                             group_by(logFC > 0) %>% 
                                                             top_n(5, abs(logFC)) %>% 
                                                             pull(feature))),
                   aes(label = feature, 
                       alpha = ifelse(feature %in% hg, 'hl', 'nm')), 
                   size = 6, fontface = 2, family = "Arial",
                   show.legend = F) +
  scale_y_log10() +
  xlim(-max(abs(deGSub$logFC)), max(abs(deGSub$logFC))) +
  scale_color_manual(breaks = c('DN', 'UP'),
                     values = c('blue', 'red'),
                     labels = c('Down-regulated genes',
                                'Up-regulated genes')) +
  scale_alpha_manual(breaks = c('hl', 'nm'),
                     values = c(1, 0.5)) +
  labs(x = '',
       y = '',
       col = '') +
  guides(color = guide_legend(override.aes = list(size=3)))+
  theme_bw() +
  theme(legend.position = "none",
        axis.text = element_blank(), 
        aspect.ratio=1)

ggsave('3.Output/CellFigure/volcanoClst0Eff_Nolabel.png', 
       width = 5, height = 5, units = "in", dpi = 300, bg = "white")

# ggsave('3.Output/Revision/volcanoClst0Eff.png', 
#        width = 8, height = 5, units = "in", dpi = 300, bg = "white")

## Cluster 7
deGSub <- deG %>% 
  filter(group == "C7: IFN-induced", 
         padj < 0.05)%>% 
  mutate(sig = -log(padj)) %>% 
  mutate(sig = ifelse(is.infinite(sig),
                      max(sig[!is.infinite(sig)]) +  runif(length(sig[is.infinite(sig)])) * 550,
                      sig))

hg = c('Ifit3', 'Gzmk', 'Cxcr3', 'Ifngr1', 'Ccl5')

deGSub %>% 
  ggplot(aes(x = logFC, 
             y = sig, 
             col = ifelse(logFC > 0, 'UP', 'DN'))) +
  geom_point(alpha = 0.3, size = 3) +
  geom_label_repel(data = deGSub %>% filter(feature %in% c(hg,
                                                           deGSub %>% filter(!feature %in% hg) %>% 
                                                             group_by(logFC > 0) %>% 
                                                             top_n(5, abs(logFC)) %>% 
                                                             pull(feature))),
                   aes(label = feature, 
                       alpha = ifelse(feature %in% hg, 'hl', 'nm')), 
                   size = 5, fontface = 2,
                   show.legend = F) +
  scale_y_log10() +
  xlim(-max(abs(deGSub$logFC)), max(abs(deGSub$logFC))) +
  scale_color_manual(breaks = c('DN', 'UP'),
                     values = c('blue', 'red'),
                     labels = c('Down-regulated genes',
                                'Up-regulated genes')) +
  scale_alpha_manual(breaks = c('hl', 'nm'),
                     values = c(1, 0.5)) +
  labs(x = bquote(bold(log^{fold~change})),
       y = expression(bold(atop(Significance, (-log^{adjusted~italic(p)~value})))),
       col = 'Regulation \nstatus') +
  guides(color = guide_legend(override.aes = list(size=3)))+
  theme_bw() +
  theme(legend.position = "bottom",
        text = element_text(family = "Arial", size = 18, face = 2), 
        aspect.ratio=1)

ggsave('3.Output/CellFigure/volcanoClst7Ifn_labeled.png', 
       width = 8, height = 5, units = "in", dpi = 300, bg = "white")

deGSub %>% 
  ggplot(aes(x = logFC, 
             y = sig, 
             col = ifelse(logFC > 0, 'UP', 'DN'))) +
  geom_point(alpha = 0.3, size = 3) +
  geom_label_repel(data = deGSub %>% filter(feature %in% c(hg,
                                                           deGSub %>% filter(!feature %in% hg) %>% 
                                                             group_by(logFC > 0) %>% 
                                                             top_n(5, abs(logFC)) %>% 
                                                             pull(feature))),
                   aes(label = feature, 
                       alpha = ifelse(feature %in% hg, 'hl', 'nm')), 
                   size = 6, fontface = 2, family = "Arial",
                   show.legend = F) +
  scale_y_log10() +
  xlim(-max(abs(deGSub$logFC)), max(abs(deGSub$logFC))) +
  scale_color_manual(breaks = c('DN', 'UP'),
                     values = c('blue', 'red'),
                     labels = c('Down-regulated genes',
                                'Up-regulated genes')) +
  scale_alpha_manual(breaks = c('hl', 'nm'),
                     values = c(1, 0.5)) +
  labs(x = '',
       y = '',
       col = '') +
  guides(color = guide_legend(override.aes = list(size=3)))+
  theme_bw() +
  theme(legend.position = "none",
        axis.text = element_blank(), 
        aspect.ratio=1)

ggsave('3.Output/CellFigure/volcanoClst7Ifn_Nolabel.png', 
       width = 5, height = 5, units = "in", dpi = 300, bg = "white")

## Cluster 4
deGSub <- deG %>% 
  filter(group == "C4: Exhuasted", 
         padj < 0.05)%>% 
  mutate(sig = -log(padj)) %>% 
  mutate(sig = ifelse(is.infinite(sig),
                      max(sig[!is.infinite(sig)]) +  runif(length(sig[is.infinite(sig)])) * 550,
                      sig))

deGSub %>% 
  ggplot(aes(x = logFC, 
             y = sig, 
             col = ifelse(logFC > 0, 'UP', 'DN'))) +
  geom_point(alpha = 0.3, size = 3) +
  geom_label_repel(data = deGSub %>% 
                     # filter(!feature %in% hg) %>% 
                     group_by(logFC > 0) %>% 
                     top_n(5, abs(logFC)),
                   aes(label = feature), 
                   size = 6, fontface = 2, alpha = 0.8,family = "Arial", 
                   show.legend = F) +
  scale_y_log10() +
  scale_color_manual(breaks = c('DN', 'UP'),
                     values = c('blue', 'red'),
                     labels = c('Down-regulated genes',
                                'Up-regulated genes')) +
  labs(x = bquote(bold(log^{fold~change})),
       y = expression(bold(atop(Significance, (-log^{adjusted~italic(p)~value})))),
       col = 'Regulation \nstatus') +
  guides(color = guide_legend(override.aes = list(size=3)))+
  xlim(-max(abs(deGSub$logFC)), max(abs(deGSub$logFC))) +
  theme_bw() +
  theme(legend.position = "bottom",
        text = element_text(family = "Arial", size = 18, face = 2), 
        aspect.ratio=1)

ggsave('3.Output/CellFigure/volcanoClst4Exh_labeled.png', 
       width = 8, height = 5, units = "in", dpi = 300, bg = "white")




deGSub %>% 
  ggplot(aes(x = logFC, 
             y = sig, 
             col = ifelse(logFC > 0, 'UP', 'DN'))) +
  geom_point(alpha = 0.3, size = 3) +
  geom_label_repel(data = deGSub %>% 
                     # filter(!feature %in% hg) %>% 
                     group_by(logFC > 0) %>% 
                     top_n(5, abs(logFC)),
                   aes(label = feature), 
                   size = 6, fontface = 2, alpha = 0.8,family = "Arial", 
                   show.legend = F) +
  scale_y_log10() +
  xlim(-max(abs(deGSub$logFC)), max(abs(deGSub$logFC))) +
  scale_color_manual(breaks = c('DN', 'UP'),
                     values = c('blue', 'red'),
                     labels = c('Down-regulated genes',
                                'Up-regulated genes')) +
  scale_alpha_manual(breaks = c('hl', 'nm'),
                     values = c(1, 0.5)) +
  labs(x = '',
       y = '',
       col = '') +
  guides(color = guide_legend(override.aes = list(size=3)))+
  theme_bw() +
  theme(legend.position = "none",
        axis.text = element_blank(), 
        aspect.ratio=1)

ggsave('3.Output/CellFigure/volcanoClst4Exh_Nolabel.png', 
       width = 5, height = 5, units = "in", dpi = 300, bg = "white")


# ggsave('3.Output/Revision/volcanoClst7Ifn.png', 
#        width = 8, height = 5, units = "in", dpi = 300, bg = "white")

# volcano plot without highlighted genes (cluster 1 3 )----
## cluster 3
deGSub <- deG %>% 
  filter(group == "C3: Proliferation", 
         padj < 0.05)%>% 
  mutate(sig = -log(padj)) %>% 
  mutate(sig = ifelse(is.infinite(sig),
                      max(sig[!is.infinite(sig)]) +  runif(length(sig[is.infinite(sig)])) * 550,
                      sig))

deGSub %>% 
  ggplot(aes(x = logFC, 
             y = sig, 
             col = ifelse(logFC > 0, 'UP', 'DN'))) +
  geom_point(alpha = 0.3, size = 3) +
  geom_label_repel(data = deGSub %>% 
                     # filter(!feature %in% hg) %>% 
                     group_by(logFC > 0) %>% 
                     top_n(5, abs(logFC)),
                   aes(label = feature), 
                   size = 6, fontface = 2, alpha = 0.8,family = "Arial", 
                   show.legend = F) +
  scale_y_log10() +
  scale_color_manual(breaks = c('DN', 'UP'),
                     values = c('blue', 'red'),
                     labels = c('Down-regulated genes',
                                'Up-regulated genes')) +
  labs(x = bquote(bold(log^{fold~change})),
       y = expression(bold(atop(Significance, (-log^{adjusted~italic(p)~value})))),
       col = 'Regulation \nstatus') +
  guides(color = guide_legend(override.aes = list(size=3)))+
  xlim(-max(abs(deGSub$logFC)), max(abs(deGSub$logFC))) +
  theme_bw() +
  theme(legend.position = "bottom",
        text = element_text(family = "Arial", size = 18, face = 2), 
        aspect.ratio=1)

ggsave('3.Output/CellFigure/volcanoClst3Prolif_labeled.png', 
       width = 8, height = 5, units = "in", dpi = 300, bg = "white")




deGSub %>% 
  ggplot(aes(x = logFC, 
             y = sig, 
             col = ifelse(logFC > 0, 'UP', 'DN'))) +
  geom_point(alpha = 0.3, size = 3) +
  geom_label_repel(data = deGSub %>% 
                     # filter(!feature %in% hg) %>% 
                     group_by(logFC > 0) %>% 
                     top_n(5, abs(logFC)),
                   aes(label = feature), 
                   size = 6, fontface = 2, alpha = 0.8,family = "Arial", 
                   show.legend = F) +
  scale_y_log10() +
  xlim(-max(abs(deGSub$logFC)), max(abs(deGSub$logFC))) +
  scale_color_manual(breaks = c('DN', 'UP'),
                     values = c('blue', 'red'),
                     labels = c('Down-regulated genes',
                                'Up-regulated genes')) +
  scale_alpha_manual(breaks = c('hl', 'nm'),
                     values = c(1, 0.5)) +
  labs(x = '',
       y = '',
       col = '') +
  guides(color = guide_legend(override.aes = list(size=3)))+
  theme_bw() +
  theme(legend.position = "none",
        axis.text = element_blank(), 
        aspect.ratio=1)

ggsave('3.Output/CellFigure/volcanoClst3Prolif_Nolabel.png', 
       width = 5, height = 5, units = "in", dpi = 300, bg = "white")


## cluster 1
# deGSub <- deG %>% 
#   filter(group == "C1: Activated", 
#          padj < 0.05)%>% 
#   mutate(sig = -log(padj)) %>% 
#   mutate(sig = ifelse(is.infinite(sig),
#                       max(sig[!is.infinite(sig)]) +  runif(length(sig[is.infinite(sig)])) * 550,
#                       sig))
# 
# deGSub %>% 
#   ggplot(aes(x = logFC, 
#              y = sig, 
#              col = ifelse(logFC > 0, 'UP', 'DN'))) +
#   geom_point(alpha = 0.3, size = 3) +
#   geom_label_repel(data = deGSub %>% 
#                      # filter(!feature %in% hg) %>% 
#                      group_by(logFC > 0) %>% 
#                      top_n(5, abs(logFC)),
#                    aes(label = feature), 
#                    size = 6, fontface = 2, alpha = 0.8,family = "Arial", 
#                    show.legend = F) +
#   scale_y_log10() +
#   scale_color_manual(breaks = c('DN', 'UP'),
#                      values = c('blue', 'red'),
#                      labels = c('Down-regulated genes',
#                                 'Up-regulated genes')) +
#   labs(x = bquote(bold(log^{fold~change})),
#        y = expression(bold(atop(Significance, (-log^{adjusted~italic(p)~value})))),
#        col = 'Regulation \nstatus') +
#   guides(color = guide_legend(override.aes = list(size=3)))+
#   xlim(-max(abs(deGSub$logFC)), max(abs(deGSub$logFC))) +
#   theme_bw() +
#   theme(legend.position = "bottom",
#         text = element_text(family = "Arial", size = 18, face = 2), 
#         aspect.ratio=1)
# 
# ggsave('3.Output/CellFigure/volcanoClst1Act_labeled.png', 
#        width = 8, height = 5, units = "in", dpi = 300, bg = "white")
# 
# 
# 
# 
# deGSub %>% 
#   ggplot(aes(x = logFC, 
#              y = sig, 
#              col = ifelse(logFC > 0, 'UP', 'DN'))) +
#   geom_point(alpha = 0.3, size = 3) +
#   geom_label_repel(data = deGSub %>% 
#                      # filter(!feature %in% hg) %>% 
#                      group_by(logFC > 0) %>% 
#                      top_n(5, abs(logFC)),
#                    aes(label = feature), 
#                    size = 6, fontface = 2, alpha = 0.8,family = "Arial", 
#                    show.legend = F) +
#   scale_y_log10() +
#   xlim(-max(abs(deGSub$logFC)), max(abs(deGSub$logFC))) +
#   scale_color_manual(breaks = c('DN', 'UP'),
#                      values = c('blue', 'red'),
#                      labels = c('Down-regulated genes',
#                                 'Up-regulated genes')) +
#   scale_alpha_manual(breaks = c('hl', 'nm'),
#                      values = c(1, 0.5)) +
#   labs(x = '',
#        y = '',
#        col = '') +
#   guides(color = guide_legend(override.aes = list(size=3)))+
#   theme_bw() +
#   theme(legend.position = "none",
#         axis.text = element_blank(), 
#         aspect.ratio=1)
# 
# ggsave('3.Output/CellFigure/volcanoClst1Act_Nolabel.png', 
#        width = 5, height = 5, units = "in", dpi = 300, bg = "white")
# 
# 
# 
# ## cluster 4
# deGSub <- deG %>% 
#   filter(group == "Exhuasted", 
#          padj < 0.05)
# 
# 
# deGSub %>% 
#   ggplot(aes(x = logFC, 
#              y = ifelse(is.infinite(-log(padj)), 
#                         max(-log(padj)[!is.infinite(-log(padj))]) + runif(length(-log(padj)[is.infinite(-log(padj))])), 
#                         -log(padj)), 
#              col = ifelse(logFC > 0, 'UP', 'DN'))) +
#   geom_point(alpha = 0.3, size = 3) +
#   geom_label_repel(data = deGSub %>% 
#                      # filter(!feature %in% hg) %>% 
#                      group_by(logFC > 0) %>% 
#                      top_n(5, abs(logFC)),
#                    aes(label = feature), 
#                    size = 5, fontface = 2, alpha = 0.8,
#                    show.legend = F) +
#   scale_y_log10() +
#   scale_color_manual(breaks = c('DN', 'UP'),
#                      values = c('blue', 'red'),
#                      labels = c('Down-regulated genes',
#                                 'Up-regulated genes')) +
#   labs(x = bquote(bold(log^{fold~change})),
#        y = expression(bold(atop(Significance, (-log^{adjusted~italic(p)~value})))),
#        col = 'Regulation \nstatus') +
#   guides(color = guide_legend(override.aes = list(size=3)))+
#   theme_bw() +
#   theme(legend.position = "bottom",
#         text = element_text(family = "Arial", size = 18, face = 2))
# 
# ggsave('3.Output/Revision/volcanoClst4Exh.png', 
#        width = 8, height = 5, units = "in", dpi = 300, bg = "white")
