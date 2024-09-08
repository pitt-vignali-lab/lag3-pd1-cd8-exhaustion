# LAG-3 and PD-1 synergy in CD8 T cell exhaustion

Analysis code for:

> **LAG-3 and PD-1 synergize on CD8<sup>+</sup> T cells to drive T cell
> exhaustion and hinder autocrine IFN-γ-dependent anti-tumor immunity**
> Andrews L.P., Butler S.C., *et al.* *Cell* **187**, 4355–4372 (2024).
> PMID [39121848](https://pubmed.ncbi.nlm.nih.gov/39121848/)
> · PMC[11323044](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11323044/)

This repository contains the single-cell RNA-seq, scTCR-seq and RNA velocity
analysis code behind **Figures 2–6** and the related supplemental figures. Tumour
growth, adoptive transfer, flow cytometry and knockdown experiments are
described in the paper and have no code here.

---

## What the analysis does

Tumour-infiltrating CD8<sup>+</sup> T cells were sorted from mice carrying
conditional deletions of *Pdcd1* (PD-1), *Lag3*, or both, in two settings — a
polyclonal compartment and antigen-specific pMEL cells — and profiled by 10x
3′ scRNA-seq with CITE-seq sample multiplexing and paired scTCR-seq.

- Both datasets are integrated by reciprocal PCA and clustered into functional
  states, then re-clustered to resolve substructure (Fig 2).
- RNA velocity shows the direction of state transitions in the overall,
  polyclonal and pMEL compartments (Fig 2e–g).
- Cluster abundance is tested on the k-NN graph with `miloR`, and pathway
  activity with GESECA.
- *Klr* family expression, in particular *Klrc1* (NKG2A), is related to clonal
  behaviour (Fig 3).
- Paired TCR data from CellRanger `vdj` and TRUST4 are combined to measure
  repertoire diversity, clonotype abundance, and clonotype sharing between
  clusters, and clonal cells are scored for proliferation, cytotoxicity and
  apoptosis programmes (Fig 4, Fig 5).

---

## Repository layout

```
1.Data/     input data (not distributed — see "Data availability")
2.Script/   analysis scripts, numbered in figure order
3.Output/   figures and tables produced by the scripts (not distributed)
```

Scripts expect the **repository root** as the R working directory:

```r
setwd("/path/to/lag3-pd1-cd8-exhaustion")
source("2.Script/01_data_integration.R")
```

Each R script carries a header giving its figure panel, inputs, outputs, the
upstream script it depends on, and its original filename in the lab archive.
The four velocity notebooks are Python; their outputs have been cleared so the
repository carries code rather than rendered figures.

---

## Figure → script map

### Integration, clustering and annotation

| Figure | Analysis | Script |
| --- | --- | --- |
| — | Reciprocal PCA integration of the polyclonal and pMEL datasets | [`01_data_integration.R`](2.Script/01_data_integration.R) |
| — | Cluster identities and ranked marker genes | [`02_cluster_labeling.R`](2.Script/02_cluster_labeling.R) |
| **2** | Re-clustering of the clusters of interest | [`03_subclustering.R`](2.Script/03_subclustering.R) |
| — | Functional-state annotation | [`04_cluster_annotation.R`](2.Script/04_cluster_annotation.R) |
| — | Marker gene expression across clusters | [`05_cluster_markers.R`](2.Script/05_cluster_markers.R) |
| **2** | Final cluster labels and marker panels | [`06_fig2_cluster_labels_markers.R`](2.Script/06_fig2_cluster_labels_markers.R) |
| **2** | Per-cluster volcano plots | [`07_fig2_de_volcano.R`](2.Script/07_fig2_de_volcano.R) |
| **2** | Key features distinguishing the genotypes | [`08_fig2_similarity_key_features.R`](2.Script/08_fig2_similarity_key_features.R) |
| **2** | PCA and k-means of sample similarity | [`09_fig2_similarity_pca.R`](2.Script/09_fig2_similarity_pca.R) |
| — | Differential expression behind the similarity structure | [`10_similarity_de_function.R`](2.Script/10_similarity_de_function.R) |
| — | The same comparison by genotype | [`11_similarity_de_by_genotype.R`](2.Script/11_similarity_de_by_genotype.R) |
| — | Pseudobulk differential expression | [`12_pseudobulk_degs.R`](2.Script/12_pseudobulk_degs.R) |
| — | GSEA of the pseudobulk comparison | [`13_gsea_pseudobulk.R`](2.Script/13_gsea_pseudobulk.R) |
| **2** | GESECA gene set co-regulation analysis | [`14_geseca_functional_enrichment.R`](2.Script/14_geseca_functional_enrichment.R) |
| **2** | GESECA follow-up and simplified GO view | [`15_geseca_functional_enrichment_2.R`](2.Script/15_geseca_functional_enrichment_2.R) |
| **2** | Differential abundance on the k-NN graph (`miloR`) | [`16_milo_differential_abundance.R`](2.Script/16_milo_differential_abundance.R) |

### RNA velocity

| Figure | Analysis | Script |
| --- | --- | --- |
| **2e–g** | Prepare loom input and export embeddings | [`17_velocity_prepare_loom.R`](2.Script/17_velocity_prepare_loom.R) |
| **2e** | Velocity of all clusters | [`18_fig2e_velocity_overall.ipynb`](2.Script/18_fig2e_velocity_overall.ipynb) |
| **2f** | Velocity of the polyclonal compartment | [`19_fig2f_velocity_polyclonal.ipynb`](2.Script/19_fig2f_velocity_polyclonal.ipynb) |
| **2g** | Velocity of antigen-specific pMEL cells | [`20_fig2g_velocity_pmel.ipynb`](2.Script/20_fig2g_velocity_pmel.ipynb) |
| **2** | Velocity on the combined object | [`21_velocity_combined.ipynb`](2.Script/21_velocity_combined.ipynb) |

### NKG2A / Klr family (Figure 3)

| Figure | Analysis | Script |
| --- | --- | --- |
| **3** | *Klr* family expression, including *Klrc1* (NKG2A) | [`22_fig3_klr_family_expression.R`](2.Script/22_fig3_klr_family_expression.R) |
| **3** | Correlation of clonal features with *Klrc1* | [`23_fig3_klrc1_correlation.R`](2.Script/23_fig3_klrc1_correlation.R) |
| **3** | Differential expression in the exhausted cluster by genotype | [`24_fig3_cluster4_analysis.R`](2.Script/24_fig3_cluster4_analysis.R) |
| — | Differential expression in clusters 1 and 3 | [`25_de_genes_cluster13.R`](2.Script/25_de_genes_cluster13.R) |

### TCR repertoire and clonal fate (Figures 4–5)

| Figure | Analysis | Script |
| --- | --- | --- |
| **4** | Merge CellRanger `vdj` and TRUST4 TCR calls | [`26_fig4_tcr_init.R`](2.Script/26_fig4_tcr_init.R) |
| **4** | Clonotype calling from TCRβ CDR3 sequences | [`27_fig4_tcr_clonotypes.R`](2.Script/27_fig4_tcr_clonotypes.R) |
| **4** | Clonotype distribution across clusters | [`28_fig4_tcr_distribution.R`](2.Script/28_fig4_tcr_distribution.R) |
| **4** | Combined TCR table across genotypes | [`29_fig4_tcr_combined.R`](2.Script/29_fig4_tcr_combined.R) |
| **4** | Repertoire summary statistics | [`30_fig4_tcr_summary.R`](2.Script/30_fig4_tcr_summary.R) |
| **4** | Clonotype abundance and diversity | [`31_fig4_tcr_abundance.R`](2.Script/31_fig4_tcr_abundance.R) |
| **5** | Proliferation and cytotoxicity module scores | [`32_fig5_tcr_proliferation_score.R`](2.Script/32_fig5_tcr_proliferation_score.R) |
| **5** | Apoptosis module scores | [`33_fig5_tcr_apoptosis_score.R`](2.Script/33_fig5_tcr_apoptosis_score.R) |
| **5** | Clonotype sharing between clusters (graph and heat map) | [`34_fig5_clonality_sharing_matrix.R`](2.Script/34_fig5_clonality_sharing_matrix.R) |

### Companion dataset and figure rendering

| Figure | Analysis | Script |
| --- | --- | --- |
| **6** | Clustering and density views of the companion CD8 TIL dataset | [`35_fig6_bulk_rnaseq_analysis.R`](2.Script/35_fig6_bulk_rnaseq_analysis.R) |
| **2a–f** | Published Figure 2 UMAP, bar and marker panels | [`36_figure_panels_fig2.R`](2.Script/36_figure_panels_fig2.R) |
| **2** | Published sub-clustering panels | [`37_figure_panels_subclustering.R`](2.Script/37_figure_panels_subclustering.R) |
| **2** | Published volcano panels | [`38_figure_panels_volcano.R`](2.Script/38_figure_panels_volcano.R) |
| **2a–d** | Labelled and unlabelled UMAP and distribution bars | [`39_figure_panels_umap_and_bars.R`](2.Script/39_figure_panels_umap_and_bars.R) |

Run `01` → `02` → `03` → `04` first. `04` writes the annotated object that most
downstream scripts read, and `26`/`29` build the combined TCR table that
Figures 3–5 depend on.

---

## Data availability

**No sequencing data, loom files or intermediate objects are included here.**

| Dataset | Source |
| --- | --- |
| scRNA-seq and scTCR-seq of tumour-infiltrating CD8<sup>+</sup> T cells | GEO [**GSE262740**](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE262740) |
| Companion CD8 TIL RNA-seq | GEO [**GSE271184**](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE271184) |

### Upstream processing happens before these scripts

- **Alignment and counting** — CellRanger 6.1.1 against mm10, with
  `cite-seq-count` 1.4.4 for the TotalSeq-C hashtags. Cells with fewer than 200
  or more than 10,000 UMIs were removed, and genes were kept only if detected at
  ≥ 1 UMI in 1 % of cells.
- **Sample demultiplexing** — hashtag counts were CLR-normalised and resolved
  with `HTODemux` (`positive.quantile = 0.999`), with k-means cutoffs per
  sample; double-positive cells were dropped.
- **TCR** — CellRanger `vdj` plus TRUST4 1.0.12; the two call sets were merged
  and duplicate sequences trimmed.
- **Velocity** — spliced / unspliced loom files were produced with Velocyto
  0.17 and are expected at `1.Data/Velo/` and `1.Data/combined.loom`.

---

## Software

| Tool | Version |
| --- | --- |
| CellRanger | 6.1.1 |
| cite-seq-count | 1.4.4 |
| Velocyto | 0.17 |
| TRUST4 | 1.0.12 |
| R | 4.1.2 |
| Seurat | 4 |
| tidyverse | 2.0.0 |
| presto (`wilcoxauc`) | 1.0.0 |
| fgsea (GESECA) | 1.28.0 |
| msigdbr | 7.5.1 |
| miloR | 1.10.0 |
| SingleCellExperiment | 1.24.0 |
| scater | 1.30.1 |
| scran | 1.30.2 |
| igraph | 2.0.1 |
| ggraph | 2.1.0 |

Differentially expressed genes were called at adjusted *P* < 0.05.

---

## Contact

**Jian Cui** — Senior Bioinformatics Scientist  
Vignali Lab, Department of Immunology  
School of Medicine, University of Pittsburgh  
<cuij@pitt.edu>

Questions about the analysis or this code are welcome.

## License

[MIT](LICENSE)
