# Generate UMAP plot for all cells (both restim/resting conditions)

# Load Packages -----------------------------------------------------------

library(Seurat)
library(SeuratDisk)
library(tidyverse)
library(patchwork)
library(ggrastr)


# Import Seurat object ----------------------------------------------------

tcells <- LoadH5Seurat("data_out/tcells_all_filt_06.h5Seurat")
tcells


# UMAP and Cluster --------------------------------------------------------

tcells <- RunPCA(tcells, npcs = 100)
tcells <- RunUMAP(tcells, dims = 1:35)
tcells <- FindNeighbors(tcells, dims = 1:35)
tcells <- FindClusters(tcells, resolution = 1.2, algorithm = 3)


# Rasterized condition UMAP plot for manuscript ---------------------------

umap_df <- tcells@reductions$umap@cell.embeddings %>%
  as.data.frame() %>%
  rownames_to_column("cell") %>%
  tibble()

umap_df <- full_join(umap_df,
                     tcells@meta.data %>%
                       as.data.frame() %>%
                       rownames_to_column("cell"),
                     by = "cell")
umap_df

set.seed(1)
umap_plot <- umap_df %>%
  sample_frac() %>%
  ggplot(aes(x = UMAP_1 * -1, y = UMAP_2, color = condition)) +
  rasterize(
    geom_point(shape = 16, size = 0.1, alpha = 0.5),
    dpi = 300) +
  scale_color_brewer(palette = "Set1", direction = -1,
                     name = "Condition") + 
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "top") +
  xlab("UMAP 1") +
  ylab("UMAP 2")
umap_plot




