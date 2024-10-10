# UMAP and cluster for downstream analysis
# Picked dimensions and clustering algorithms by examining outputs of 
# various numbers of dimensions and clustering algorithms

# NoStim and Stim will be split for seperate analyses

# Load Packages -----------------------------------------------------------

library(Seurat)
library(SeuratDisk)
library(scico)
library(tidyverse)

# Import Seurat Object -----------------------------------------------------

tcells <- LoadH5Seurat("data_out/tcells_all_filt_06.h5Seurat")

# Split into two object based on conditions -------------------------------

nostim <- subset(tcells, condition == "Nostim")
stim <- subset(tcells, condition == "Stim")
rm(tcells)

# UMAP and Cluster --------------------------------------------------------

nostim <- RunPCA(nostim, npcs = 50, ndims.print = 1:20)
nostim <- RunUMAP(nostim, dims = 1:20)
nostim <- FindNeighbors(nostim, dims = 1:20)
nostim <- FindClusters(nostim, resolution = 0.5, algorithm = 3)

stim <- RunPCA(stim, npcs = 50, ndims.print = 1:20)
stim <- RunUMAP(stim, dims = 1:20)
stim <- FindNeighbors(stim, dims = 1:20)
stim <- FindClusters(stim, resolution = 0.4, algorithm = 3)


# Plots -------------------------------------------------------------------

UMAPPlot(nostim, label = T)
UMAPPlot(stim, label = T)

# Combine similar clusters ------------------------------------------------

# Combine Stim clusters 2 and 15. In downstream analysis it was shown that
# these two clusters were very similar. Cluster 2 represents the bulk of cells
# with negative regulator guides, and cluster 15 is exclusively MUC1 guides

# Thread on merging clusters https://github.com/satijalab/seurat/issues/3202

new.cluster.ids <- c(0:14,2) %>% as.character()
names(new.cluster.ids) <- levels(stim)
new.cluster.ids
stim <- RenameIdents(stim, new.cluster.ids)
UMAPPlot(stim, label = TRUE)

# Join "15" and "2" in the metadata
stim@meta.data$seurat_clusters <- fct_collapse(stim@meta.data$seurat_clusters, 
                                               `2` = c("2","15"))

# Save data ---------------------------------------------------------------

SaveH5Seurat(
  nostim, 
  filename = "data_out/Nostim.h5Seurat")

SaveH5Seurat(
  stim, 
  filename = "data_out/Stim.h5Seurat")


