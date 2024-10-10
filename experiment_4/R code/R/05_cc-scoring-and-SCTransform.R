# Cell cycle regression and Normalize/scale scRNA-seq matrix with SCTransform

# https://satijalab.org/seurat/articles/sctransform_vignette.html
# https://satijalab.org/seurat/articles/cell_cycle_vignette.html

# Will follow order proposed in this comment for addressing SCTransform and cell 
# cycle scoring: https://github.com/satijalab/seurat/issues/1679#issuecomment-508744726

# Load Packages -----------------------------------------------------------

library(Seurat)
library(SeuratDisk)
library(tidyverse)
library(patchwork)


# Import Seurat Object ----------------------------------------------------

tcells <- LoadH5Seurat("data_out/tcells_all_filt_04.h5Seurat")

# Transform Data ---------------------------------------------------------

# Run SCTransform
tcells <- SCTransform(tcells, 
                      vars.to.regress = "percent.mt", 
                      verbose = TRUE)

# Cell Cycle Scoring ------------------------------------------------------

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

tcells <- CellCycleScoring(tcells, 
                           s.features = s.genes, 
                           g2m.features = g2m.genes, 
                           set.ident = TRUE)


# Run SCTransform: Regress S.Score and G2M.Score --------------------------

tcells <- SCTransform(tcells, 
                      assay = 'RNA',
                      new.assay.name = 'SCT', # Overwrite original SCT
                      vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"))


# Save Seurat Object ------------------------------------------------------

SaveH5Seurat(tcells, 
             filename = "data_out/tcells_all_filt_05.h5Seurat",
             overwrite = T)
