# Generate 2D contour plots of sgRNA distributions in UMAP space

# %% Load Packages -----------------------------------------------------------

library(Seurat)
library(SeuratDisk)
library(tidyverse)
library(patchwork)


# Import Seurat Objects ---------------------------------------------------

dat <- tibble(
  condition = c("Resting", "Re-stimulated"),
  object = list(
    LoadH5Seurat("data/HuTcellsCRISPRaPerturbSeq_Resting.h5Seurat"),
    LoadH5Seurat("data/HuTcellsCRISPRaPerturbSeq_Re-stimulated.h5Seurat")
  )
)
dat


# Re-stim Contour plots ----------------------------------------------------

n_bins <- 12
plot_limits_x <- c(-8.5, 9.5)
plot_limits_y <- c(-5.5, 6.5)


# Contour plot of all cells (re-stimulated)
dat$object[[2]] %>%
  Embeddings("umap") %>%
  as.data.frame() %>%
  ggplot(aes(x = UMAP_1, y = UMAP_2)) +
  geom_density_2d(aes(color = ..level..),
                  bins = n_bins,
                  size = 0.8,
                  show.legend = F) +
  scale_color_viridis_c() +
  expand_limits(x = plot_limits_x, y = plot_limits_y) +
  theme_bw()


# Function to generate a contour plot for a given gene or category
make_contour_plot <- function(obj, 
                              Gene = NULL, 
                              category = NULL, 
                              control_backdrop = NULL) {
  
  # Subset data by gene or gene category
  if(!is.null(Gene)) {
    dat <- subset(obj, subset = gene == Gene) 
  } 
  if(!is.null(category)) {
    dat <- subset(obj, subset = gene_category == category) 
  }
  
  # Build plot w control layer
  if(!is.null(control_backdrop)) {
    plot <- control_backdrop +
      geom_density_2d(aes(color = ..level..),
                      data =  dat %>%
                        Embeddings("umap") %>%
                        as.data.frame(),
                      bins = n_bins, 
                      size = 0.8,
                      show.legend = FALSE) +
      scale_color_viridis_c() +
      ggtitle(label = paste0(Gene,category))
    return(plot)
  } 
  
  # Build plot without control layer
  if(is.null(control_backdrop)) {
    plot <- dat %>%
      Embeddings("umap") %>%
      as.data.frame() %>%
      ggplot(aes(x = UMAP_1, y = UMAP_2)) +
      geom_density_2d(aes(color = ..level..), bins = n_bins, size = 0.8,
                      show.legend = FALSE) +
      scale_color_viridis_c() +
      expand_limits(x = plot_limits_x, y = plot_limits_y) +
      theme_bw() +
      ggtitle(label = paste0(Gene,category))
    return(plot)
  }
  
}

# Generate contour of only NT control cells for backdrop comparison
control_contour <- subset(dat$object[[2]], gene == "NO-TARGET") %>%
  Embeddings("umap") %>%
  as.data.frame() %>%
  ggplot(aes(x = UMAP_1, y = UMAP_2)) +
  stat_density2d(aes(fill = ..level..), 
                 geom = "polygon", 
                 n = 100,
                 bins = n_bins,
                 show.legend = FALSE) +
  scale_fill_gradient(low = "#f0f0f0", high = "#969696") +
  expand_limits(x = plot_limits_x, y = plot_limits_y) +
  theme_bw()
control_contour

# Test function
make_contour_plot(dat$object[[2]], "MAP4K1")
make_contour_plot(dat$object[[2]], "MAP4K1", control_backdrop = control_contour)


# Apply function to all sgRNA targets
gene_contour_plots <- dat$object[[2]]@meta.data$gene %>%
  as.character() %>%
  unique() %>%
  as.list() %>%
  map(function(x) make_contour_plot(dat$object[[2]], 
                                    Gene = x,
                                    control_backdrop = control_contour))


# Plots for manuscript main figure

ALL_CELLS <- control_contour +
  geom_density_2d(aes(color = ..level..),
                  data =  dat$object[[2]] %>%
                    Embeddings("umap") %>%
                    as.data.frame(),
                  bins = n_bins, 
                  size = 0.8,
                  show.legend = FALSE) +
  scale_color_viridis_c() +
  ggtitle("All Cells")
ALL_CELLS

CONTROL <- make_contour_plot(dat$object[[2]], 
                             Gene = "NO-TARGET",
                             control_backdrop = control_contour) +
  ggtitle("Non-Targeting Controls")
CONTROL

PERTURBED <- control_contour +
  geom_density_2d(aes(color = ..level..),
                  data =  dat$object[[2]] %>%
                    subset(gene != "NO-TARGET") %>%
                    Embeddings("umap") %>%
                    as.data.frame(),
                  bins = n_bins, 
                  size = 0.8,
                  show.legend = FALSE) +
  scale_color_viridis_c() +
  ggtitle("Perturbed Cells")
PERTURBED


VAV1 <-make_contour_plot(dat$object[[2]], 
                         Gene = "VAV1",
                         control_backdrop = control_contour)
MAP4K1 <-make_contour_plot(dat$object[[2]], 
                           Gene = "MAP4K1",
                           control_backdrop = control_contour)
FOXQ1 <- make_contour_plot(dat$object[[2]], 
                           Gene = "FOXQ1",
                           control_backdrop = control_contour)
GATA3 <- make_contour_plot(dat$object[[2]], 
                           Gene = "GATA3",
                           control_backdrop = control_contour)
IL1R1 <- make_contour_plot(dat$object[[2]], 
                           Gene = "IL1R1",
                           control_backdrop = control_contour)
TNFRSF1A <- make_contour_plot(dat$object[[2]], 
                              Gene = "TNFRSF1A",
                              control_backdrop = control_contour)
TBX21 <- make_contour_plot(dat$object[[2]], 
                           Gene = "TBX21",
                           control_backdrop = control_contour)


CONTROL + PERTURBED + MAP4K1 + 
  VAV1 + TNFRSF1A + IL1R1 +
  FOXQ1 + GATA3 + TBX21 +
  plot_layout(ncol = 3) &
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) 



# Resting cells -----------------------------------------------------------

# Contour plots -----------------------------------------------------------

n_bins <- 12
plot_limits_x <- c(-7.3, 10.5)
plot_limits_y <- c(-7.5, 7.8) 


# Contour plot of all cells
dat$object[[1]] %>%
  Embeddings("umap") %>%
  as.data.frame() %>%
  ggplot(aes(x = UMAP_1, y = UMAP_2)) +
  geom_density_2d(aes(color = ..level..), 
                  bins = n_bins, 
                  size = 0.8,
                  show.legend = F) +
  scale_color_viridis_c() +
  expand_limits(x = plot_limits_x, y = plot_limits_y) +
  theme_bw()


# Backdrop control contour
control_contour <- subset(dat$object[[1]], gene == "NO-TARGET") %>%
  Embeddings("umap") %>%
  as.data.frame() %>%
  ggplot(aes(x = UMAP_1, y = UMAP_2)) +
  stat_density2d(aes(fill = ..level..), 
                 geom = "polygon", 
                 n = 100,
                 bins = n_bins,
                 show.legend = FALSE) +
  scale_fill_gradient(low = "#f0f0f0", high = "#969696") +
  expand_limits(x = plot_limits_x, y = plot_limits_y) +
  theme_bw()
control_contour


# Apply function to guide genes
gene_contour_plots <- dat$object[[1]]@meta.data$gene %>%
  as.character() %>%
  unique() %>%
  as.list() %>%
  map(function(x) make_contour_plot(dat$object[[1]], 
                                    Gene = x,
                                    control_backdrop = control_contour))


# Plots for manuscript main figure


ALL_CELLS <- control_contour +
  geom_density_2d(aes(color = ..level..),
                  data =  dat$object[[1]] %>%
                    Embeddings("umap") %>%
                    as.data.frame(),
                  bins = n_bins, 
                  size = 0.8,
                  show.legend = FALSE) +
  scale_color_viridis_c() +
  ggtitle("All Cells")
ALL_CELLS

CONTROL <- make_contour_plot(dat$object[[1]], 
                             Gene = "NO-TARGET",
                             control_backdrop = control_contour) +
  ggtitle("Non-Targeting Controls")
CONTROL

PERTURBED <- control_contour +
  geom_density_2d(aes(color = ..level..),
                  data =  dat$object[[1]] %>%
                    subset(gene != "NO-TARGET") %>%
                    Embeddings("umap") %>%
                    as.data.frame(),
                  bins = n_bins, 
                  size = 0.8,
                  show.legend = FALSE) +
  scale_color_viridis_c() +
  ggtitle("Perturbed Cells")
PERTURBED


VAV1 <-make_contour_plot(dat$object[[1]], 
                         Gene = "VAV1",
                         control_backdrop = control_contour)
MAP4K1 <-make_contour_plot(dat$object[[1]], 
                           Gene = "MAP4K1",
                           control_backdrop = control_contour)
FOXQ1 <- make_contour_plot(dat$object[[1]], 
                           Gene = "FOXQ1",
                           control_backdrop = control_contour)
GATA3 <- make_contour_plot(dat$object[[1]], 
                           Gene = "GATA3",
                           control_backdrop = control_contour)
IL1R1 <- make_contour_plot(dat$object[[1]], 
                           Gene = "IL1R1",
                           control_backdrop = control_contour)
TNFRSF1A <- make_contour_plot(dat$object[[1]], 
                              Gene = "TNFRSF1A",
                              control_backdrop = control_contour)
TBX21 <- make_contour_plot(dat$object[[1]], 
                           Gene = "TBX21",
                           control_backdrop = control_contour)


CONTROL + PERTURBED + MAP4K1 + 
  VAV1 + TNFRSF1A + IL1R1 +
  FOXQ1 + GATA3 + TBX21 +
  plot_layout(ncol = 3) &
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) 




