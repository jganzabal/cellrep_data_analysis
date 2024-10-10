# Generate log2(CD4/CD8) expression score for each cell


# Load Packages -----------------------------------------------------------

library(Seurat)
library(SeuratDisk)
library(tidyverse)
library(patchwork)
# Install the package if it's not already installed
if (!requireNamespace("scico", quietly = TRUE)) {
  install.packages("scico")
}
library(scico)


# Import Seurat Object ----------------------------------------------------

tcells <- LoadH5Seurat("data_out/tcells_all_filt_05.h5Seurat")


# Get CD4/CD8 score for each cell -----------------------------------------

# Extract dataframe of CD4 and CD8A/CD8B gene expression
CD4_8_df <- FetchData(tcells, c("CD4", "CD8A", "CD8B"))
CD4_8_df %>% head()

# For CD8 average CD8A and CD8B expression for each cell
CD4_8_df <- CD4_8_df %>%
  mutate(CD8 = select(CD4_8_df, CD8A, CD8B) %>% rowMeans())
CD4_8_df %>% head()

CD4_8_df$CD4[CD4_8_df$CD4 != 0] %>% min()
# Min CD4 value != 0 is 0.69

CD4_8_df$CD8[CD4_8_df$CD8 != 0] %>% min()
# Min CD8 value != 0 is 0.35

# Calculate CD4/CD8 score (log2FC). Use minimum value/2 as pseudocount
pseudocount <- CD4_8_df$CD8[CD4_8_df$CD8 != 0] %>% 
  min() %>%
  (function(x) x/2)

CD4_8_df <- CD4_8_df %>%
  mutate(CD4 = CD4 + pseudocount,
         CD8 = CD8 + pseudocount,
         `CD4.CD8.Score` = log2(CD4/CD8))
CD4_8_df %>% head()

# Add to CD4.CD8.Score to seurat object metadata
tcells <- Seurat::AddMetaData(tcells, 
                    select(CD4_8_df, CD4.CD8.Score))
tcells@meta.data %>% head()

# Plot CD4/CD8 Scores --------------------------------------------------

# Distribution of scores
tcells@meta.data %>%
  ggplot(aes(x = CD4.CD8.Score)) +
  geom_density() +
  facet_grid(condition~donor)


# Set cutoffs for CD4/CD8 assignment
tcells@meta.data %>%
  ggplot(aes(x = CD4.CD8.Score)) +
  geom_density() +
  facet_grid(~condition) +
  geom_vline(xintercept = -0.9, linetype = "dashed") +
  geom_vline(xintercept = 1.4, linetype = "dashed") +
  theme_bw() +
  theme(panel.grid = element_blank())


# Add CD4/CD8 assignments to metadata -------------------------------------

tcells@meta.data <- tcells@meta.data %>%
  mutate(CD4_or_CD8 = if_else(CD4.CD8.Score > 1.4, "CD4", "unassigned"),
         CD4_or_CD8 = if_else(CD4.CD8.Score < -0.9, "CD8", CD4_or_CD8))

tcells %>% head()


# Save Seurat Object ------------------------------------------------------

SaveH5Seurat(
  tcells,
  filename = "data_out/tcells_all_filt_06.h5Seurat",
  overwrite = TRUE)

