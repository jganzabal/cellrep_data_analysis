# Get top markers for each cluster

# %%Load Packages -----------------------------------------------------------

library(Seurat)
library(SeuratDisk)
library(tidyverse)


# Import Seurat Objects ---------------------------------------------------

stim <- LoadH5Seurat("data/HuTcellsCRISPRaPerturbSeq_Re-stimulated.h5Seurat")
resting <- LoadH5Seurat("data/HuTcellsCRISPRaPerturbSeq_Resting.h5Seurat")


# Generate tibble with list column for testing ----------------------------

dat <- tibble(
  condition = c("Resting", "Stim"),
  object = list(resting, stim)
)
dat

rm(resting, stim) # Free up memory


# Get Cluster Markers -----------------------------------------------------

dat <- dat %>%
  mutate(markers = map(object, FindAllMarkers))

dat

dat %>%
  select(-object) %>%
  unnest(cols = markers)
  
# Save markers
dat %>%
  select(-object) %>%
  unnest(cols = markers) %>%
  mutate(condition = if_else(condition == "Stim", "Re-stimulated", condition)) %>% 
  write_tsv("data_out/clusterMarkers.txt")


