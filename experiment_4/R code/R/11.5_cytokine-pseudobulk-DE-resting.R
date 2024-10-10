# Pseudobulk differential expression analysis of sgRNA targets versus NT control
# cells for each cytokine (or other effector) gene - Resting condition

# Load Packages -----------------------------------------------------------

library(Seurat)
library(SeuratDisk)
library(tidyverse)
library(patchwork)

# Import Seurat Objects ---------------------------------------------------

resting <- LoadH5Seurat("data/HuTcellsCRISPRaPerturbSeq_Resting.h5Seurat")

# Gene Lists --------------------------------------------------------------

cytokine_genes <- read_tsv("data/GO_0005125_Cytokines.txt", col_names = F)
cytokine_genes <- cytokine_genes$X3 %>% unique() %>% sort()
cytokine_genes

granzyme_genes <- c("GZMA", "GZMB", "GZMH", "GZMK", "GZMM", "PRF1", "GNLY", 
                    "THBS1")
# Why TSP1: https://science.sciencemag.org/content/368/6493/897

detected_features <- resting@assays$SCT@counts %>% rownames()

cytokines_in_data <- cytokine_genes[cytokine_genes %in% detected_features]
cytokines_in_data
cytokine_genes[!cytokine_genes %in% detected_features]

granzyme_genes_in_data <- granzyme_genes[granzyme_genes %in% detected_features]
granzyme_genes_in_data

test_genes <- c(cytokines_in_data, granzyme_genes_in_data)


# Find Differentially Expressed Genes -------------------------------------

# Versus non-targeting control cells, genes in library set

# Get sgRNA target genes from library
targetgenes <- resting@meta.data %>%
  dplyr::filter(gene != "NO-TARGET") %>%
  pull(gene) %>%
  as.character() %>%
  unique()
targetgenes

# Function to get differential expression for target gene of interest versus
# non-targeting control cells (in a given condition)
## On setting new identies for FindMarkers(): 
##  https://github.com/satijalab/seurat/issues/252
target_gene_overexpression <- function(obj, 
                                       Gene, 
                                       GeneSet = test_genes) {
  obj <- subset(obj, gene %in% c("NO-TARGET", Gene))
  Idents(obj) <- obj@meta.data$'crispr'
  markers <- FindMarkers(obj, 
                         ident.1 = "perturbed", 
                         ident.2 = "NT",
                         features = GeneSet,
                         logfc.threshold = 0,
                         min.pct = 0)
  return(markers)
}

# Dataframe for differential expression testing
de_df <- bind_rows(tibble(gene = targetgenes))
de_df


# CD8 DE 
CD8 <- resting %>% subset(CD4.or.CD8 == "CD8")

de_resting_CD8 <- de_df %>%
  mutate(DiffExpr = map(gene,
                        function(x) target_gene_overexpression(CD8, x)))

de_resting_CD8 <- de_resting_CD8 %>%
  mutate(DiffExpr = map(DiffExpr, rownames_to_column)) %>%
  unnest(DiffExpr) %>%
  dplyr::rename(crispr_gene = gene, mRNA_gene = rowname)
de_resting_CD8

# CD4 DE 
CD4 <- resting %>% subset(CD4.or.CD8 == "CD4")

de_resting_CD4 <- de_df %>%
  mutate(DiffExpr = map(gene,
                        function(x) target_gene_overexpression(CD4, x)))

de_resting_CD4 <- de_resting_CD4 %>%
  mutate(DiffExpr = map(DiffExpr, rownames_to_column)) %>%
  unnest(DiffExpr) %>%
  dplyr::rename(crispr_gene = gene, mRNA_gene = rowname)
de_resting_CD4



# Get Cytokine / CRISPRa categories ---------------------------------------

# Add guide and cytokine categories 

guide_categories <- read_tsv("data/guide-target-categories.txt")

# Switch to factors and rearrange levels
guide_categories$gene_functional_category %>% unique()
guide_categories$gene_functional_category <- factor(
  guide_categories$gene_functional_category,
  levels = c("Negative Regulator", "Control", "TCR Stim/Co-Stim", 
             "GTPase signaling", "NF-kB Signaling", "TNF receptor superfamily", 
             "Other Cytokine Receptor", "Other Signal Transduction",
             "Transcription Factor", "Other / Unknown")
)
guide_categories

# Add some line-breaks to categories for more compact plotting
guide_categories$gene_functional_category <- guide_categories$gene_functional_category %>%
  fct_recode("Negative\nRegulator" = "Negative Regulator",
             "TCR\nStim/Co-Stim"= "TCR Stim/Co-Stim",
             "GTPase\nsignaling" = "GTPase signaling",
             "NF-kB\nSignaling" ="NF-kB Signaling",
             "TNF\nreceptor\nsuperfamily" = "TNF receptor superfamily",
             "Other\nCytokine\nReceptor" = "Other Cytokine Receptor",
             "Other\nSignal\nTransduction" = "Other Signal Transduction",
             "Transcription\nFactor" = "Transcription Factor",
             "Other /\nUnknown" = "Other / Unknown")


# Do the same for cytokine categories
cytokine_categories <- read_tsv("data/cytokine_categories.txt")
cytokine_categories
cytokine_categories$Category %>% unique()
cytokine_categories <- cytokine_categories %>%
  dplyr::filter(!is.na(Category))
cytokine_categories

cytokine_categories$Category %>% unique()
cytokine_categories$Category <- factor(
  cytokine_categories$Category,
  levels = c("Th1","Th2","Chemokine","Growth Factor","Apoptosis","Other",
             "Granzyme/Cytolytic")) 
cytokine_categories

cytokine_categories$Category <- cytokine_categories$Category %>%
  fct_recode("Growth\nFactor" = "Growth Factor",
             "Granzyme/\nCytolytic" = "Granzyme/Cytolytic")





# Plots -------------------------------------------------------------------

max.lfc <- 1
min.lfc <- -1
min.padj <- 1e-10

# Add categories to dataframes for plotting
CD4_plot_df <- left_join(rename(cytokine_categories, mRNA_gene = Cytokine),
                         de_resting_CD4,
                         by = "mRNA_gene") %>%
  left_join(rename(guide_categories, crispr_gene = gene), 
            by = "crispr_gene") %>%
  rename(cytokine_category = Category)
CD4_plot_df

CD8_plot_df <- left_join(rename(cytokine_categories, mRNA_gene = Cytokine),
                         de_resting_CD8,
                         by = "mRNA_gene") %>%
  left_join(rename(guide_categories, crispr_gene = gene), 
            by = "crispr_gene") %>%
  rename(cytokine_category = Category)
CD8_plot_df

# Filter out targets with no significant DE expressed genes 
genes2keep_CD4CD8 <- bind_rows(CD4_plot_df, CD8_plot_df) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  pull(crispr_gene) %>%
  unique()
genes2keep_CD4CD8

# Plots
CD4_plot <- CD4_plot_df %>%
  dplyr::filter(crispr_gene %in% genes2keep_CD4CD8) %>%
  mutate(avg_log2FC = if_else(avg_log2FC > max.lfc, max.lfc, avg_log2FC),
         avg_log2FC = if_else(avg_log2FC < min.lfc, min.lfc, avg_log2FC),
         p_val_adj = if_else(p_val_adj < min.padj, min.padj, p_val_adj),
         significant = if_else(p_val_adj < 0.05, TRUE, FALSE)) %>%
  ggplot(aes(x = crispr_gene, y = mRNA_gene, fill = avg_log2FC, 
             size = -log10(p_val_adj))) +
  geom_point(aes(color = significant), shape = 21) +
  scale_color_manual(values = c("white", "black")) +
  facet_grid(gene_functional_category~cytokine_category, 
             scales = "free", space = "free") +
  scale_size_area(max_size = 2.5) +
  scale_fill_distiller(palette = "PiYG", limits = c(min.lfc, max.lfc)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6)) +
  theme(axis.text.y = element_text(size = 6)) +
  theme(strip.text.y = element_text(angle = 0, size = 6)) +
  theme(strip.text.x = element_text(size = 6, angle = 45)) +
  theme(panel.spacing = unit(0, "lines")) +
  xlab("CRISPRa Target") +
  ylab("Cytokine Gene") +
  coord_flip() +
  theme(legend.position = "bottom") +
  ggtitle("CD4+ T cells")


CD8_plot <- CD8_plot_df %>%
  dplyr::filter(crispr_gene %in% genes2keep_CD4CD8) %>%
  mutate(avg_log2FC = if_else(avg_log2FC > max.lfc, max.lfc, avg_log2FC),
         avg_log2FC = if_else(avg_log2FC < min.lfc, min.lfc, avg_log2FC),
         p_val_adj = if_else(p_val_adj < min.padj, min.padj, p_val_adj),
         significant = if_else(p_val_adj < 0.05, TRUE, FALSE)) %>%
  ggplot(aes(x = crispr_gene, y = mRNA_gene, fill = avg_log2FC, 
             size = -log10(p_val_adj))) +
  geom_point(aes(color = significant), shape = 21) +
  scale_color_manual(values = c("white", "black")) +
  facet_grid(gene_functional_category~cytokine_category, 
             scales = "free", space = "free") +
  scale_size_area(max_size = 2.5) +
  scale_fill_distiller(palette = "PiYG", limits = c(min.lfc, max.lfc)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6)) +
  theme(axis.text.y = element_text(size = 6)) +
  theme(strip.text.y = element_text(angle = 0, size = 6)) +
  theme(strip.text.x = element_text(size = 6, angle = 45)) +
  theme(panel.spacing = unit(0, "lines")) +
  xlab("CRISPRa Target") +
  ylab("Cytokine Gene") +
  coord_flip() +
  theme(legend.position = "bottom") +
  ggtitle("CD8+ T cells")

CD4_plot+CD8_plot

# ggsave("cytokines-resting-DE.pdf",
#        width = 8.8,
#        height = 7)
