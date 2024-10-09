# %%
# # Cell cycle regression and Normalize/scale scRNA-seq matrix
# We'll use scanpy for the Python equivalent of Seurat's operations
# Note: There's no direct equivalent to SCTransform in scanpy, so we'll use standard normalization and scaling

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
# %%
# Import AnnData Object ----------------------------------------------------

tcells = sc.read_h5ad("data/tcells_all_filt.h5ad")
tcells
# %%
tcells.obs
# %%
#  Transform Data ---------------------------------------------------------

# Normalize and log-transform the data
sc.pp.normalize_total(tcells, target_sum=1e4)
sc.pp.log1p(tcells)
# %%
# Regress out effects of total counts per cell and the percentage of mitochondrial genes
sc.pp.regress_out(tcells, ['total_counts', 'pct_counts_mt'])
# %%
tcells
# %%
# Scale the data
sc.pp.scale(tcells, max_value=10)

# %%
#  Cell Cycle Scoring ------------------------------------------------------
# Esto lo agregue yo por que no ecnontre los archivos
s_genes = ['MCM5', 'PCNA', 'TYMS', 'FEN1', 'MCM2', 'MCM4', 'RRM1', 'UNG', 'GINS2', 'MCM6', 'CDCA7', 'DTL', 'PRIM1', 'UHRF1', 'CENPU', 'HELLS', 'RFC2', 'RPA2', 'NASP', 'RAD51AP1', 'GMNN', 'WDR76', 'SLBP', 'CCNE2', 'UBR7', 'POLD3', 'MSH2', 'ATAD2', 'RAD51', 'EXO1', 'TIPIN', 'DSCC1', 'BLM', 'CASP8AP2', 'USP1', 'CLSPN', 'POLA1', 'CHAF1B', 'BRIP1', 'E2F8']
g2m_genes = ['HMGB2', 'CDK1', 'NUSAP1', 'UBE2C', 'BIRC5', 'TPX2', 'TOP2A', 'NDC80', 'CKAP2L', 'CKAP2', 'AURKB', 'BUB1', 'KIF11', 'ANP32E', 'TUBB4B', 'GTSE1', 'KIF20B', 'HJURP', 'CDCA3', 'CDC20', 'TTK', 'CDC25C', 'KIF2C', 'RANGAP1', 'NCAPD2', 'DLGAP5', 'CDCA2', 'NUSAP1', 'UBE2C', 'TPX2', 'TOP2A', 'AURKB', 'BUB1', 'KIF11']

# Load cell cycle genes (you might need to prepare this list separately)
# s_genes = pd.read_csv('path_to_s_genes.csv', header=None)[0].tolist()
# g2m_genes = pd.read_csv('path_to_g2m_genes.csv', header=None)[0].tolist()

# Score cell cycle
sc.tl.score_genes_cell_cycle(tcells, s_genes=s_genes, g2m_genes=g2m_genes)

# %%
# Regress out cell cycle scores --------------------------

sc.pp.regress_out(tcells, ['S_score', 'G2M_score'])
# %%
# Scale the data again after regression
sc.pp.scale(tcells, max_value=10)

# Save AnnData Object ------------------------------------------------------
# %%
tcells.write_h5ad("data/tcells_all_filt_processed.h5ad")
