# %%
import scanpy as sc
# %%
adata = sc.read_10x_mtx(
    '../data/GSE190604/',
    prefix='GSE190604_',
    cache=True,
    gex_only=False, # Only gene expression
    # make_unique=True
)
# %%
[g for g in adata.var.index if g.startswith('CD4')]

# %%
