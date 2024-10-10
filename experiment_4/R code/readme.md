# How to test this
Open this folder with cursos or vscode with de devcontainer for dependencies to work easily


what does sc.tl.rank_genes_groups does?

Identifies marker genes (or differentially expressed genes) between different groups of cells.

Rather than telling you the likelihood of a gene belonging to a group, it tells you how much more (or less) expressed a gene is in one group compared to others, which helps in identifying marker genes for that group. These marker genes are important because they are more characteristic of a particular group of cells (e.g., a cluster).

'wilcoxon': Wilcoxon rank-sum test (non-parametric).

The p-value in the context of sc.tl.rank_genes_groups() (and statistical tests in general) indicates the **likelihood** that the observed **difference in gene expression** between groups occurred **by chance**.

A small p-value alone doesn't tell you about the magnitude of the difference, only the likelihood that it's statistically significant. To understand the biological importance of the difference, you would also need to look at metrics like the log fold change, which tells you how large the expression difference is between groups.