# CRISPR activation and interference screens decode stimulation responses in primary human T cells
## Results replication
- Coding replication from raw data to some figures and tables
- Mainly in python

# Data
### Experiment 1: 
#### Genome-wide CRISPRa screens identify regulators of IL-2 and IFN-γ production in T cells
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE174292:
- Raw count (IL2, IFNg): https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE174255

### Experiment 2: 
#### Complementary CRISPRa and CRISPRi screens comprehensively reveal circuits of cytokine production in T cells:
- Raw count CD4 Supplementary(IL2, IFNg, TFNa): https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE190846

### Experiment 3: 
#### Arrayed characterization of selected CRISPRa screen hits
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE174292:
- Raw count (4 donors stim vs no stim, FOXQ1): https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE174284


### Experiment 4: 
#### CRISPRa Perturb-seq characterizes the molecular phenotypes of cytokine regulators
- Raw count: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE190604
- R Repository and processed Data: https://zenodo.org/records/5784651

# Folders
## experiment_4/R code
You can run the original R code with a devcontainer

## R_translated
R code from zenodo experiment 4 translated to python

# Python Environment
```bash
conda create -n cellrep python=3.7
conda activate cellrep
conda install -c bioconda mageck=0.5.9.2
pip install -r requirements.txt
```


# Cell activation scores in R code explained:

1. Differential expression Analysis for NO-TARGET
Group: 'condition' (Stim vs Nostim)
Result:
GEN --> LFC is the expression metric in Stim gropu
Keep only relevant genes: pvalues > 0.001
This is called **gene_weights**: It is what genes express when they are stimulated relative to NoStim (Only NO-TARGET)

2. Get average expression values for these genes in NT control cells
From scRNA-seq Stim matix after SCTransform (normalization):
- For each gene in **gene_weights**, calculate average counts from only NT cells
- To log
- The are called **NT_AVG**
Expression average per gene in Stim

3. Cell Activation Score
- For each cell we calculate it for every gene
- value * gene_weights / NT_AVG
Where value is the count in the matrix (expression)
cell_act_score = Add all genes for each cell

It is like multiplying the LFC for the "value" normalized
The "physical units" are in LFC





