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


# Calculo de cell activation scores:

1. Differential expression Analysis para NO-TARGET
Group: 'condition' (Stim vs Nostim)
Resultado:
GEN --> LFC es metrica de expresión en el grupo Stim
Se queda con los genes relevantes: pvalues > 0.001
A esto lo defino como **gene_weights**
Es "lo que expresan" los genes cuando estimulo sin target relativo a Nostim

2. Get average expression values for these genes in NT control cells
De matriz scRNA-seq Stim despues de aplicar SCTransform (normalizar):
- Para cada GEN en **gene_weights**, calcula el promedio de counts tomando los celulas NT
- Los pasa a log
- los llama **NT_AVG**
Promedio de expresión por gen en Stim

3. Cell Activation Score
- Para cada célula calculo sobre cada gen:
value * gene_weights / NT_AVG
Donde value es el valor de la matriz
cell_act_score = sumo todos los genes para cada célula
Es como que multiplico al LFC por el "valor" normalizado
A nivel unidades me quedan unidades de LFC





