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





