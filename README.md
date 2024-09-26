# MISF
## 1.Introduction
MISF is a Multimodal data Integration algorithm based on adaptive Similarity network learning and matrix Factorization. MISF adaptively integrates multimodal data and learns the lower-dimensional representations of cells and genes.
## 2.Requirements:
   Matlab2023
## 3.Quick start
### 3.1 Prepare data
The inputs include multi-omics data or spatial transcriptome data. 

    load(mESC.mat)
    X1=RNA;  
    X2=DNA;
    true_label=real_label;
### 3.2 Run MISF
Setting the number of clusters K value.
run "main_MISF.m"
## 4.Notes
Multi-omics data application notes:

1. For algorithm implementation, please run the file “main_MISF.m”.
2. The number of data clusters is set: 3 for scCancer, 2 for mESC, 9 for Pbmc and 14 for Kidney.
3. The RNA modality is labeled as “RNA”, the ATAC modality as “ATAC”, and the DNA modality as “DNA” for all data sets. Different datasets have different modal types, so when introducing the dataset, please note that
When introducing datasets, pay attention to the assignment of X1 and X2, and pay special attention to the assignment of X2 to DNA for mESC.
4. scCancer and mESC datasets have real_label; Pbmc and Kidney datasets don't have real_label, so we should pay attention to the operation of true_label assignment when introducing data.
5. scCancer, mESC two dataset evaluation indicators are used in external clustering evaluation indicators, the results can be obtained directly run; Pbmc, Kidney dataset indicators using internal clustering evaluation indicators, the results need to be saved “F”.
For the Pbmc,Kidney dataset, the indicators use internal clustering evaluation indicators, and the results can be obtained by running the “internal_evaluate.m” file after saving the “F” and “prel” variables. In addition, when applying the Pbmc and Kidney datasets, it is necessary to hide the external evaluation indicator codes in the “main_MISF” file.
In addition, when applying the Pbmc and Kidney datasets, it is necessary to hide the external evaluation index code in the “main_MISF” file to avoid reporting errors.
6. The “r” in the “main_MISF.m” file should be larger than the number of target clusters in the application dataset.
7. The iteration threshold sita of clustering_MISF function should be set according to the demand, generally set to 1.0e-3 or 1.0e-6.

Spatial transcriptome data application notes:

1. Please run the file “main_MISF.m” to realize the algorithm.
2. Set the number of data clusters: 6 for Cortex, 9 for HPOA, and 7 for Prostate.
3. The RNA modal label of all data sets is “gene1_1”, and the spatial coordinate modal label is “spatial1_1”.
4. Cortex and HPOA datasets have real labels “type1_1”; Prostate dataset does not have real labels, so we should pay attention to whether the true_label assignment operation should be hidden when introducing data.
5. Cortex, HPOA two dataset evaluation indicators are used in external clustering evaluation indicators, the results can be obtained directly run; Prostate dataset indicators using internal clustering evaluation indicators, the results need to be saved “F”.
For the Prostate dataset, the indicators use internal clustering evaluation indicators, and the results can be obtained by running the “internal_evaluate.m” file after saving the “F” and “prel” variables. In addition, when applying the Prostate dataset, it is necessary to hide the external evaluation indicator codes in the “main_MISF” file.
file when applying the Prostate dataset.
6. The “r” in line 35 of the “main_MISF.m” file needs to be larger than the number of target clusters in the application dataset.
7. The iteration threshold sita of the clustering_MISF function should be set according to the demand, and is usually set to 1.0e-3 or 1.0e-6.


    result <- scTsI(data_sc,data_bulk)
