Type: Package

Title: Identification of lncRNA-disease association using bi-random walks
====================

Files:

1.Dataset

disease_Name.txt：The name of all the diseases

lncRNA_Name.txt：The name of all the lncRNA-disease

interMatrix.mat：Known lncRNA-disease associations

disSim_Jaccard.mat：Disease similairty calculate by gene ontology annotations information

2.Code

BiWalkLDA.m: BiwalkLDA framework to predict potential  lncRNA-disease association

birandom.m: bi-random walk algorithm

Similarity.m: Gaussian interaction profile kernel similarity for disease and lncRNA

normFun.m: Laplacian normalization
