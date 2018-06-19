Type: Package
Title: Prediction of lncRNA-disease associations based on Bi-random walk
=================
Files:
1.Dataset
disease_Name.txt store the name of all the diseases
lncRNA_Name.txt store the name of all the lncRNA-disease
interMatrix.mat store known lncRNA-disease associations
disSim_Jaccard.mat sore disease similairty calculate by gene ontology annotations information

2.Code
BiWalkLDA.m: BiwalkLDA framework to predict potential  lncRNA-disease association
birandom.m: bi-random walk algorithm
Similarity.m: Gaussian interaction profile kernel similarity for disease and lncRNA
normFun.m: Laplacian normalization
