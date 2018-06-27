License
====================
PROGRAM: BiwalkLDA

AUTHOR: Jialu Hu and Yiqun Gao

EMAIL: jhu@nwpu.edu.cn, yiqun.gao@nwpu-bioinformatics.com

Copyright (MTLAB) <2018> 

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses/.

Title: Identification of lncRNA-disease association using bi-random walks
====================
###########################################################################################
This section introduces functions of different files. 

1.Dataset

disease_Name.txt：The name of all the diseases

lncRNA_Name.txt：The name of all the lncRNA-disease

interMatrix.mat：Known lncRNA-disease associations

disSim_Jaccard.mat：Disease similairty calculate by gene ontology annotations information

2.Code

BiWalkLDA.m: BiwalkLDA framework to predict potential lncRNA-disease association

birandom.m: Bi-random walk algorithm

Similarity.m: Calculate Gaussian interaction profile kernel similarity for disease and lncRNA

normFun.m: Laplacian normalization
###########################################################################################
This section introduces how to run BiWalkLDA to make prediction

1.Download souce code from https://github.com/screamer/biwalk

2.Open MATLAB and enter the folder ./Code

3.Run command ‘predict_result = BiWalkLDA(alpha, beta, left, right)’. The higher the value of predict_result(i, j), the higher the possibility of the existence of potential association between lncRNA l(i) and disease d(j).

Notice: the path of the association data can be change in BiWalkLDA.m 