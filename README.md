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
Description This file is README file of the package of NetCoffee2, developed by Yiqun Gao.
###########################################################################################
This section introduces functions of different files. 

Files:

1.Dataset

disease_Name.txt：The name of all the diseases

lncRNA_Name.txt：The name of all the lncRNA-disease

interMatrix.mat：Known lncRNA-disease associations

disSim_Jaccard.mat：Disease similairty calculate by gene ontology annotations information

2.Code

BiWalkLDA.m: BiwalkLDA framework to predict potential lncRNA-disease association

birandom.m: bi-random walk algorithm

Similarity.m: Gaussian interaction profile kernel similarity for disease and lncRNA

normFun.m: Laplacian normalization
###########################################################################################