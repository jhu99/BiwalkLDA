The code contains three .m files:

    1. Similarity: construct the adjacency matrix for the microbe-disease association network and calculate microbe similarity and disease similarity

    2. normFun: Laplacian normalization for microbe similarity and disease similarity.

    3. BiRWHMDA: bi-random walk on the heterogeneous network to calculate association scores for each microbe-disease pair. 



Instructions:

    step 1: run Similarity.m to construct the adjacency matrix for the microbe-disease association network and calculate microbe similarity and disease similarity
              
               input: knowndiseasemicrobeinteraction.txt

               output: interaction.mat  microbesimilarity.mat  diseasesimilarity.mat

    step 2: run BiRWHMDA.m to calculate association scores for each microbe-disease pair 

               input: interaction.mat  microbesimilarity.mat  diseasesimilarity.mat

               output: Result.mat


All files of data and code should be stored in the same folder.