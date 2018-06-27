% BiWalkLDA: Prediction of lncRNA-disease associations based on bi-random
% walk
%  
%%%%%%%%%%
%   interMatrix.mat: an n*m association matrix between lncRNAs and diseases, n is 
%the number of lncRNAs, and m is the number of diseases
%   disSim_Jaccard.mat: an m*m similarity matrix of disease

function [Biwalk_recover] = BiWalkLDA(alpha, beta, left, right)
%% load data
    LD=importdata('../Datasets/Dataset1/interMatrix.mat');    
    dissim=importdata('../Datasets/Dataset1/disSim_Jaccard.mat');
    interaction = LD;
    
%% Calculate lncRNA similarity and disease similarity
    [nl,nd]=size(LD);
    [kl, kd] = Similarity(nl, nd, interaction);
    kd = dissim * alpha +(1-alpha)*kd;
    lncSim = kd;
    
%% complete interaction information for a new lncRNA
    for i=1:nl
        if length(find(interaction(i,:)))==0
            rowVec=lncSim(i,:); 
            rowVec(i)=0;
            simNeighbors=find(rowVec>=mean(mean(lncSim))); 
            if length(simNeighbors)
                new_row=zeros(1,nd);
                for l=1:length(simNeighbors)
                   new_row=new_row+interaction(simNeighbors(l),:);       
                end
                new_row=new_row/length(simNeighbors);      
                interaction(i,:)=new_row;     
            end
        end
    end
    %% extract feature vectors of lncRNAs and diseases
    Biwalk_recover = birandom(interaction,kd, kl, beta, left, right); 
end
   