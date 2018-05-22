% SMFLDA: Prediction of lncRNA-disease associations based on inductive matrix completion
%  
%%%%%%%%%%
%   interMatrix.mat: an n*m association matrix between lncRNAs and diseases, n is 
%the number of lncRNAs, and m is the number of diseases
%   lncSim.mat: an n*n sequence similarity matrix of lncRNAs
%   disSim_Jaccard.mat: an m*m similarity matrix of disease

%% configuration
addpath('SIMC');

%% load data
LD=importdata('../Datasets/Dataset1/interMatrix.mat');
lncSim=importdata('../Datasets/Dataset1/lncSim.mat');    
dissim=importdata('../Datasets/Dataset1/disSim_Jaccard.mat');

Pint = find(LD); % pair of interaction
Nint = length(Pint);
%posFilt = crossvalind('Kfold', Nint, nFold);
over_auc = 0;
for foldID = 1 : Nint
    interaction = LD;
    %for i = 1:Nint
        %if(posFilt(i)== foldID)
            %interaction(Pint(i)) = 0;
        %end
    %end
    interaction(Pint(foldID)) = 0;
    %% computing Gaussian interaction profile kernel of lncRNAs
    [nl,nd]=size(interaction);
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
    
    [LL,~]=gKernel(nl,nd,interaction); 
    %% extract feature vectors of lncRNAs and diseases
    lnc_feature=pca_energy(LL,0.8);
    dis_feature=pca_energy(dissim,0.6);
    %% using inductive matrix completion to complete the association matrix of lncRNA-disease
    Omega=find(interaction==1);    
    M_recover=SIMC(interaction,Omega,lnc_feature,dis_feature); 
    %% calculate AUC
    [row_m, col_m] = size(M_recover);
    clear result
    for i = 1:row_m
        for j = 1:col_m
            index = (i-1)* col_m + j;
            result(index,1) = index;
            result(index,2) =  M_recover(i,j);
            result(index,3) =  0;
        end
    end
    for i = 1:Nint
            [I, J] = ind2sub(size(interaction), Pint(i));
            index = (I-1)*col_m + J;
            if(i == foldID)
                result(index,3) =  1;
            else
                result(index,3) =  -1;
            end
    end
    [len1, ~] = size(result);
    sum_0 = 0;
    sum_1 = 0;
    for i = 1:len1
        if(result(i,3)==0)
            sum_0 = sum_0 + 1;
        end
        if(result(i,3)==1)
            sum_1 = sum_1 + 1;
        end
    end
    score_final = sortrows(result, 2);
    [len1, len2] = size(score_final);
    cnt_0 = 0;
    cnt_1 = 0;
    k = 1;
    for i = len1:-1:1
        if(score_final(i,3) ==0)
            cnt_0 = cnt_0 + 1;
            fp(k,1) = cnt_0;
            fp(k,2) = cnt_1;
            k = k+1;
        end
        if(score_final(i,3) ==1)
            cnt_1 = cnt_1 + 1;
            fp(k,1) = cnt_0;
            fp(k,2) = cnt_1;
            k = k+1;
        end
    end
    for i = 1:k-1
        if(fp(i,1)~=0)
            fp(i, 1) = fp(i, 1)/sum_0;     
        end
        if(fp(i,2)~=0)
            fp(i, 2) = fp(i, 2)/sum_1;     
        end      
    end
    plot(fp(:,1),fp(:,2))
    clear area;
    auc = 0;
    for i = 2:k-1
        auc = auc + (fp(i,1) - fp(i-1,1))*fp(i,2);
    end
    cv
    auc
    over_auc = over_auc + auc;
end
over_auc = over_auc/Nint;

