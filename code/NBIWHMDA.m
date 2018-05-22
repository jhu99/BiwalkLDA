A=textread('../data/knowndiseasemicrobeinteraction.txt');
% nd:the number of diseases
% nm:the number of microbe
% pp:the number of known diseae-microbe associations
nd=max(A(:,1)); %第一列的最大值
nm=max(A(:,2)); %第二列的
[pp,qq]=size(A);
%interaction: adjacency matrix for the disease-microbe association network
%interaction(i,j)=1 means microbe j is related to disease i
for i=1:pp
        interaction(A(i,1),A(i,2))=1;
end
save ./temp/interaction.txt interaction -ascii;

nFold = 5;
Pint = find(interaction);
Nint = length(Pint);
posFilt = crossvalind('kfold', Nint, nFold);
overauc = 0;
for foldID = 1:nFold
    interaction_new = load('./temp/interaction.txt');
    for i = 1:Nint
        if(posFilt(i) == foldID)
            interaction_new(Pint(i)) = 0;
        end
    end
   %%check 0
   for i=1:nd
       if length(find(interaction_new(i,:)))==0
           i
       end
   end
   sim_dis_path = load("../data/sim_dis_path.txt");
   sim_dis_chemcial = load("../data/sim_dis_chemical.txt");
   for i = 1:nd
       sim_dis_path(i,:) = sim_dis_path(i,:)/sum(sim_dis_path(i,:)); 
   end
   mean(sim_dis_path(7,:))
   mean()
   Zscore = birandom(interaction_new, 0.4, 2, 2);
   score = zeros(nm*nd,3);
   for i = 1:nd
       for j = 1:nm
           score((i-1)*nm + j, 1) = (i-1) * nm +j;
           score((i-1)*nm + j, 2) = Zscore(i,j);
           score((i-1)*nm + j, 3) = 0;
       end
   end
   for i = 1:Nint
       [I, J] = ind2sub(size(interaction), Pint(i));
       index = (I - 1)*nm + J;
       if(posFilt(i) == foldID)
           score(index,3) = 1;
       else
           score(index, 3) = -1;
       end
   end
   
   score = sortrows(score,2);
   [len1, len2] = size(score);
   sum_0 = 0;
   sum_1 = 0;
   for i = 1:len1
       if(score(i, 3) == 0)
           sum_0 = sum_0 + 1;
       end
       if(score(i, 3) == 1)
           sum_1 = sum_1 + 1;
       end
   end
    cnt_0 = 0;
    cnt_1 = 0;
    k = 1;
    for i = len1:-1:1
        if(score(i,3) ==0)
            cnt_0 = cnt_0 + 1;
            fp(k,1) = cnt_0;
            fp(k,2) = cnt_1;
            k = k+1;
        end
        if(score(i,3) ==1)
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
    overauc = overauc + auc;
    disp(foldID)
    disp(auc)
end
    overauc = overauc/nFold