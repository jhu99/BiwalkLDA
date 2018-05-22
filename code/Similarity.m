function [result_dis,result_mic]=gKernel(nd,nm,inter_lncdis)
%gKernel compute Gaussian interaction profile kernel
%   Usage:  [result_lnc,result_dis]=gKernel(nl,nd,inter_lncdis)
%	Inputs:
%			nl: the number of lncRNAs
%			nd:	the number of diseases
%			inter_lncdis: an nl*nd association matrix between lncRNAs and diseases
%
%	Outputs:
%			result_lnc: Gaussian interaction profile kernel of lncRNAs
%			result_dis: Gaussian interaction profile kernel of diseases

    for i=1:nd
        sl(i)=norm(inter_lncdis(i,:))^2;
    end
    gamal=nd/sum(sl')*1;
    for i=1:nd
        for j=1:nd
            pkl(i,j)=exp(-gamal*(norm(inter_lncdis(i,:)-inter_lncdis(j,:)))^2);
        end
    end        
    for i=1:nm
        sd(i)=norm(inter_lncdis(:,i))^2;
    end
    gamad=nm/sum(sd')*1; 
    for i=1:nm
        for j=1:nm
            pkd(i,j)=exp(-gamad*(norm(inter_lncdis(:,i)-inter_lncdis(:,j)))^2);
        end
    end 
 
    result_dis=pkl;
    result_mic=pkd;
end
   

