function [Rt] = birandom(interaction, alpha,l,r)
    [nd, nm] = size(interaction);
   clear temp
  
    [kd, km] = Similarity(nd, nm, interaction);
    normWdd = normFun(kd);
    normWrr = normFun(km);
    R0 = interaction/sum(interaction(:));
    Rt = R0;
    for t=1:max(l,r)
      ftl = 0;
      ftr = 0;
      %random walk on the microbe similarity network
      if(t<=l)
         nRtleft = alpha * Rt * normWrr + (1-alpha)*R0;
         ftl = 1;
      end   
      %random walk on the disease similarity network
      if(t<=r)
         nRtright = alpha *  normWdd * Rt + (1-alpha)*R0;
         ftr = 1;
      end
      %Rt: predictive association scores between each microbe-disease pair
      Rt =  (ftl*nRtleft + ftr*nRtright)/(ftl + ftr);
    end
end

