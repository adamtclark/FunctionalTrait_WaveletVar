function [pd] = trait_kernel(trp0, D, h, a)
%Takes arguments for a quadratic kernel and returns
%probability densities for each species
%Inputs:
%trp0: trait value around which to expand kernelo
%D trait list for all species
%h and a - constants for kernel shape
%Returns:
%pd, a vector of probability densities for each species
  
    G=1;
    pd = (1/G).*(1-((trp0-D)./h).^2).^a;
    pd(trp0==D)=0;
    pd(pd<0) = 0;
    %hlst = h;
    %Log weighting
    %if(nanmax(Navg)>100)
    %    pd=pd.*log10(Navg);
    %else
    %    pd=pd.*log10(Navg.*50e4);
    %end
    %pd(Navg==0)=0;
end