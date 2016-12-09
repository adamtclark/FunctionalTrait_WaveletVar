function [wts] = getweights(weighttype, Nsave)
    if(strcmp(weighttype, 'none'))
        wts=ones(size(Nsave));
    elseif(strcmp(weighttype, 'lin'))
        wts=Nsave;
    else
        wts=log(Nsave+1);
    end
    wts=wts/sum(wts);
end