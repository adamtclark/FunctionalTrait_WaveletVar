function [dLL, LLcmp, AICcmp, meanmat, meanmatperm] = do_lik_boot_cat(E, dispbin, dispnamesbin, scale, Nsave, niter, weighttype, subs2)
wts=getweights(weighttype, Nsave);
%How we implement weights:
%1) get weighted mean for Etmp
%2) get weighted mean for LL (since we are down-weighting things in the
%mean, we should also down-weight them in the likelihood test)

subs = ~isnan(mean(E,1))&(nansum(dispbin,2)>0)'&(Nsave>10)';
if(exist('subs2'))
    subs=subs&subs2;
end

scsubs=(scale<=120);
scale=scale(scsubs);
E = E(scsubs,subs);
dispbin=dispbin(subs,:);
wts=wts(subs);
wts=wts/sum(wts);

if(~exist('niter'))
    niter=20000;
end

meanmatperm=zeros(size(E,1),niter,size(dispbin,2));
meanmat=zeros(size(E,1),niter,size(dispbin,2));

%categories groupe together
clear LL meanmatper;
for(i=1:niter)
    bootsub=1:size(E,2);
    
    disbinperm = dispbin(randsample(bootsub, size(dispbin,1), 'false'),:);
    
    deltaE = zeros(size(E,1), size(E,2));
    n=1;
    for(ibin=1:size(dispbin,2))
        Etmp=E(:,disbinperm(:,ibin)==1);
        
        wts_ii=wts(disbinperm(:,ibin)==1);
        wts_ii=wts_ii/sum(wts_ii);
        
        meanmatperm(:,i,ibin) = nansum(repmat(wts_ii', size(Etmp,1), 1).*Etmp,2);
        
        deltaE(:,n:(size(Etmp,2)+n-1)) = Etmp - repmat(meanmatperm(:,i,ibin), 1, size(Etmp,2));
        
        n=size(Etmp,2)+n;
    end
    
    covsub=cov(deltaE');
    try
        LL(i)=sum(mvnpdf_log(deltaE', zeros(1,size(deltaE,1)), covsub).*wts);
    catch
        LL(i)=NaN;
    end
end

%separate
clear LLS meanmat;
for(i=1:niter)
    %Note - first element is non-permuted
    if(i==1)
        bootsub=1:size(E,2);
    else
        bootsub=randsample(1:size(E,2), size(E,2), 'true');
    end
        
    deltaE = zeros(size(E,1), size(E,2));
    n=1;
    for(ibin=1:size(dispbin,2))
        bootsub2=bootsub(ismember(bootsub, find(dispbin(:,ibin))));
        
        Etmp=E(:,bootsub2);
        
        wts_ii=wts(bootsub2);
        wts_ii=wts_ii/sum(wts_ii);
        
        meanmat(:,i,ibin) = nansum(repmat(wts_ii', size(Etmp,1), 1).*Etmp,2);
        
        deltaE(:,n:(size(Etmp,2)+n-1)) = Etmp - repmat(meanmat(:,i,ibin), 1, size(Etmp,2));
        
       	n=size(Etmp,2)+n;
    end
    
    covsub=cov(deltaE');
    
    try
        LLS(i)=sum(mvnpdf_log(deltaE', zeros(1,size(deltaE,1)), covsub).*wts);%mean(mvnpdf_log(deltaE', zeros(1,size(deltaE,1)), covsub));
    catch
        LLS(i)=NaN;
    end
end

n=size(E,1);
k=1;
k=size(dispbin,2);
AICT=2*k - 2*LL + 2*k*(k+1)/(n-k-1);
AICS=2*k - 2*LLS + 2*k*(k+1)/(n-k-1);

plotcat(dispbin, meanmat, scale, dispnamesbin)
LLcmp = [LL', LLS'];
AICcmp = [AICT', AICS'];

dLL=LLS(1)-LL;
end