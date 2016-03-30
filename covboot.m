function [LL, AIC] = covboot(predmat, obsmat, k)
%covboot = @(fE)  mvnpdf(nanmean(fE,2),zeros(43,1),cov(fE'));
%Rows are observations, columns are variables
%Eest = fE';
subs=logical(~isnan(sum(predmat))&~isnan(sum(obsmat)));
diffmat=obsmat(:,subs)-predmat(:,subs);
zeromat = zeros(size(predmat(:,subs),1),1);

n=size(diffmat,1);

sig=cov(diffmat');

%wts=log(Navg(subs))/sum(nansum(Navg(subs)));
%wts(Navg==0)=0; wts=wts/sum(wts);

if(sum(eig(sig)<0)>0)
    LL=NaN;
else
    LL=mean(mvnpdf_log(diffmat', zeromat', sig));
    %LL=sum(wts.*mvnpdf_log(diffmat', zeromat', sig));
end

%LL=0;
%for(ii=1:size(diffmat,2))
%    sig=2*Aw.*(obsmat(:,ii)*obsmat(:,ii)');
%    LL=LL+mvnpdf_log(predmat(:,ii), predmat(:,ii), sig);
%end
%LL=LL./size(diffmat,2);

AIC= 2*k - 2*LL + 2*k*(k+1)/(n-k-1);
end