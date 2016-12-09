function [LL, AIC] = covboot(predmat, obsmat, k, wts)
%covboot = @(fE)  mvnpdf(nanmean(fE,2),zeros(43,1),cov(fE'));
%Rows are observations, columns are variables

subs=logical(~isnan(sum(predmat))&~isnan(sum(obsmat)));
diffmat=obsmat(:,subs)-predmat(:,subs);
zeromat = zeros(size(predmat(:,subs),1),1);
wts=wts(subs);
wts=wts/sum(wts);

n=size(diffmat,1);

sig=cov(diffmat');

try
    LL=sum(wts.*mvnpdf_log(diffmat', zeromat', sig));
catch
    LL=NaN;
end

AIC= 2*k - 2*LL + 2*k*(k+1)/(n-k-1);
end