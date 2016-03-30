function [dLL, LLcmp, AICcmp, meanmat, meanmatperm] = do_lik_boot_cat(E, dispbin, dispnamesbin, scale, niter)
if(~exist('niter'))
    niter=20000;
end

meanmatperm=zeros(size(E,1),niter,size(dispbin,2));
meanmat=zeros(size(E,1),niter,size(dispbin,2));
%total

clear LL meanmatper;
%obsPERM=zeros(size(E,1), size(E,2), niter);
%predPERM=zeros(size(E,1), size(E,2), niter);
for(i=1:niter)
    %if(i==1)
        bootsub=1:size(E,2);
    %else
    %    bootsub=randsample(1:size(E,2), size(E,2), 'true');
    %end
    
    disbinperm = dispbin(randsample(bootsub, size(dispbin,1), 'false'),:);
    
    deltaE = zeros(size(E,1), size(E,2));
    n=1;
    for(ibin=1:size(dispbin,2))
        Etmp=E(:,disbinperm(:,ibin)==1);        
        
        meanmatperm(:,i,ibin) = mean(Etmp');
        
        deltaE(:,n:(size(Etmp,2)+n-1)) = Etmp - repmat(meanmatperm(:,i,ibin), 1, size(Etmp,2));
        
        %obsPERM(:,n:(size(Etmp,2)+n-1),i)=Etmp;
        %predPERM(:,n:(size(Etmp,2)+n-1),i)=repmat(meanmatperm(:,i,ibin), 1, size(Etmp,2));
        
        n=size(Etmp,2)+n;
    end
    
    covsub=cov(deltaE');    
    LL(i)=mean(mvnpdf_log(deltaE', zeros(1,size(deltaE,1)), covsub));
end

%separate
clear LLS meanmat;
%obsBOOT=zeros(size(E,1), size(E,2), niter);
%predBOOT=zeros(size(E,1), size(E,2), niter);
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
        
        meanmat(:,i,ibin) = mean(Etmp');
        
        deltaE(:,n:(size(Etmp,2)+n-1)) = Etmp - repmat(meanmat(:,i,ibin), 1, size(Etmp,2));
        
        %obsBOOT(:,n:(size(Etmp,2)+n-1),i)=Etmp;
        %predBOOT(:,n:(size(Etmp,2)+n-1),i)=repmat(meanmatperm(:,i,ibin), 1, size(Etmp,2));
        
       	n=size(Etmp,2)+n;
    end
    
    covsub=cov(deltaE');    
    LLS(i)=mean(mvnpdf_log(deltaE', zeros(1,size(deltaE,1)), covsub));
end

%logLlstPERM=zeros(size(E,2), niter);
%logLlstBOOT=zeros(size(E,2), 1);
%zeromat=zeros(size(E,1), niter);
%zeromatsm=zeros(size(E,1), 1);
%for(i=1:size(E,2))
%    deltamat=predPERM(:,i,:)-obsPERM(:,i,:);
%    deltamat=deltamat(:,:);
%    
%    deltamatsm=predBOOT(:,i,1)-obsBOOT(:,i,1);
%    cvmat=cov(deltamat');
%
%    logLlstPERM(i,:)=mvnpdf_log(deltamat', zeromat', cvmat);
%    logLlstBOOT(i)=mvnpdf_log(deltamatsm, zeromatsm, cvmat); 
%end

%LL=mean(logLlstPERM);
%LLS=mean(logLlstBOOT);


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