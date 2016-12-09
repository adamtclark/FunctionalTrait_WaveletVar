function [modfit]=do_lik_boot_cont_fitting(E, Nsave, h0, Dmat, nbootstraps, Awfit, weighttype, subs2)
%Does a weighted analysis of parameters based on a continous trait
%h0 is the width of the kernel, or the initial guess if useN is true
%useN=1 means that h will be adjusted to keep nspc species in each bin
%nsplits is the number of splits to make along the trait
%Dmat is a list of traits, 1 row per species
%Set nsplits = 0 to estimate wavelet variance for each observed species trait.

wts=getweights(weighttype, Nsave);

h0start=h0;

if(size(E,1)==28)
    load('scale_5m_xmap2.mat')
    scale=scale*5;
    L=[ 101   201];
else
    load('scale.mat')
    L=[ 1000   500];
end

scsubs=(scale<=120);
scale=scale(scsubs);

% clean up
keepsp=(mean(Nsave,2)>10)&(~isnan(Dmat))&(isfinite(Dmat));
if(exist('subs2'))
    keepsp=keepsp&subs2';
end

Nsave=Nsave(keepsp,:);
E=E(scsubs,keepsp,:);
Dmat=Dmat(keepsp);
wts=wts(keepsp);
wts=wts/sum(wts);

EwtsvPERM=NaN(size(E,1), length((Dmat)), nbootstraps);
EwtsvBOOT=NaN(size(E,1), length((Dmat)), nbootstraps);

diffmatPERM=NaN(size(E,1), nbootstraps);
diffmatBOOT=NaN(size(E,1), nbootstraps);

obsmatPERM=NaN(size(E,1), size(E,2), nbootstraps);
predmatPERM=NaN(size(E,1), size(E,2), nbootstraps);
trtmatPERM=NaN(size(E,2), nbootstraps);

obsmatBOOT=NaN(size(E,1), size(E,2), nbootstraps);
predmatBOOT=NaN(size(E,1), size(E,2), nbootstraps);
trtmatBOOT=NaN(size(E,2), nbootstraps);

for(permtype=1:2)
    %Run permutation test for permtype=1; run bootstrapping for permtype=2
for(nbootst=1:nbootstraps)
    if(permtype==1)
        bootsubs=1:size(Dmat,1);
        subs=randsample(bootsubs, size(Dmat,1));
    else
        if(nbootst==1)
            subs=1:size(Dmat,1);
        else
            subs=sort(randsample(1:size(Dmat,1), size(Dmat,1), 'true'));
        end
    end
    
    R=nanmean(Dmat(subs),2);
    R=R(~isnan(R));
    R(R==0)=[];

    trp = sort(unique((Dmat(subs))));

    if(h0start==0)
       h0= 1.06*std(R)*length(R)^(-1/5)*2.2; %Rule of thumb for smoother kernel
    end

    if(permtype==1)
        Em = nanmean(E(:,bootsubs), 3);
        Em(isnan(Em(:)))=0;
    else
        Em = nanmean(E(:,subs), 3);
        Em(isnan(Em(:)))=0;
    end
        
    pdlst = zeros(length(trp), length(subs));
    for i=1:length(trp)
        pdlst(i,:)=trait_kernel(trp(i), Dmat(subs), h0, 1, wts(subs));
        pdlst(i,:)=pdlst(i,:)/nansum(pdlst(i,:));
    end
    pdlst(isnan(pdlst))=0;

    Ewt= Em*pdlst';

    predmat = NaN(length(scale), length(subs));
    Ewt(:,sum(Ewt)==0)=repmat(nanmean(Ewt(:,sum(Ewt)~=0),2), 1, sum(sum(Ewt)==0));
    for(i=1:length(trp))
        matchlst=find(Dmat(subs)==trp(i));
        predmat(:,matchlst)=repmat(Ewt(:,i),1,length(matchlst));
    end
    warning('off', 'all')
    dfUsedtmp = trace(predmat/Em);
    warning('on', 'all')

    [~,ord]=sort(Dmat(subs));
    
    %%%%% RUN FITTING ALGORITHM
    if(permtype==2 & nbootst==1)
        [totfit]=getcontfit(predmat, eye(size(predmat,2)), scale, 0, Awfit);
    end
    %%%%%
    
    if(permtype==1)
        [logLlstPERM(nbootst), AICPERM(nbootst)]=covboot(predmat, Em, dfUsedtmp, wts);
        
        predmatsv=predmat(:,ord);
        EwtsvPERM(:,:,nbootst)=predmatsv;
        
        diffmatPERM(:,nbootst)=mean((Em-predmat).^2');
        obsmatPERM(:,:,nbootst)=Em(:,ord);
        predmatPERM(:,:,nbootst)=predmat(:,ord);
        trtmatPERM(:,nbootst)=Dmat(subs(ord));
    else
        [logLlstBOOT(nbootst), AICBOOT(nbootst)]=covboot(predmat, Em, dfUsedtmp, wts);
        
        predmatsv=predmat(:,ord);
        EwtsvBOOT(:,:,nbootst)=predmatsv;
        
        diffmatBOOT(:,nbootst)=mean((Em-predmat).^2');
        obsmatBOOT(:,:,nbootst)=Em(:,ord);
        predmatBOOT(:,:,nbootst)=predmat(:,ord);
        trtmatBOOT(:,nbootst)=Dmat(subs(ord));
    end
end
end

modfit = struct('logLlstPERM', logLlstPERM, 'aicPERM', AICPERM, 'logLlstBOOT', logLlstBOOT, 'aicBOOT', AICBOOT, ...
    'diffmatPERM', diffmatPERM, 'obsmatPERM', obsmatPERM, 'predmatPERM', predmatPERM, 'trtmatPERM', trtmatPERM, ...
    'diffmatBOOT', diffmatBOOT, 'obsmatBOOT', obsmatBOOT, 'predmatBOOT', predmatBOOT, 'trtmatBOOT', trtmatBOOT, ...
    'h0', h0, 'totfit', totfit, 'Dmatnew', Dmat);
end




