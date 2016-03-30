function [modfit]=do_lik_boot_cont_fitting(E, h0, Dmat, nbootstraps, plotxname, unlog, Awfit)
%Does a weighted analysis of parameters based on a continous trait
%h0 is the width of the kernel, or the initial guess if useN is true
%useN=1 means that h will be adjusted to keep nspc species in each bin
%nsplits is the number of splits to make along the trait
%Dmat is a list of traits, 1 row per species
%Set nsplits = 0 to estimate wavelet variance for each observed species trait.

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
keepsp=(~isnan(Dmat))&(isfinite(Dmat));%&(E(end,:)<100)';

E=E(scsubs,keepsp,:);
Dmat=Dmat(keepsp);

%logLlstPERM=zeros(nbootstraps, 1);
%AICPERM=zeros(nbootstraps, 1);

%logLlstBOOT=zeros(nbootstraps, 1);
%AICBOOT=zeros(nbootstraps, 1);

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
        %if(nbootst==1)
            bootsubs=1:size(Dmat,1);
        %else
        %    bootsubs=sort(randsample(1:size(Dmat,1), size(Dmat,1), 'true'));
        %end
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
        pdlst(i,:)=trait_kernel(trp(i), Dmat(subs), h0, 1);
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
        [logLlstPERM(nbootst), AICPERM(nbootst)]=covboot(predmat, Em, dfUsedtmp);
        
        predmatsv=predmat(:,ord);
        EwtsvPERM(:,:,nbootst)=predmatsv;
        
        diffmatPERM(:,nbootst)=mean((Em-predmat).^2');
        obsmatPERM(:,:,nbootst)=Em(:,ord);
        predmatPERM(:,:,nbootst)=predmat(:,ord);
        trtmatPERM(:,nbootst)=Dmat(subs(ord));
    else
        [logLlstBOOT(nbootst), AICBOOT(nbootst)]=covboot(predmat, Em, dfUsedtmp);
        
        predmatsv=predmat(:,ord);
        EwtsvBOOT(:,:,nbootst)=predmatsv;
        
        diffmatBOOT(:,nbootst)=mean((Em-predmat).^2');
        obsmatBOOT(:,:,nbootst)=Em(:,ord);
        predmatBOOT(:,:,nbootst)=predmat(:,ord);
        trtmatBOOT(:,nbootst)=Dmat(subs(ord));
    end
end
end

%logLlstPERM=zeros(size(E,2), nbootstraps);
%logLlstBOOT=zeros(size(E,2), 1);
%zeromat=zeros(size(E,1), nbootstraps);
%zeromatsm=zeros(size(E,1), 1);

%for(i=1:size(E,2))
%    deltamat=predmatPERM(:,i,:)-obsmatPERM(:,i,:);
%    deltamat=deltamat(:,:);
%    
%    deltamatsm=predmatBOOT(:,i,1)-obsmatBOOT(:,i,1);
%    
%    cvmat=cov(deltamat');
%
%    logLlstPERM(i,:)=mvnpdf_log(deltamat', zeromat', cvmat);
%    logLlstBOOT(i)=mvnpdf_log(deltamatsm, zeromatsm, cvmat);
%    
%    pr=predmatPERM(:,i,:);
%    ob=obsmatPERM(:,i,:);
%    dfUsedtmp = trace(pr(:,:)/ob(:,:));
%end

if(nbootstraps>0) %Plot output
    trp=sort((Dmat));
    if(unlog)
        trp=10.^trp;
    end
    Ewt=repmat(quantile(nanmean(EwtsvBOOT,3),0.5,2), 1, length(trp));
    
    CR1=zeros(size(Ewt));
    CR1(quantile(EwtsvBOOT,0.025,3)>Ewt)=1;
    CR2=zeros(size(Ewt));
    CR2(quantile(EwtsvBOOT,0.975,3)<Ewt)=1;

    Ewt_int=nanmean(EwtsvBOOT,3);

    
    [~,unq]=unique(trp);
    %pcolor(scale,trp(unq), (Ewt_int(:,unq)./Ewt(:,unq))')
    pcolor(scale,trp(unq), ((Ewt_int(:,unq)./Ewt(:,unq))'-1)*100)
    
    %set(gcf, 'colormap', gray) 
    set(gca, 'xscale', 'log')
    if(unlog)
        set(gca, 'yscale', 'log')
    end
    
    shading flat
    hold on
    
    %CR1smooth=crsmooth(CR1);
    %CR2smooth=crsmooth(CR2);
    CR1smooth=medfilt2(CR1(:,unq),[5 5]);
    CR2smooth=medfilt2(CR2(:,unq),[5 5]);
    
    if(sum(CR1smooth(:))>10)
        [c1 h1]=contourf(scale,trp(unq),CR1smooth', 'LevelList', 1);
    end
    if(sum(CR2smooth(:))>10)
        [c2 h2]=contourf(scale,trp(unq),CR2smooth', 'LevelList', 1);
    end
    hold off
    
    if(sum(CR1smooth(:))>10)
        hpatch1=findobj(h1, 'Type', 'patch');
        hh1=hatchfill(hpatch1, 'cross', 45, 7);
        set(hh1, 'Color', 'black', 'LineWidth', 1)
        set([h1], 'LineStyle', 'none')
    end

    if(sum(CR2smooth(:))>10)
        hpatch2=findobj(h2, 'Type', 'patch');
        hh2=hatchfill(hpatch2, 'cross', 45, 7);
        set(hh2, 'Color', 'white', 'LineWidth', 1)
        set([h2], 'LineStyle', 'none')
    end
    
    zmin = max([quantile((Ewt_int(:)./Ewt(:)-1)*100, 0.025), -2]);
    zmax = min([quantile((Ewt_int(:)./Ewt(:)-1)*100, 0.975), 2]);
    dz = mean([1-zmin, zmax-1]);
    caxis([1-dz 1+dz]) %Set reasonable z axes limits
    colorbar
    ylabel(plotxname,'fontsize',15)
    xlabel('scale (m)','fontsize',15)
    xlim([2, max(scale)])
    
    
    NumTicks = 5;
    L = get(gca,'XLim');
    tcks=unique(round(exp(log(L(1)):(log(L(2))-log(L(1)))/NumTicks:log(L(2)))./5).*5);
    set(gca,'XTick',tcks)

    if(unlog)
        L = get(gca,'YLim');
        tcks=unique(round(exp(log(L(1)):(log(L(2))-log(L(1)))/NumTicks:log(L(2)))./5).*5);
        set(gca,'YTick',tcks)
    else
        L = get(gca,'YLim');
        tcks=unique(round(((L(1)):((L(2))-(L(1)))/NumTicks:(L(2)))./5).*5);
        set(gca,'YTick',tcks)
    end
end

modfit = struct('logLlstPERM', logLlstPERM, 'aicPERM', AICPERM, 'logLlstBOOT', logLlstBOOT, 'aicBOOT', AICBOOT, ...
    'diffmatPERM', diffmatPERM, 'obsmatPERM', obsmatPERM, 'predmatPERM', predmatPERM, 'trtmatPERM', trtmatPERM, ...
    'diffmatBOOT', diffmatBOOT, 'obsmatBOOT', obsmatBOOT, 'predmatBOOT', predmatBOOT, 'trtmatBOOT', trtmatBOOT, ...
    'h0', h0, 'totfit', totfit, 'Dmatnew', Dmat);
end




