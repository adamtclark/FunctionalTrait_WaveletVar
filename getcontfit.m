function[totmat, b1mat, b2mat, AICmat] = getcatfit(E, dispbin, scale, doplot, Awfit)
    if(~exist('doplot'))
       doplot=1; 
    end
    
    if(~exist('doempdist'))
        doempdist=0;
    end
    
    %Fit mean parameters
    scalefit=scale(1:40);
    if(~exist('Awfit'))
        [~,~,~,Awfit]=spec_analitycal(inline('1+K*0','K','hd'),1000,500,'dx',1,'param',[8 1],...
            'scale',scalefit);
    end
    
    colvar=size(dispbin,2);
    collst=0:(1/colvar):1;
    collst=collst(2:end);
    revcollst=collst(end:-1:1);

    clear b1mat b2mat AIC1mat AIC2mat
    for(ii = 1:size(dispbin,2))
        Ematfit=E(:,logical(dispbin(:,ii)));
        Ematfit=Ematfit(1:40,~isnan(sum(Ematfit)));
        
        y=nanmean(Ematfit, 2);
        SIGMA=2*Awfit.*(y*y');

        beta0=log(4);
        LO=-1;UP=5;

        [b1,~,~,CovB] = gnlinfit(scalefit,y,SIGMA,@mod_1,beta0,LO,UP);
        b1_sd=sqrt(diag(CovB));
        
        p1=mod_1(b1, scalefit);
        
        beta0=[b1 0 -2.3];
        
        LO=[0 -8.2 -6.6];
        UP=[5 1.4078 -1.4];
        
        [b2,~,~,CovB] = gnlinfit(scalefit,y,SIGMA,@mod_2,beta0,LO,UP);
        b2_sd=sqrt(diag(CovB));
        p2=mod_2(b2, scalefit);
        
        if(doplot==1)
            hold all
            loglog(scalefit, p2, 'linestyle', '-', 'color', [collst(ii) 0 revcollst(ii)], 'linewidth', 2);
        end
        

        Lb1=mvnpdf_log(y, p1, SIGMA);
        Lb2=mvnpdf_log(y, p2, SIGMA);
        k1=1; k2=3;

        AIC1=2*k1-2*(Lb1)+2*k1*(k1+1)/(size(scalefit,1)-k1-1);
        AIC2=2*k2-2*(Lb2)+2*k2*(k2+1)/(size(scalefit,1)-k2-1);
        
        b1mat(ii,:,1)=b1;
        b2mat(ii,:,1)=b2;
        
        b1mat(ii,:,2)=b1_sd;
        b2mat(ii,:,2)=b2_sd;
        
        AICmat(ii,1)=AIC1;
        AICmat(ii,2)=AIC2;
    end
    totmat=[[b1mat(:,:,1), zeros(size(b1mat,1), 2)]', ...
            [b1mat(:,:,2), zeros(size(b1mat,1), 2)]', ...
            b2mat(:,:,1)', ...
            b2mat(:,:,2)', ...
            [AICmat'; zeros(1,size(b1mat,1))]];
end