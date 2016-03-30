function C1 = mod_1(beta,x)

hd=exp(beta(:,1));%.*(sqrt(2)./2);
th=2*pi./x;

D= exp(-th.^2*hd^2/2); %Gaussian dispersal kernel
C1(:,1) = (1 - D).^-1;