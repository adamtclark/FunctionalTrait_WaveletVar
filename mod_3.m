function C1 = mod_3(beta,x)

hd=exp(beta(1)).*(sqrt(2)./2);
hk=ilogit(beta(2)).*25.*sqrt(6);
hh=ilogit(beta(3)).*25.*sqrt(6);
P1=ilogit(beta(4)).*5;
P2=beta(5);


th=2*pi./x(:,1);
Cij=x(:,2);

D=(1+(hd*th).^2).^-(1);
K=8*besselj(2,hk*th).*((hk*th).^(-2));
H=exp(-((hh*th).^2)/2);

C1(:,1) = (1-P2*H.*Cij)./(1-D + P1*K);