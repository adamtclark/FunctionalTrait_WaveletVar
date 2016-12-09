function C1 = mod_2(beta,x)

hd=exp(beta(1)).*(sqrt(2)./2);
hk=ilogit(beta(2)).*25.*sqrt(6);
P1=ilogit(beta(3)).*5;

th=2*pi./x;
D=(1+(hd*th).^2).^-1;
K= 8*besselj(2,hk*th).*((hk*th).^(-2));
C1(:,1) = (1-D+P1*K).^-1;