load -ascii outdat40.dat
nu=outdat40(:,1);
X2r40=outdat40(:,2);
X2i40=outdat40(:,3);
a2240=outdat40(:,8);
b2240=-outdat40(:,9);
dd40=2;
m2240=dd40;
rao40=(X2r40+complex(0,1)*X2i40)./(dd40-nu*m2240-nu.*a2240+complex(0,1)*nu.*b2240);
% choose D=15 m (and L=30 m), calculate frequency in seconds
om=sqrt(nu)’*sqrt(15/9.81);
% choose mean wave period T2=6 sec of the Jonswap spectrum
% mean wave period T1=1.086 T2
T2=6;
omt2=om*T2;
omt1=omt2*1.086;
% Calculate the Jonswap spectrum
NN=40;
for nn=1:NN;
    sigma(nn)=0.07;
    if omt1(nn) > 5.24;
    sigma(nn)=0.09;
    continue
    end
end
YY=exp(-((0.191*omt1-1)./(sqrt(2)*sigma)).^2);
SS=(155./(omt1.^4.*omt2)).*exp(-944./(omt1.^4)).*(3.3).^YY;
hold on
axis([0 2 0 2.5])
g040=plot(nu,abs(rao40),’k -’)
get(g040)
g5=plot(nu,T2*10*SS,’k -.’);
get(g5);
set(g5,’MarkerSize’,[10])
set(gca,’FontSize’,20)


intresp=0.0;
for nn=1:NN-1;
    aa=0.5*(SS(nn)+SS(nn+1));
    bb=0.5*(abs(rao40(nn))^2+abs(rao40(nn+1))^2);
    dom=om(nn+1)-om(nn);
    intresp=intresp+aa*bb*dom*T2;
end
sqrt(intresp)