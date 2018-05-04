nu=1.2; 
% characteristics of geometry 
load -ascii box2.dat xm=box2(:,1);
%load -ascii box1.dat xm=box(:,1); %D=1,L=2
%load -ascii box2.dat xm=box(:,1); %D=1,L=1
%load -ascii box3.dat xm=box(:,1); %D=1;L=0.1
xm=box2(:,1); 
ym=box2(:,2);
xp=box2(:,3);
yp=box2(:,4);

NN=30; 
in=linspace(1,NN,NN);

dx=xp-xm; %step in x
dy=yp-ym; %step in y
ds=((dx).^2+(dy).^2).^(1/2); %step in both dimentions
bx=0.5*(xm+xp); %midpoint
by=0.5*(ym+yp); 
n1=-(yp-ym)./ds; %unity vektor
n2=(xp-xm)./ds; 
% points for Gauss integration on each segment 
xg1=-0.5*dx/sqrt(3)+bx; 
xg2=0.5*dx/sqrt(3)+bx; 
yg1=-0.5*dy/sqrt(3)+by; 
yg2=0.5*dy/sqrt(3)+by; 
% contributions to integral equation, rhs stores rhs, lhs stores lhs 
for i=1:NN 
    for j=1:NN 
% rhs, log(r) term with 2pts Gauss quadrature 
        xa1=xg1(j)-bx(i); 
        xa2=xg2(j)-bx(i); 
        ya1=yg1(j)-by(i); 
        ya2=yg2(j)-by(i); 
        ra1=sqrt(xa1*xa1+ya1*ya1); 
        ra2=sqrt(xa2*xa2+ya2*ya2); 
        g0=(log(ra1)+log(ra2))*0.5;
        % all other terms with midpoint rule 
        xa=bx(j)-bx(i); 
        yb=by(j)+by(i); 
        rb=sqrt(xa*xa+yb*yb); 
        g1=-log(rb); 
        zz=nu*(yb-complex(0,1)*xa); 
        f1=-2*exp(zz)*(expint(zz)+log(zz)-log(-zz)); 
        f2=2*pi*exp(zz); 
        g2=real(f1)+complex(0,1)*real(f2); 
        gg(i,j)=(g0+g1+g2)*ds(j);
    % lhs different for phi2
        arg0=imag(log((xm(j)-bx(i)+complex(0,1)*(ym(j)-by(i)))/ (xp(j)-bx(i)+complex(0,1)*(yp(j)-by(i))))); 
        if j-i == 0 
            arg0=-pi; 
        end
        arg1=imag(log((xm(j)-bx(i)+complex(0,1)*(ym(j)+by(i)))/ (xp(j)-bx(i)+complex(0,1)*(yp(j)+by(i))))); 
        help1=(n1(j)* (imag(f1)+complex(0,1)*imag(f2)) +n2(j)*(real(f1)+complex(0,1)*real(f2)) )*nu*ds(j);
        ss(i,j)=(arg0+arg1+help1); 
    end
end
rhs=gg*n2; 
phi2=ss\rhs;
hold on 
plot(in,real(phi2), 'k -') 
plot(in,imag(phi2),'k +');  
legend('Re(\phi_2)','Im(\phi_2)')
set(gca,'FontSize',20)



