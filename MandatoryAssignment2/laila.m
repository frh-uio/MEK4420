clear;
NN=30; 
in=linspace(1,NN,NN);
nu=linspace(0.01,2,NN); %nu=1.2; %K

% characteristics of geometry 
load -ascii box1.dat; %D=1,L=2
%load -ascii box2.dat xm=box(:,1); %D=1,L=1
%load -ascii box3.dat xm=box(:,1); %D=1;L=0.1
xm=box1(:,1); 
ym=box1(:,2);
xp=box1(:,3);
yp=box1(:,4);


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
sff22=0;
for k=1:NN
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
            zz=nu(k)*(yb-complex(0,1)*xa); 
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
            help1=(n1(j)* (imag(f1)+complex(0,1)*imag(f2)) +n2(j)*(real(f1)+complex(0,1)*real(f2)) )*nu(k)*ds(j);
            ss(i,j)=(arg0+arg1+help1); 
        end
    end
    rhs=gg*n2; 
    phi2=ss\rhs;
    ff22=phi2.*n2.*ds;
    sff22 = [sff22,sum(ff22)];
    
end


a22=real(sff22);
a22=a22(2:end);
b22=-imag(sff22);
b22=b22(2:end);

% %print pretty table:
% tab=[0,0.5,1,1.5,2];
% for g=0:lenght(tab)
%     fprintf('%8.2f %8.3f\n', [a22(g),b22(g)]')
% end
% 

plot(nu,a22, 'k -');
hold on
plot(nu,b22, 'k .-');
legend('a_{22}/\rho D^2','b_{22}/\rho \omega D^2');
set(gca,'FontSize',20)

% %added mass and damping
% ff22=phi2.*n2.*ds;
% %KD = D*w^2/g = 1.2 -> w = sqrt(11.7722), (g=9.81)
% omega = sqrt(11.7722);
% rho=1;
% D=1;
% sff22=sum(ff22);
% a22=rho*D^2*real(sff22);
% b22=-omega*rho*D^2*imag(sff22)
% 
% %damping from the energy balance
% phi0=exp(nu*(by-complex(0,1)*bx)); 
% AM2=complex(0,1)*(phi2.*(nu*n2-nu*complex(0,1)*n1)-n2).*phi0.*ds; 
% AP2=complex(0,1)*(phi2.*(nu*n2+nu*complex(0,1)*n1)-n2).*conj(phi0).*ds; 
% sAM2=sum(AM2); 
% sAP2=sum(AP2);
% b22e = 0.5*omega*( (abs(sAM2))^2 + (abs(sAP2))^2 )
% 
% %exiting force
% g=9.81;
% L=2;
% K=nu;
% XFK_rho1= g*L*exp(-K*D) * (sin(K*L/2.0))/(K*L/2);
% %XFK_rho2= g*L*exp(-K*D); %narrow box section
% 
% %damping from the Haskins relations
% b22H = omega*rho*((abs(XFK_rho1))^2/(rho*g)^2)