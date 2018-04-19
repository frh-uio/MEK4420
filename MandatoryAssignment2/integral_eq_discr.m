clear;
clc;

nu=1.2;

% characteristics of geometry
load -ascii box1.dat
% load -ascii box2.dat
% load -ascii box3.dat

xm=box1(:,1);
ym=box1(:,2);
xp=box1(:,3);
yp=box1(:,4);

% xm=box2(:,1);
% ym=box2(:,2);
% xp=box2(:,3);
% yp=box2(:,4);

% xm=box3(:,1);
% ym=box3(:,2);
% xp=box3(:,3);
% yp=box3(:,4);

N=30;
in=linspace(1,N,N);

dx=xp-xm;
dy=yp-ym;
ds=((dx).^2+(dy).^2).^(1/2);

% Midpoints
xbar=0.5*(xm+xp);
ybar=0.5*(ym+yp);

n1=-(dy)./ds;
n2=(dx)./ds;

% points for Gauss integration on each segment
xg1=-0.5*dx/sqrt(3)+xbar;
xg2=0.5*dx/sqrt(3)+xbar;
yg1=-0.5*dy/sqrt(3)+ybar;
yg2=0.5*dy/sqrt(3)+ybar;

% incoming wave potential
phi0=exp(nu*(ybar-complex(0,1)*xbar));
phi0n=nu*(n2-complex(0,1)*n1).*phi0;

% contributions to integral equation, rhs stores rhs, lhs stores lhs
for i=1:N
    for j=1:N
        % rhs, log(r) term with 2pts Gauss quadrature
        xa1=xg1(j)-xbar(i);
        xa2=xg2(j)-xbar(i);
        ya1=yg1(j)-ybar(i);
        ya2=yg2(j)-ybar(i);
        ra1=sqrt(xa1*xa1+ya1*ya1);
        ra2=sqrt(xa2*xa2+ya2*ya2);
        g0=(log(ra1)+log(ra2))*0.5;
        
        % all other terms with midpoint rule
        xa=xbar(j)-xbar(i);
        yb=ybar(j)+ybar(i);
        rb=sqrt(xa*xa+yb*yb);
        g1=-log(rb);
        zz=nu*(yb-complex(0,1)*xa);
        f1=-2*exp(zz)*(expint(zz)+log(zz)-log(-zz));
        f2=2*pi*exp(zz);
        g2=real(f1)+complex(0,1)*real(f2);
        gg(i,j)=(g0+g1+g2)*ds(j);
        
        % lhs
        arg0=imag(log((xm(j)-xbar(i)+complex(0,1)*(ym(j)-ybar(i)))/...
            (xp(j)-xbar(i)+complex(0,1)*(yp(j)-ybar(i)))));
        
        if j-i == 0
            arg0=pi;
        end
        
        arg1=imag(log((xm(j)-xbar(i)+complex(0,1)*(ym(j)+ybar(i)))...
            /(xp(j)-xbar(i)+complex(0,1)*(yp(j)+ybar(i)))));

        arg2=(n1(j)*(imag(f1)+complex(0,1)*imag(f2))+n2(j)...
            *(real(f1)+complex(0,1)*real(f2)) )*nu*ds(j);

        ss(i,j)=(arg0+arg1+arg2);
    end
end
rhs=gg*phi0n;
lhs=ss*phi0;
hold on
plot(in, real(rhs), 'k -', 'LineWidth',2)
g0=plot(in, real(lhs), 'r .', 'MarkerSize',[13]);
plot(in, imag(rhs), 'b -', 'LineWidth',2)
g1=plot(in, imag(lhs), 'r +', 'MarkerSize',[8]);
title('Real and imaginary parts of RHS and LHS of eq. 16', 'FontSize', 14)
xlabel('Section number m', 'FontSize', 14) % x-axis label
legend({'Re\{RHS\}','Re\{LHS\}','Im\{RHS\}','Im\{LHS\}'})
set(gca,'FontSize',14)


% Numerical solution to integral eq.
% lhs
% arg0=imag(log((xm(j)-bx(i)+complex(0,1)*(ym(j)-by(i)))...
% /(xp(j)-bx(i)+complex(0,1)*(yp(j)-by(i)))));
% if j-i == 0
%     arg0=-pi;
% end
% arg1=imag(log((xm(j)-bx(i)+complex(0,1)*(ym(j)+by(i)))...
% /(xp(j)-bx(i)+complex(0,1)*(yp(j)+by(i)))));
% help1=(n1(j)* (imag(f1)+complex(0,1)*imag(f2))+n2(j)...
% *(real(f1)+complex(0,1)*real(f2)) )*nu*ds(j);
% ss(i,j)=(arg0+arg1+help1);
% end
% end
% rhs=gg*n2;
% phi2=ss\rhs;
% 
% ff22=phi2.*n2.*ds;
% sff22=sum(ff22);
% 
% phi0=exp(nu*(by-complex(0,1)*bx));
% 
% AM2=complex(0,1)*(phi2.*(nu*n2-nu*complex(0,1)*n1)-n2).*phi0.*ds;
% AP2=complex(0,1)*(phi2.*(nu*n2+nu*complex(0,1)*n1)-n2).*conj(phi0).*ds;
% 
% sAM2=sum(AM2);
% sAP2=sum(AP2);