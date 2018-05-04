clear;
clc;
% characteristics of geometry

% load -ascii box1.dat
% load -ascii box2.dat
load -ascii box3.dat

% xm=box1(:,1);
% ym=box1(:,2);
% xp=box1(:,3);
% yp=box1(:,4);

% xm=box2(:,1);
% ym=box2(:,2);
% xp=box2(:,3);
% yp=box2(:,4);

xm=box3(:,1);
ym=box3(:,2);
xp=box3(:,3);
yp=box3(:,4);

L = 0.1;
D = 1;


N=25;
M=60;
in=linspace(1,N,N);
K=linspace(0.01,2,M);

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

sff22 = 0;
sAM2 = 0;
sAP2 = 0;
bb_22_approx = L^2;
% contributions to integral equation, rhs stores rhs, lhs stores lhs
for k=1:M
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
            zz=K(k)*(yb-complex(0,1)*xa);
            f1=-2*exp(zz)*(expint(zz)+log(zz)-log(-zz));
            f2=2*pi*exp(zz);
            g2=real(f1)+complex(0,1)*real(f2);
            gg(i,j)=(g0+g1+g2)*ds(j);

            % lhs
            arg0=imag(log((xm(j)-xbar(i)+complex(0,1)*(ym(j)-ybar(i)))/...
                (xp(j)-xbar(i)+complex(0,1)*(yp(j)-ybar(i)))));

            if j-i == 0
                arg0=-pi;
            end

            arg1=imag(log((xm(j)-xbar(i)+complex(0,1)*(ym(j)+ybar(i)))...
                /(xp(j)-xbar(i)+complex(0,1)*(yp(j)+ybar(i)))));

            arg2=(n1(j)*(imag(f1)+complex(0,1)*imag(f2))+n2(j)...
                *(real(f1)+complex(0,1)*real(f2)) )*K(k)*ds(j);

            ss(i,j)=(arg0+arg1+arg2);
        end
    end
    rhs=gg*n2;
    phi2=ss\rhs;
    
    % calculation of added mass a22 and damping b22
    ff22=phi2.*n2.*ds;
    sff22=[sff22,sum(ff22)];
    
    phi0=exp(K(k)*(ybar-complex(0,1)*xbar));

    AM2=complex(0,1)*(phi2.*(K(k)*n2-K(k)*complex(0,1)*n1)-n2).*phi0.*ds;
    AP2=complex(0,1)*(phi2.*(K(k)*n2+K(k)*complex(0,1)*n1)-n2).*conj(phi0).*ds;
   
    sAM2=[sAM2,sum(AM2)];
   
    sAP2=[sAP2,sum(AP2)];
    
    % approximate solution using Froude-Krylov approximation
    b_22_approx = (4*exp(-2*K(k)*D)*(sin((K(k)*L)/2))^2 )/(K(k)*K(k));
    bb_22_approx = [bb_22_approx, b_22_approx];
end

a22=real(sff22);
a22=a22(2:end);
b22=-imag(sff22);
b22=b22(2:end);

% b22_2=0.5*( (abs(sAP2)).^2 + (abs(sAM2)).^2 );
% b22_2_full = b22_2(2:end);
% 
% bb_22_approx = bb_22_approx(2:end); 

figure()
plot(K, a22, 'b -', 'LineWidth', 2)
ylim([-0.0 6])
hold on
plot(K, b22, 'r .', 'LineWidth', 2)
% hold on
% plot(K, b22_2_full, 'k +', 'LineWidth', 1)
% hold on
% plot(K, bb_22_approx, 'g -', 'LineWidth', 2)

%title('Added mass and damping', 'FontSize', 16)
title('dampings coefficient b_{22} vs wavenumber K', 'FontSize', 16)
xlabel('\omega^2 D / g', 'FontSize', 20) % x-axis label
%ylabel('\phi_2', 'FontSize', 18) % x-axis label
legend('a_{22}/ \rho D^2','b_{22}/ \rho \omega D^2')
%legend('b_{22}/ \rho \omega D^2','b_{22}^E/ \rho \omega D^2', 'b_{22}^{FK}/ \rho \omega D^2')
set(gca,'FontSize',14)
