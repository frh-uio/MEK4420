clear;
clc;
% characteristics of geometry

% load -ascii box1.dat
load -ascii box2.dat
% load -ascii box3.dat

% xm=box1(:,1);
% ym=box1(:,2);
% xp=box1(:,3);
% yp=box1(:,4);

xm=box2(:,1);
ym=box2(:,2);
xp=box2(:,3);
yp=box2(:,4);

% xm=box3(:,1);
% ym=box3(:,2);
% xp=box3(:,3);
% yp=box3(:,4);

L = 1;
D = 1;


N=30;
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

X2_direct = 0;
X2_FK = 0;
X2_HK1 = 0;
X2_HK2 = 0;
sff22 = 0;
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
    phi0=exp(K(k)*(ybar-complex(0,1)*xbar));
    rhs=gg*n2;
    phi2=ss\rhs;
    
    rhsD=-2*pi*phi0;
    phiD=ss\rhsD;
    
    % calculation of added mass a22 and damping b22
    X22_direct=phiD.*n2.*ds;
    X2_direct=[X2_direct,abs(sum(X22_direct))];
   
    AM2=complex(0,1)*(phi2.*(K(k)*n2-K(k)*complex(0,1)*n1)-n2).*phi0.*ds;

    X22_FK = (L*exp(-K(k)*D))*sin(K(k)*L/2.0)/(K(k)*L/2.0);
    X2_FK = [X2_FK, abs(sum(X22_FK))];
    
    X22_HK1 = -phi0.*(n2+complex(0,1)*phi2*K(k)*-K(k).*phi2.*n2).*ds;
    X2_HK1 = [X2_HK1, abs(sum(real(X22_HK1)))];
    
    X22_HK2 = AM2;
    X2_HK2 = [X2_HK2, abs(sum(X22_HK2))];
    
    %mass-damping force:
    ff22=phi2.*n2.*ds;
    sff22 = [sff22,sum(ff22)];
end
X2_direct = X2_direct(2:end);
X2_FK = X2_FK(2:end);
X2_HK1 = X2_HK1(2:end);
X2_HK2 = X2_HK2(2:end);
sff22 = sff22(2:end);

% figure()
% plot(K, X2_direct, 'b -', 'LineWidth', 2)
% hold on
% plot(K, X2_FK, 'r -', 'LineWidth', 2)
% hold on
% plot(K, X2_HK1, 'k -', 'LineWidth', 2)
% hold on
% plot(K, X2_HK2, 'g -', 'LineWidth', 2)
% 

% xlabel('\omega^2 D / g', 'FontSize', 20) % x-axis label
% ylabel('|X_2| / \rho g', 'FontSize', 20) % x-axis label
% legend('X_2^{Direct}','X_2^{FK}','X_2^{HK1}','X_2^{HK2}')
% set(gca,'FontSize',14)

% plot of body response
xi2 = abs(X2_direct./(L - K.*L -K.*sff22 ));
b22f=(abs(X2_FK)).^2; % Damping from FK
xiFK = abs(X2_FK./(L - K.*L + K.*complex(0,1).*b22f));
a22=real(sff22);
xia22= abs(X2_FK./(L - K.*L + K.*complex(0,1).*b22f) - a22);

figure(1);
plot(K, xi2, 'k -', 'LineWidth', 2);
hold on
plot(K, xiFK, 'r .', 'markersize', [20]);
plot(K, xia22, 'b -', 'LineWidth', 2);
title('Response as a function of wavenumber', 'FontSize', 20)
legend('|\xi_2|/ A', '|\xi_2^{FK}|/ A', '|\xi_2 (a_{22})|/ A')
xlabel('\omega^2 D /g', 'FontSize', 20)
set(gca,'FontSize',20)

