clear;
clc;

%K=1.2;

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

N=25;
M=60;
in=linspace(1,N,N);
K=linspace(0.01,2,M);

% characteristics of geometry
dx=xp-xm;
dy=yp-ym;
ds=((dx).^2+(dy).^2).^(1/2);

% Midpoints
xbar=0.5*(xm+xp);
ybar=0.5*(ym+yp);

n1=-(dy)./ds;
n2=(dx)./ds;

% incoming wave potential
X2 = 0;
% contributions to integral equation, rhs stores rhs, lhs stores lhs
for k=1:M
    for i=1:N
        for j=1:N
        % rhs
        xa=xbar(j)-xbar(i);
        yb=ybar(j)+ybar(i);
        zz=K(k)*(yb-complex(0,1)*xa);
        f1=-2*exp(zz)*(expint(zz)+log(zz)-log(-zz));
        f2=2*pi*exp(zz);
        % lhs
        arg0=imag(log((xm(j)-xbar(i)+complex(0,1)*(ym(j)-ybar(i)))/(xp(j)-xbar(i)+complex(0,1)*(yp(j)-ybar(i)))));
        if j-i == 0
            arg0=-pi;
        end
        arg1=imag(log((xm(j)-xbar(i)+complex(0,1)*(ym(j)+ybar(i)))/(xp(j)-xbar(i)+complex(0,1)*(yp(j)+ybar(i)))));
        help1=(n1(j)* (imag(f1)+complex(0,1)*imag(f2))+n2(j)*(real(f1)+complex(0,1)*real(f2)) )*K(k)*ds(j);
        ss(i,j)=(arg0+arg1+help1);
        end
    end
    phi0=exp(K(k)*(ybar-complex(0,1)*xbar));
    rhsD=-2*pi*phi0;
    phiD=ss\rhsD;
    
    %exciting force
    XX2=phiD.*n2.*ds;
    sXX2=sum(XX2);

    X2 = [X2, abs(sXX2)];
end

X2 = X2(2:end);
plot(K, X2, 'r -', 'LineWidth', 1)
title('Exciting force X_2 in heave, L/D=0.1', 'FontSize', 14)
%xlabel('Number of segments', 'FontSize', 14) % x-axis label
xlabel('\omega^2 D / g', 'FontSize', 22) % x-axis label
%ylabel('\phi_D', 'FontSize', 18) % x-axis label
ylabel('|X_2| / \rho g', 'FontSize', 22)
%legend({'Re\{\phi_D\}','Im\{\phi_D\}'})
set(gca,'FontSize',14)
