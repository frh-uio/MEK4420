clear;
clc;

% box of draught D and length (width) L, with D as unit length
% left vertical side
D=1;
L=0.1;
Nside=10;
Nbott=5;
N=Nside+Nbott+Nside;
dy=D/Nside;
dx=L/Nbott;
xp(1)=-L/2;
xm(1)=-L/2;
yp(1)=-dy;
ym(1)=-dy*(1-1);

coord=[xm(1),ym(1),xp(1),yp(1)];

for i=2:Nside;
    xp(i)=-L/2;
    xm(i)=-L/2;
    yp(i)=-dy*i;
    ym(i)=-dy*(i-1);
    coord2=[xm(i),ym(i),xp(i),yp(i)];
    coord=[coord;coord2];
end
for i=1+Nside:Nside+Nbott;
    i1=i-Nside;
    xp(i)=-L/2+dx*i1;
    xm(i)=-L/2+dx*(i1-1);
    yp(i)=-D;
    ym(i)=-D;
    coord2=[xm(i),ym(i),xp(i),yp(i)];
    coord=[coord;coord2];
end
for i=1+Nside+Nbott:Nside+Nbott+Nside;
    i1=i-Nside-Nbott;
    xp(i)=L/2;
    xm(i)=L/2;
    yp(i)=-D+dy*i1;
    ym(i)=-D+dy*(i1-1);
    coord2=[xm(i),ym(i),xp(i),yp(i)];
    coord=[coord;coord2];
end
save -ascii box3.dat coord;
% plot geometry
xa=xm;
xa=[xm,xp(N)];
ya=ym;
ya=[ym,yp(N)];
hold on
axis ([-1.1 1.1 -1.1 0.1])
h1=plot(xa,ya,'k .','MarkerSize',[25])
plot(xa,ya,'k -')
title('Discretization of a rectangular geometry', 'FontSize',20)
xlabel(['Length L = ' num2str(L)],'FontSize',20) % x-axis label
ylabel(['Draught D = ' num2str(D)], 'FontSize',20) % y-axis label
% axis off
% get(h1);
% set(h1,'MarkerSize',[25])
% set(gca,'FontSize',20)