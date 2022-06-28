clc, clear, close all

load('output.mat')

Fx = F(1,:).*sin(beta).*cos(phi);
Fy = F(1,:).*sin(beta).*sin(phi);
Fz = F(1,:).*cos(beta);
% s = convhull(Fx,Fy,Fz);

figure
ax = axes();
xlim(ax, [-70 70]);
ylim(ax, [-70 70]);
zlim(ax, [-10 80]);
axis equal
hold(ax, 'on')
plot3(Fx,Fy,Fz,'.r')
% plot3(Fx,-Fy,Fz,'.r')
xlabel('Fx [N]')
ylabel('Fy [N]')
zlabel('Fz [N]')

figure
ax1 = axes();
xlim(ax1, [-70 70]);
ylim(ax1, [-70 70]);
zlim(ax1, [-10 80]);
axis equal
hold(ax1, 'on')
% trisurf(s,Fx,Fy,Fz)%,'.r','FaceColor','Cyan')
hold on, grid on
% surf(Fx(s),-Fy(s),Fz(s))%,'.r')
% plot3(Fx(surf),Fy(surf),Fz(surf))
% plot3(Fx(surf),-Fy(surf),Fz(surf))
xlabel('Fx [N]')
ylabel('Fy [N]')
zlabel('Fz [N]')