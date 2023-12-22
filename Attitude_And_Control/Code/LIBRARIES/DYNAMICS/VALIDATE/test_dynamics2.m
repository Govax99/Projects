%% clear
clear;
clc;
close all;
mode = 'wheel';
%% data
Ix = 0.07;
Iy = 0.0109;
Iz = 0.0504;
Mx = 0;
My = 0;
Mz = 0;
J = [Ix 0 0; 0 Iy 0; 0 0 Iz];
invJ = inv(J);
Ir = 0.005;
omegar = 0;
%% general motion
omegax0 = 0.45;
omegay0 = 0.52;
omegaz0 = 0.55;
omega0 = [omegax0, omegay0, omegaz0]';
%verified
%% spin stability 1
Ix = 0.01;
Iy = 0.05;
Iz = 0.07;
J = [Ix 0 0; 0 Iy 0; 0 0 Iz];
invJ = inv(J);
omegax0 = 2*pi;
omegay0 = 0.52;
omegaz0 = 0.55;
omega0 = [omegax0, omegay0, omegaz0]';
%verified + ellipsoid
figure;
T0 = 0.5 * omega0' * J * omega0;
aT = sqrt(2*T0/Ix);
bT = sqrt(2*T0/Iy);
cT = sqrt(2*T0/Iz);
[XT,YT,ZT]=ellipsoid(0,0,0,aT,bT,cT,80);
surf(XT,YT,ZT,'FaceAlpha',0.7,'FaceColor','red','EdgeColor','none');
hold on;
hn = norm(J * omega0);
aH = hn/Ix;
bH = hn/Iy;
cH = hn/Iz;
[XH,YH,ZH]=ellipsoid(0,0,0,aH,bH,cH,80);
surf(XH,YH,ZH,'FaceAlpha',0.7,'FaceColor','blue','EdgeColor','none');
legend({'Kinetic energy ellipsoid','Momentum ellipsoid'});
axis equal;
%% spin stability 2
Ix = 0.01;
Iy = 0.05;
Iz = 0.07;
J = [Ix 0 0; 0 Iy 0; 0 0 Iz];
invJ = inv(J);
omegax0 = 0.45;
omegay0 = 2*pi;
omegaz0 = 0.55;
omega0 = [omegax0, omegay0, omegaz0]';
%verified
figure;
T0 = 0.5 * omega0' * J * omega0;
aT = sqrt(2*T0/Ix);
bT = sqrt(2*T0/Iy);
cT = sqrt(2*T0/Iz);
[XT,YT,ZT]=ellipsoid(0,0,0,aT,bT,cT,80);
surf(XT,YT,ZT,'FaceAlpha',0.7,'FaceColor','red','EdgeColor','none');
hold on;
hn = norm(J * omega0);
aH = hn/Ix;
bH = hn/Iy;
cH = hn/Iz;
[XH,YH,ZH]=ellipsoid(0,0,0,aH,bH,cH,80);
surf(XH,YH,ZH,'FaceAlpha',0.7,'FaceColor','blue','EdgeColor','none');
legend({'Kinetic energy ellipsoid','Momentum ellipsoid'});
axis equal;
%% spin stability 3
Ix = 0.01;
Iy = 0.05;
Iz = 0.07;
J = [Ix 0 0; 0 Iy 0; 0 0 Iz];
invJ = inv(J);
omegax0 = 0.45;
omegay0 = 0.52;
omegaz0 = 2*pi;
omega0 = [omegax0; omegay0; omegaz0];
%verified
figure;
T0 = 0.5 * omega0' * J * omega0;
aT = sqrt(2*T0/Ix);
bT = sqrt(2*T0/Iy);
cT = sqrt(2*T0/Iz);
[XT,YT,ZT]=ellipsoid(0,0,0,aT,bT,cT,80);
surf(XT,YT,ZT,'FaceAlpha',0.7,'FaceColor','red','EdgeColor','none');
hold on;
hn = norm(J * omega0);
aH = hn/Ix;
bH = hn/Iy;
cH = hn/Iz;
[XH,YH,ZH]=ellipsoid(0,0,0,aH,bH,cH,80);
surf(XH,YH,ZH,'FaceAlpha',0.7,'FaceColor','blue','EdgeColor','none');
legend({'Kinetic energy ellipsoid','Momentum ellipsoid'});
axis equal;
%% dual spin
omegax0 = 1e-6;
omegay0 = 1e-6;
omegaz0 = 0.02;
omega0 = [omegax0, omegay0, omegaz0]';
if strcmp(mode,'wheel')
    omegar = 2*pi;
end
T_w = [0 0 0.0004]';