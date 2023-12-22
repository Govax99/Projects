clear; close all; clc;
mu = astroConstants(13);
%% DATA
[R0,V0] = kep2car(6773,0,0,0,0,0); %need kep2car and astroConstants to work
[a,en,i,OM,om,M] = car2kep(R0,V0);
E = fzero(@(E) E-en*sin(E)-M, M); %in the case is not circular
theta0 = 2*atan(sqrt((1+en)/(1-en))*tan(E/2));
%% Example 1
Ix = 0.04;
Iy = 0.06;
Iz = 0.08;
%orient_model(Ix,Iy,Iz);
J = [Ix 0 0; 0 Iy 0; 0 0 Iz];
invJ = inv(J);
n = sqrt(mu/norm(R0)^3);
omega0 = [1e-6 1e-6 n];
Tsim = 2*pi/n;
% compute initial attitude
A = [1 0 0; 0 1 0; 0 0 1];
qin = 0.5*sqrt(abs(1+A(1,1)+A(2,2)+A(3,3)));
q0 = [0.25/qin*(A(2,3)-A(3,2)),0.25/qin*(A(3,1)-A(1,3)),0.25/qin*(A(1,2)-A(2,1)),qin];
e_angles0 = [atan(-A(2,1)/A(2,2)),asin(A(2,3)),atan(-A(1,3)/A(3,3))];

%% Example 2
Ix = 0.06;
Iy = 0.08;
Iz = 0.04;
%orient_model(Ix,Iy,Iz);
J = [Ix 0 0; 0 Iy 0; 0 0 Iz];
invJ = inv(J);
n = sqrt(mu/norm(R0)^3);
omega0 = [1e-6 1e-6 n];
Tsim = 2*pi/n;
% compute initial attitude
A = [1 0 0; 0 1 0; 0 0 1];
qin = 0.5*sqrt(abs(1+A(1,1)+A(2,2)+A(3,3)));
q0 = [0.25/qin*(A(2,3)-A(3,2)),0.25/qin*(A(3,1)-A(1,3)),0.25/qin*(A(1,2)-A(2,1)),qin];
e_angles0 = [atan(-A(2,1)/A(2,2)),asin(A(2,3)),atan(-A(1,3)/A(3,3))];