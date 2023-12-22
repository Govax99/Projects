clear; close all; clc;
om_Earth = 0.000072921;
mu = 398600;
Rt = 6378.1363;
%% Data (aero-torque cuboid l = 10)
[h0,hmax,rho0,H] = rho_atmo();
rho = 1.215;
Cd = 4;

R0 = [7000 0 0]';
V0 = [0 7.546 0]';

Ix = 0.04;
Iy = 0.06;
Iz = 0.08;
J = [Ix 0 0; 0 Iy 0; 0 0 Iz];
invJ = inv(J);

omega0 = [1e-6 1e-6 0.5]';
A0 = [1 0 0; 0 1 0; 0 0 1];

%% Data (magnetic torque)
k = 5;
alfaG = 0;
a = 6371.2; % reference earth radius
m = [0.1 0.1 0.1]';
[g,h] = ghnorm(2015);
%% Data (Solar pressure)
omega0 = [0.45 0.52 0.55]';
A0 = [1 0 0; 0 1 0; 0 0 1];
% dot product can be element wise
N_sur = 10;
N_i = [1 0 0; 0 1 0; 0 0 1; -1 0 0; 0 -1 0; 0 0 -1; 1 0 0; -1 0 0; 1 0 0; -1 0 0]';
A_i = [600 600 400 600 600 400 1200 1200 1200 1200]'*10^-4;
r_i = [10 0 0; 0 10 0; 0 0 15; -10 0 0; 0 -10 0; 0 0 -15; 0 60 0; 0 60 0; 0 -60 0; 0 -60 0]'*10^-2;
Ix = 100.9*10^-2;
Iy = 25.1*10^-2;
Iz = 91.6*10^-2;
% ATT: for rotating solar panel learn how to rotate inertia matrix and sum
% Jtot = Jbody + A_IP * Jpanel * A_IP euler transport?
J = [Ix 0 0; 0 Iy 0; 0 0 Iz];
invJ = inv(J);
T_sun = 365.3*24*3600;
n_sun = 2*pi/T_sun;
eps = 23.45*pi/180;
R_s = 149600000;
Fe = 1358;
c = 3*10^8;
P = Fe/c; %transform in km?
rho_s = ones(10,1)*0.1;
rho_d = ones(10,1)*0.1;
rho_s(7:10,1) = rho_s(7:10,1)*5;