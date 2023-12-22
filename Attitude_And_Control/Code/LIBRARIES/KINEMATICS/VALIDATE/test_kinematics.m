clear; close all; clc;
addpath('my_math','dynamics','kinematic_transform','kinematics');
%% Initialize variables
Ix = 0.07;
Iy = 0.0504;
Iz = 0.0109;
J = [Ix 0 0; 0 Iy 0; 0 0 Iz];
invJ = inv(J);
C = pi/4; % C defined in [0.2 2pi]
omega0 = [C 0.1 0.1];
h0 = [Ix*omega0(1) Iy*omega0(2) Iz*omega0(3)]';
x0 = h0/norm(h0);
v = [0 0 1];
y0 = cross(v,x0)/norm(cross(v,x0));
z0 = cross(x0,y0)/norm(cross(x0,y0));
A0t = [x0' ;y0; z0];
A = A0t';
qin = 0.5*sqrt(abs(1+A(1,1)+A(2,2)+A(3,3)));
q0 = [0.25/qin*(A(2,3)-A(3,2)),0.25/qin*(A(3,1)-A(1,3)),0.25/qin*(A(1,2)-A(2,1)),qin];
e_angles0 = [atan(-A(2,1)/A(2,2)),asin(A(2,3)),atan(-A(1,3)/A(3,3))];