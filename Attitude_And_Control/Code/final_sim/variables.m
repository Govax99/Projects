clear;
clc;

%% initalization of constants of the s/c
tic

m = 54; % mass of the s/c [kg]
l = 0.5; % lenght of the cube [m]

% with solar panels closed
Ixx = 1/6 * m * l^2 + 0.25; % [kg*m^2]
Iyy = 1/6 * m * l^2; % [kg*m^2]
Izz = 1/6 * m * l^2 + 0.5; % [kg*m^2]

J = [Ixx 0 0; 0 Iyy 0; 0 0 Izz];
invJ = inv(J);

Ir = 0.05; %[kg*m^2]

Jr = [Ir 0 0; 0 Ir 0; 0 0 Ir]; %[kg*m^2]
invJr = inv(Jr);

%% initialize environment
% AERO TORQUE
N_sur = 10;
N_i = [1 0 0; 0 1 0; 0 0 1; -1 0 0; 0 -1 0; 0 0 -1; 1 0 0; -1 0 0; 1 0 0; -1 0 0]'; %normal to surfaces
A_i = [2500 2500 2500 2500 2500 2500 2500 2500 2500 2500]'*10^-4; %area of surfaces
r_i = [25 0 0; 0 25 0; 0 0 25; -25 0 0; 0 -25 0; 0 0 -25; -25 50 0; -25 50 0; -25 -50 0; -25 -50 0]'*10^-2; %distance cp from cg to surfaces
Across = 0.825;
Cd = 2.1; %coefficent of drag (from orbital data)

% MAGNETIC
k = 13;
alfaG = 0;

Rt = 6378.1363;
[h0,hmax,rho0,H] = rho_atmo();

om_Earth = 7.2921e-5;
p = [0.05 0.05 0.05]'; %dipole moment, guess
[g,h] = ghnorm(2020); %if we use my function

%the date of the mission must be chosen for sun direction and alfaG0

date = date2jd([2021,2,23,0,0,0]); %provvisory

%% Initialize sensors and attitude determination
%star sensor
SstarsRef=[1/sqrt(3) 1 0; 1/sqrt(3) 0 1; 1/sqrt(3) 0 0];
% rng(1);
% SstarsRef=rand(3,3);
%normalize matrix columns
% for i=1:3
%     SstarsRef(:,i)=SstarsRef(:,i)/sqrt(dot(SstarsRef(:,i),SstarsRef(:,i)));
% end
resStar =0.001;
alfaX =deg2rad(0.001);
alfaY =deg2rad(0.001);
alfaZ =deg2rad(0.001);
maxStar = 0.01;
varStar = maxStar^2*ones(1,3);
%sun sensor
resSun = 0.01;
angleSun = deg2rad(0.01);
varSun = angleSun^2;
%attitude determination
varVec = [varSun,varStar];
a_vec=1./varVec;

%% initial conditions of the dynamics
omegax0 = 0;
omegay0 = 0;
omegaz0 = 0;
omega0 = [omegax0 omegay0 omegaz0]';

omegaxr0 = 0;
omegayr0 = 0;
omegazr0 = 0;
omegar0 = [omegaxr0 omegayr0 omegazr0]';

m_MT = 200; %[A*m]
cutoff_omega = 0.1;

%% inital conditions of the kinematics
q0 = [0 0 0 1]';

%% initial condition of the orbit of the s/c
% it all depends on how the inital conditions are given
[R0,V0] = kep2car(0.8123e4,0.1789,50.3529*360/(2*pi),pi,pi,pi); % need kep2car and astroConstants to work
[a,en,i,OM,om,M] = car2kep(R0,V0);

E = fzero(@(E) E-en*sin(E)-M, M); % elliptical case
theta0 = 2*atan(sqrt((1+en)/(1-en))*tan(E/2)); % true anomaly

mu = astroConstants(13);
T = 2*pi*sqrt(a^3/mu);
%% Earth pointing
R1 = [cos(OM) sin(OM) 0
    -sin(OM) cos(OM) 0
    0 0 1];
R2 = [1 0 0
    0 cos(i) sin(i)
    0 -sin(i) cos(i)];
R3 = [cos(om) sin(om) 0
    -sin(om) cos(om) 0
    0 0 1];

A = R3 * R2 * R1;

% q0 = dcm2quat(A); % q = a + bi + cj + dk
% q0 = [q0(2) q0(3) q0(4) q0(1)]; % q = ai + bj + ck + d

%% initial conditions of the PID controller
kp = 20;
ki = 0;
kd = 0;
kpp = 10;
kip = 1;
kdp = 0;
kp1 = 0.2;
kp2 = 2;
komega = 0.01;

toc