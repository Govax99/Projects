% Spacecraft Guidance and Navigation (2021/2022)
% Assignment # 2
% Author: Davide Zamblera

%% Section 1
clearvars; close all; clc;
cspice_furnsh('assignment02.tm')
addpath("simulator_ex3\")
set(0, 'defaulttextinterpreter','latex');  
set(0, 'defaultAxesTickLabelInterpreter','latex');  
set(0, 'defaultLegendInterpreter','latex');

% define target
l = 10;
target.l = l;
h = 5;
target.h = h;
d = 3;
target.d = d;
m = 213000;
J = m/12*diag([d^2 + h^2, l^2 + h^2, l^2 + d^2]);
target.J = J;
R=diag([10 10 10]);
target.R = R;

% define camera model
foc=30; % mm
dens=54; %pix/mm
b=1;
p0=[960;600]; % center pixel location;
Cframe=[1,0,0;
        0,0,-1;
        0,1,0];
target.Cframe = Cframe;
R=10;

Cam.f=foc;
Cam.d=dens;
Cam.p0=p0;
Cam.b=b;
Cam.Cframe=Cframe;
Cam.R=R;


% mean motion of GEO orbit computation:
T = 86400;
n=(2*pi)/T;
target.n = n;



%% Compute all measurements:

% initial time definition:
t0='2023-04-01T14:55:12.023';
et0=cspice_str2et(t0);
Tobs = 86400;
t_samp = 1;
tt_et = et0:t_samp:(et0 + Tobs);
tt = 0:t_samp:Tobs;

r0 = [12 -60 0]'; %nominal-reference trajectory, true initial state
v0 = [1e-4, -2*n*r0(1), -1.2e-3]';
q0 = [0.674156352338764;
      0.223585877389611;
      0.465489474399161;
      0.528055032413102];
om0 = [-0.001262427155865;
        0.001204540074343;
       -0.000039180139156];

x0 = [r0; v0];
att0 = [q0; om0];
target.attitude0 = att0;

% try out new stm

xnom = zeros(length(tt_et),6);
xnom(1,:) = x0';
for i = 2:length(tt_et)
    PHI = STM(et0,tt_et(i),target);
    xnom(i,:) = (PHI*x0)';
end

opt = odeset('RelTol',1e-12,'AbsTol',1e-12);
% [~,xnom] = ode113(@ode_approach,tt_et,x0,opt,target);
[~,att_nom] = ode113(@ode_attitude,tt_et,att0,opt,target);

figure;
plot3(0,0,0,'*');
hold on;
plot3(xnom(:,1),xnom(:,2),xnom(:,3));
title('Plot of the relative position of the chaser w.r.t the target in LVLH frame','Interpreter','latex');
xlabel('X [km]','Interpreter','latex');
ylabel('Y [km]','Interpreter','latex');
zlabel('Z [km]','Interpreter','latex');
grid on;
axis equal;
legend({'Target position','Chaser trajectory'},'Interpreter','latex');

y = cell(length(tt),1);
visible = cell(length(tt),1);

vertices_hist_hor = NaN(length(tt),6);
vertices_hist_vert = NaN(length(tt),6);
%vertices_hist_disparity = zeros(length(tt),6);

for i = 1:length(tt)
    r = xnom(i,1:3)';
    q = att_nom(i,1:4)';
    meas = meas_sim_pvt(n,r,q,tt(i),et0,Cam);


    y{i} = meas.y;
    visible{i} = meas.visible;

    vertices_hist_hor(i,meas.visible) = meas.y(1,:);
    vertices_hist_vert(i,meas.visible) = meas.y(2,:);
end
meas.y = y;
meas.visible = visible;

figure;
plot(vertices_hist_hor,vertices_hist_vert,'.')
grid on;
title('Camera frame, pixel measurement history','Interpreter','latex');
xlabel('Horizontal coordinate [px]','Interpreter','latex');
ylabel('Vertical coordinate [px]','Interpreter','latex');

str = strings(1,8);
for i = 1:8
    str(i) = sprintf('Vertex %d', i);
end
legend(str,'Interpreter','latex');

nvert = zeros(length(visible),1);
% Show number of visible vertices for various times
for i = 1:length(visible)
    nvert(i) = length(visible{i});
end

figure;
plot(tt_et/cspice_spd,nvert,'k.','MarkerSize',0.1);
ylabel('$n^{\circ}$ visible vertices','Interpreter','latex');
xlabel('t [MJD2000]','Interpreter','latex');
grid on;

%% Apply Extended Kalman Filter
mu0 = [15.792658268071492; -59.044939772661586; 3.227106250277039;
       -0.053960274403210; -0.053969644762889; -0.089140748762173];
P0 = diag([10 10 10 0.1 0.1 0.1]);
hist_EKF = EKF(mu0, P0, tt, meas, target, Cam);

errorEKF = (squeeze(hist_EKF.mu))' - xnom;

%extract 3sigma errors
boundaries_EKF = cell(6,1);
for i = 1:6
    boundaries_EKF{i} = zeros(length(tt),1);
    for j = 1:length(tt)
        boundaries_EKF{i}(j) = 3*sqrt(hist_EKF.P(i,i,j));
    end
end
plotlabel = {'error x [m]','error y [m]','error z [m]','error $v_x$ [m/s]','error $v_y$  [m/s]', 'error $v_z$  [m/s]'};
for i = 1:6
    figure(3+i);
    plot(tt_et/cspice_spd(),errorEKF(:,i),'r')
    hold on;
    plot(tt_et/cspice_spd(),boundaries_EKF{i},'k--')
    plot(tt_et/cspice_spd(),-boundaries_EKF{i},'k--')
    xlabel('t [MJD2000]','Interpreter','latex');
    ylabel(plotlabel{i},'Interpreter','latex');
    title('Error EKF Filter, LVLH frame @target')
    grid on;
    legend({'EKF filter','3$\sigma_{EKF}$'},'Interpreter','latex')
end

%% Apply Unscented Kalman Filter
hist_UKF = UKF(mu0, P0, tt, meas, target, Cam);


errorUKF = (squeeze(hist_UKF.mu))' - xnom;

%extract 3sigma errors
boundaries_UKF = cell(6,1);
for i = 1:6
    boundaries_UKF{i} = zeros(length(tt),1);
    for j = 1:length(tt)
        boundaries_UKF{i}(j) = 3*sqrt(hist_UKF.P(i,i,j));
    end
end

for i = 1:6
    figure(9+i);
    plot(tt_et/cspice_spd(),errorUKF(:,i),'b')
    hold on;
    plot(tt_et/cspice_spd(),boundaries_UKF{i},'k--')
    plot(tt_et/cspice_spd(),-boundaries_UKF{i},'k--')
    grid on;
    xlabel('t [MJD2000]','Interpreter','latex');
    ylabel(plotlabel{i},'Interpreter','latex');
    title('Error UKF Filter, LVLH frame @target')
    grid on;
    legend({'UKF filter','3$\sigma_{UKF}$'},'Interpreter','latex')
end


%% End of File functions

function f = ode_approach(t,s,target)
%ODE_APPROACH Clohessy Wiltshire rendez vous equations rhs
%
% PROTOTYPE:
%     f = ode_approach(t,s,target)
%
% INPUT:
%    t[dim]            time [s]
%    s[dim]            state [km; km/s]
%    target[struct]    structure with relevant informations about the
%                      target
%                       -target.n (mean motion of target)
%
% OUTPUT:
%    f[dim]           rhs of the system [-]
%
% CONTRIBUTORS:
%   Davide Zamblera - Politecnico di Milano - Space eng.
%
    n = target.n; %mean motion of target
    r = s(1:3);
    v = s(4:6);

    f = zeros(6,1);

    f(1:3) = v;
    f(4) = 3*n^2*r(1) + 2*n*v(2);
    f(5) = -2*n*v(1);
    f(6) = -n^2*r(3);
end

function PHI = STM(t0,t,target)
%STM Clohessy Wiltshire state transition matrix
%
% PROTOTYPE:
%     PHI = STM(t0,t,target)
%
% INPUT:
%    t0[dim]           initial time [s]
%    t[dim]            time [s]
%    target[struct]    structure with relevant informations about the
%                      target
%                       -target.n (mean motion of target)
%
% OUTPUT:
%    PHI[6x6]          state transition matrix [-]
%
% CONTRIBUTORS:
%   Davide Zamblera - Politecnico di Milano - Space eng.
%
    n = target.n;

    A = zeros(6);
    A(1:3,4:6) = eye(3);
    A(4,1) = 3*n^2;
    A(6,3) = -n^2;
    A(4:6,4:6) = [0 2*n 0;
                 -2*n 0 0;
                  0 0 0];
    PHI = expm(A*(t-t0));
end

function f = ode_attitude(t,s,target)
%ODE_ATTITUDE rigid body motion attitude rhs
%
% PROTOTYPE:
%     f = ode_attitude(t,s,target)
%
% INPUT:
%    t[dim]            time [s]
%    s[dim]            state [km; km/s]
%    target[struct]    structure with relevant informations about the
%                      target
%                       -target.n (mean motion of target)
% OUTPUT:
%    f[dim]           rhs of the system [-]
%
% CONTRIBUTORS:
%   Davide Zamblera - Politecnico di Milano - Space eng.
%
    J = target.J;
    q = s(1:4);
    om = s(5:7);

    f = zeros(7,1);

    M = [0     -om(1) -om(2) -om(3);
         om(1)  0      om(3) -om(2);
         om(2) -om(3)  0      om(1);
         om(3)  om(2)  -om(1)  0];
    f(1:4) = 1/2*M*q;
    f(5:7) = J\(cross(-om,J*om));

end

function hist = EKF(mu0, P0, t, meas, target, Cam)
%EKF extended kalman filter main function
%
% PROTOTYPE:
%     hist = EKF(mu0, P0, t, meas, target, Cam)
%
% INPUT:
%    mu0[6x1]          initial mean state [-]
%    P0[6x6]           initial covariance [-]
%    t[dim]            time from initial epoch [s]
%    meas[cell]        measurements and visibility for each epoch [-]
%    target[struct]    structure with relevant informations about the
%                      target
%
%    Cam[struct]    structure with relevant informations about the
%                   pinhole camera
%
% OUTPUT:
%    hist[struct]         estimation history
%                      -hist.mu: mean state [3rd dimension time]
%                      -hist.P: covariances [3rd dimension time]
%
% CONTRIBUTORS:
%   Davide Zamblera - Politecnico di Milano - Space eng.
%
    hist.mu(:,:,1) = mu0;
    hist.P(:,:,1) = P0;
    attitude0 = target.attitude0;
    R = target.R;

    %compute all attitude in advance (from link with target)
    opt = odeset('RelTol',1e-12,'AbsTol',1e-12);
    [~,attitude] = ode113(@ode_attitude,t,attitude0,opt,target);
    q = attitude(:,1:4);
    %om = attitude(:,5:7);


    for i = 1:length(t)-1
        % Prediction step
        PHI = STM(t(i),t(i+1),target);
        x_min = PHI*hist.mu(:,:,i);
        P_min = PHI*hist.P(:,:,i)*PHI';

        % Measurement step
        n_vis = length(meas.visible{i+1});
        y_meas = reshape(meas.y{i+1}, [3*n_vis 1]);
        [y_min, H] = measurement_function(t(i+1),x_min,q(i+1,:),meas.visible{i+1},target,Cam);

        Rmeas = zeros(3*n_vis);
        for j = 1:n_vis
            Rmeas(1 + 3*(j-1):3*j,1 + 3*(j-1):3*j) = R;
        end


        % update step
        K = P_min*H'*inv(H*P_min*H' + Rmeas);
        x_plus = x_min + K*(y_meas - y_min);
        P_plus = (eye(6) - K*H)*P_min;
        P_plus = (P_plus + P_plus')/2; %enforce symmetry of covariance matrix

        hist.mu(:,:,i+1) = x_plus;
        hist.P(:,:,i+1) = P_plus;

    end


end

function [measure_sim, H] = measurement_function(time,state,q,visible,target,Cam)
%MEASUREMENT_FUNCTION compute measurements given the state (optionally
%compute the jacobian of the measurement function)
%
% PROTOTYPE:
%     [measure_sim, H] = measurement_function(time,state,q,visible,target,Cam)
%
% INPUT:
%    time[dim]         epoch of the measurement [s]
%    state[6x1]        state at the measurement time [-]
%    q[dim]            attitude at the measurement time [-]
%    visible[dim]      information about the visibility of the vertices [-]
%    target[struct]    structure with relevant informations about the
%                      target
%
%    Cam[struct]    structure with relevant informations about the
%                   pinhole camera
%
% OUTPUT:
%    measure_sim[3*n_vis x 1]  pixel coordinates and disparity (every column is a visible vertex)
%    H[3*n_vis x 3*n_vis]      jacobian of the measurement function [-]
%
% CONTRIBUTORS:
%   Davide Zamblera - Politecnico di Milano - Space eng.
%
    H = [];
    l = target.l;
    d = target.d;
    h = target.h;
    n = target.n;
    C_CL = Cam.Cframe;
    u0 = Cam.p0(1);
    v0 = Cam.p0(2);
    foc = Cam.f;
    dens = Cam.d;
    b = Cam.b;


    pol = [l/2 -d/2 -h/2;
         l/2 d/2 -h/2;
         l/2 d/2 h/2;
         l/2 -d/2 h/2;
        -l/2 -d/2 -h/2;
        -l/2 d/2 -h/2;
        -l/2 d/2 h/2;
        -l/2 -d/2 h/2];
    pol = pol'; %every column is a vertex in target frame @target




    C_TI = quat2dcm(q);
    C_LI = [cos(n*time), sin(n*time) 0;
            -sin(n*time), cos(n*time) 0
            0,            0,         1]; %matrix from target frame to lvlh - unknown
    pol_I = C_TI\pol;
    pol_LVLH = C_LI*pol_I;
    pol_LVLH_at_chaser = pol_LVLH - state(1:3);
    V = C_CL*pol_LVLH_at_chaser;
    X = V(1,visible);
    Y = V(2,visible);
    Z = V(3,visible);
    n_vis = length(visible);
    measure_sim = [u0 - dens*foc*Y./Z; v0 + dens*foc*X./Z; b*dens*foc./Z];
    measure_sim = reshape(measure_sim,[3*n_vis 1]);

    if nargout > 1
        % Compute jacobian of the measurement
        dhdX = zeros(3*n_vis);
        dXdx = zeros(3*n_vis,6);
        for j = 1:n_vis
            dhdX(1 + 3*(j-1):3*j,1 + 3*(j-1):3*j) = ...
                   [0             -dens*foc/Z(j)     dens*foc*Y(j)/Z(j)^2;
                    dens*foc/Z(j)  0                -dens*foc*X(j)/Z(j)^2;
                    0              0                -b*dens*foc/Z(j)^2];
            dXdx(1 + 3*(j-1):3*j,1:3) = -C_CL;
        end
        H = dhdX*dXdx;
    end

end

function [x,wm,wc] = sigma_points(mu,P,param)
% Compute sigma points for generalized unscented transform
%
% PROTOTYPE:
%   [x,wm,wc] = sigma_points(mu,P,param)
%
% INPUT:
%    mu           mean value
%    P            covariance matrix
%    param        (alfa,beta,k) choose collocation of sigma points
%
% INTERMEDIATE PARAMETERS
%    n            dimension of mu
%
% OUTPUT:
%	 x 	       sigma points colums of matrix x[n,2*(n+1)]
%    wm        weights for mean of propagated sigma values
%    wc        weights for covariance of propagated sigma values
%
% CONTRIBUTORS:
%   Davide Zamblera - Politecnico di Milano - Space eng.
%

    n = length(mu);
    if nargin < 3
        k = 0;
        alfa = 1e-3;
        beta = 2;
    else
        alfa = param.alfa;
        k = param.k;
        beta = param.beta;
    end
    c = alfa^2*(n+k);

    x = zeros(length(mu),2*n+1); %every column is a sigma-point
    x(:,1) = mu;
    L = sqrtm(c*P);
    x(:,1) = mu;
    for i = 1:2*n
        if i <= n
            x(:,i+1) = mu + L(:,i);
        else
            x(:,i+1) = mu - L(:,i-n);
        end
    end

    wm = zeros(1,2*n+1);
    wc = zeros(1,2*n+1);
    wm(1) = 1 - n/c;
    wc(1) = (2 - alfa^2 + beta) - n/c;
    wm(2:end) = 1/(2*c);
    wc(2:end) = 1/(2*c);


end


function hist = UKF(mu0, P0, t, meas, target, Cam)
%UKF unscented kalman filter main function
%
% PROTOTYPE:
%     hist = UKF(mu0, P0, t, meas, target, Cam)
%
% INPUT:
%    mu0[6x1]          initial mean state [-]
%    P0[6x6]           initial covariance [-]
%    t[dim]            time from initial epoch [s]
%    meas[cell]        measurements and visibility for each epoch [-]
%    target[struct]    structure with relevant informations about the
%                      target
%                      -target.attitude0: attitude at initial epoch
%                      -target.R: covariance of the measurements
%
%    Cam[struct]    structure with relevant informations about the
%                   pinhole camera
%
% OUTPUT:
%    hist[struct]         estimation history
%                      -hist.mu: mean state [3rd dimension time]
%                      -hist.P: covariances [3rd dimension time]
%
% CONTRIBUTORS:
%   Davide Zamblera - Politecnico di Milano - Space eng.
%
    hist.mu(:,:,1) = mu0;
    hist.P(:,:,1) = P0;
    attitude0 = target.attitude0;
    R = target.R;
    dim_state = length(mu0);

    %compute all attitude in advance (from link with target)
    opt = odeset('RelTol',1e-12,'AbsTol',1e-12);
    [~,attitude] = ode113(@ode_attitude,t,attitude0,opt,target);
    q = attitude(:,1:4);
    %om = attitude(:,5:7);


    for i = 1:length(t)-1
        % Prediction step
        [chi,wm,wc] = sigma_points(hist.mu(:,:,i),hist.P(:,:,i));
        n_sigma_points = size(chi,2);

        PHI = STM(t(i),t(i+1),target);
        for j = 1:n_sigma_points
            chi(:,j) = PHI*chi(:,j);
        end

        n_vis = length(meas.visible{i+1});
        dim_meas = 3*n_vis;
        gamma = zeros(3*n_vis,n_sigma_points);
        for j = 1:n_sigma_points
            measure_sim = measurement_function(t(i+1),chi(:,j),q(i+1,:),meas.visible{i+1},target,Cam);
            gamma(:,j) = measure_sim;
        end

        % Apply UT weights to compute means and covariances
        y_meas = reshape(meas.y{i+1}, [3*n_vis 1]);
        x_min = zeros(dim_state,1);
        y_min = zeros(dim_meas,1);
        P_min = zeros(dim_state);
        P_ee = zeros(dim_meas);
        P_xy = zeros(dim_state,dim_meas);
        for j = 1:n_sigma_points
            x_min = x_min + wm(j)*chi(:,j);
            y_min = y_min + wm(j)*gamma(:,j);
        end

        for j = 1:n_sigma_points
            P_min = P_min + wc(j)*(chi(:,j) - x_min)*(chi(:,j) - x_min)';
            P_ee = P_ee + wc(j)*(gamma(:,j) - y_min)*(gamma(:,j) - y_min)';
            P_xy = P_xy + wc(j)*(chi(:,j) - x_min)*(gamma(:,j) - y_min)';
        end

        Rmeas = zeros(3*n_vis);
        for j = 1:n_vis
            Rmeas(1 + 3*(j-1):3*j,1 + 3*(j-1):3*j) = R;
        end
        P_ee = P_ee + Rmeas;


        % Update state
        K = P_xy*inv(P_ee);
        x_plus = x_min + K*(y_meas - y_min);
        P_plus = P_min - K*P_ee*K';
        P_plus = (P_plus + P_plus')/2; %enforce symmetry of covariance matrix

        hist.mu(:,:,i+1) = x_plus;
        hist.P(:,:,i+1) = P_plus;

    end


end

