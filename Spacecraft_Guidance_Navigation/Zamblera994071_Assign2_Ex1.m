% Spacecraft Guidance and Navigation (2020/2021)
% Assignment # 2
% Author: Davide Zamblera

%% Section 1
clearvars; close all; clc;

cspice_furnsh('assignment02.tm');
%data definition
t_utc = '2022-11-11T19:08:49.824';

r0 = [6054.30795817484, -3072.03883303992, -133.115352431876]';
v0 = [4.64750094824087, 9.18608475681236, -0.62056520749034]';
y0 = [r0; v0];
mu0 = y0;

P0 = [+5.6e-3 +3.5e-3 -7.1e-4 0 0 0;
      +3.5e-3 +9.7e-3 +7.6e-4 0 0 0;
      -7.1e-4 +7.6e-4 +8.1e-4 0 0 0;
      0 0 0                   +2.8e-7 0 0;
      0 0 0                   0 +2.7e-7 0;
      0 0 0                   0 0 +9.6e-8];

GM = cspice_bodvrd('EARTH','GM',1);

% Obtain orbit period - definition of interval of integration

a = - GM/2 * 1/(0.5*norm(v0)^2 - GM/norm(r0));
h = cross(r0,v0);
e = cross(v0,h)/GM - r0/norm(r0);
e = norm(e);
rp = a*(1 - e);
ra = a*(1 + e);
T = 2*pi*sqrt(a^3/GM);

t0 = 0;
tspan = [t0 4.1*T];

% find extremes events trigger
opts = odeset('RelTol',1e-12,'AbsTol',1e-12,'Events',@extremeEvents);
[t,y,te,ye,ie] = ode113(@TBP,tspan,y0,opts,GM);

% disregard first value since it is the starting point
te = te(2:end);
ye = ye(2:end,:);
ne = length(te);

% Display of the orbit
figure;
plotEarth();
hold on; grid on; axis equal;
plot3(y(:,1),y(:,2),y(:,3))

% Display of the events
figure
plot(t/T,vecnorm(y(:,1:3),2,2),'k')
hold on;
plot(t/T,dot(y(:,1:3),y(:,4:6),2),'r')
plot(te/T,zeros(ne),'ko')
xlabel('t [T]','Interpreter','latex');
grid on;
legend('Radius magnitude [km]','Event Function value','Event points','Interpreter','latex');
% pericenter at increasing values of eventFun
% apocenter at decreasing values of eventFun


% linCov approach, solutions at various times are stored in 3rd dimension
mu_lin = zeros(6,1,ne);
P_lin = zeros(6,6,ne);

for i = 1:ne
    [mu, P] = linCov(mu0,P0,t0,te(i));
    mu_lin(:,:,i) = mu;
    P_lin(:,:,i) = P;
end

% determine sqrt(tr(P)) values, first column 'rr', second 'vv'
trP_lin = zeros(ne,2);

for i = 1:ne
    trP_lin(i,1) = sqrt(trace(P_lin(1:3,1:3,i)));
    trP_lin(i,2) = sqrt(trace(P_lin(4:6,4:6,i)));
end

% unscentedTransform approach, solutions at various times are stored in 3rd dimension

mu_ut = zeros(6,1,ne);
P_ut = zeros(6,6,ne);

for i = 1:ne
    [mu, P] = unscentedTransform(mu0,P0,t0,te(i));
    mu_ut(:,:,i) = mu;
    P_ut(:,:,i) = P;
end

trP_ut = zeros(ne,2);

for i = 1:ne
    trP_ut(i,1) = sqrt(trace(P_ut(1:3,1:3,i)));
    trP_ut(i,2) = sqrt(trace(P_ut(4:6,4:6,i)));
end

% monteCarlo simulation, solutions at various times are stored in 3rd dimension
rng('default');
mu_mc = zeros(6,1,ne);
P_mc = zeros(6,6,ne);
pop_mc = zeros(100,6,ne);

for i = 1:ne
    [mu, P, ~, pop] = monteCarlo(mu0, P0, t0, te(i));
    pop_mc(:,:,i) = pop;
    mu_mc(:,:,i) = mu;
    P_mc(:,:,i) = P;
end

trP_mc = zeros(ne,2);

for i = 1:ne
    trP_mc(i,1) = sqrt(trace(P_mc(1:3,1:3,i)));
    trP_mc(i,2) = sqrt(trace(P_mc(4:6,4:6,i)));
end

% plot of the results for various time instants (plot of pericenter and
% apocenter on 1 graph)

opts = odeset('RelTol',1e-12,'AbsTol',1e-12);
for i = 1:ne/2
    j1 = 2*(i-1) + 1; %index of apocenters
    figure;
    str_t = strcat('Orbit $n^{\circ}$',int2str(i),' - ECI-J2000 Frame');
    [t,y] = ode113(@TBP,[0 T],y0,opts,GM);
    y0 = y(end,:);
    h = plot(y(:,1),y(:,2),'k');


    hold on;

    %plot apocenter ellipses for all methods
    mu = mu_lin(1:2,1,j1);
    P = P_lin(1:2,1:2,j1);
    h1 = ellipse_visual(mu,P,95,'Color',[0.85,0.33,0.10]);
    mu = mu_ut(1:2,1,j1);
    P = P_ut(1:2,1:2,j1);
    h2 = ellipse_visual(mu,P,95,'Color',[0.93,0.69,0.13]);
    mu = mu_mc(1:2,1,j1);
    P = P_mc(1:2,1:2,j1);
    h3 = ellipse_visual(mu,P,95,'Color',[0.49,0.18,0.56]);

    hmc = plot(pop_mc(:,1,j1),pop_mc(:,2,j1),'r.');

    j2 = 2*i; %index of pericenter

    %plot pericenter ellipses for all methods
    mu = mu_lin(1:2,1,j2);
    P = P_lin(1:2,1:2,j2);
    ellipse_visual(mu,P,95,'Color',[0.85,0.33,0.10]);
    mu = mu_ut(1:2,1,j2);
    P = P_ut(1:2,1:2,j2);
    ellipse_visual(mu,P,95,'Color',[0.93,0.69,0.13]);
    mu = mu_mc(1:2,1,j2);
    P = P_mc(1:2,1:2,j2);
    ellipse_visual(mu,P,95,'Color',[0.49,0.18,0.56]);

    plot(pop_mc(:,1,j2),pop_mc(:,2,j2),'r.')
    axis equal; grid on;
    title(str_t,'Interpreter','latex')
    xlabel('x [km]','Interpreter','latex')
    ylabel('y [km]','Interpreter','latex')
    legend([h h1 h2 h3 hmc],'Orbit', 'First Order', 'Unscented Transform', 'Monte Carlo', 'MC Points','Interpreter','latex')
end

cspice_kclear();

%% End of File Functions

function f = TBP(~,s,mu)
%TBP two body problem rhs
%
% PROTOTYPE:
%     f = TBP(~,s,mu)
%
% INPUT:
%    t[dim]            time [s]
%    s[dim]            state [km; km/s]
%    mu[dim]           system gravitational constant [km^3/s^2]
%
% OUTPUT:
%    f[dim]           rhs of the system [-]
%
% CONTRIBUTORS:
%   Davide Zamblera - Politecnico di Milano - Space eng.
%
    r = s(1:3);
    v = s(4:6);

    rnrm = norm(r);

    f = zeros(6,1);
    f(1:3) = v;
    f(4:6) = -mu*r/rnrm^3;

end

function plotEarth()
%PLOTEARTH display a sphere the size of the Earth (radius in km)
%
% PROTOTYPE:
%     plotEarth()
%
% INPUT:
%
% OUTPUT:
%
% CONTRIBUTORS:
%   Davide Zamblera - Politecnico di Milano - Space eng.
%
    [X, Y, Z] = sphere;
    EarthRadius = cspice_bodvrd('EARTH','RADII',3);
    hSurface = surf(X*EarthRadius(1), Y*EarthRadius(1), Z*EarthRadius(1));
    set(hSurface,'FaceColor',[0.5 0.5 0.5])
end

function [position,isterminal,direction] = extremeEvents(t,y,mu)
  position = dot(y(1:3),y(4:6));
  isterminal = 0;  % Halt integration 
  direction = 0;   % The zero can be approached from either direction
end

function f = ode_stm(t,s,mu)
%ode_stm rhs of the stm variational equation, applied to keplerian dynamics
%
% PROTOTYPE:
%     f = ode_stm(t,s,mu)
%
% INPUT:
%    t[dim]            time [s]
%    s[dim]            state [km; km/s]
%    mu[dim]           system gravitational constant [km^3/s^2]
%
% OUTPUT:
%    f[dim]           rhs of the system [-]
%
% CONTRIBUTORS:
%   Davide Zamblera - Politecnico di Milano - Space eng.
%
    r = s(1:3);
    rx = r(1);
    ry = r(2);
    rz = r(3);
    v = s(4:6);
    PHI = s(7:42);
    PHI = reshape(PHI,[6 6]);
    rn = sqrt(rx^2+ry^2+rz^2);

    f = zeros(42,1);
    f(1:3) = v;
    f(4:6) = -mu/rn^3*r;
    A = zeros(6);
    A(1:3,4:6) = eye(3);
    A(4:6,1) = [(3*mu*rx^2)/rn^5 - mu/rn^3;
                (3*mu*rx*ry)/rn^5;
                (3*mu*rx*rz)/rn^5];

    A(4:6,2) = [(3*mu*rx*ry)/rn^5;
                (3*mu*ry^2)/rn^5 - mu/rn^3;
                (3*mu*ry*rz)/rn^5];

    A(4:6,3) = [(3*mu*rx*rz)/rn^5;
                (3*mu*ry*rz)/rn^5;
                (3*mu*rz^2)/rn^5 - mu/rn^3];

    M = A*PHI;
    f(7:42) = reshape(M,[36, 1]);
end

function [phi, x] = STM(x0,t0,t,opts)
%STM return state transition matrix computed with variational approach
%
% PROTOTYPE:
%     [phi, x] = STM(x0,t0,t,opts)
%
% INPUT:
%    x0[dim]           initial state [unit]
%    t0[dim]           initial time [unit]
%    t[dim]            final time [unit]
%    opts[dim]         options of the integrator scheme
%
% OUTPUT:
%    x[dim]            state at final time
%    phi[dim]          State Transition Matrix from t0 to t [unit]
%
% CONTRIBUTORS:
%   Davide Zamblera - Politecnico di Milano - Space eng.
%
    mu = cspice_bodvrd('EARTH','GM',1);
    x = [];
    if nargin < 4
        opts = [];
    end
    I = eye(6);
    x0 = [x0; reshape(I, [36, 1])];

    [~,s] = ode113(@ode_stm,[t0 t],x0,opts,mu);
    s = s(end,:)';
    phi = reshape(s(7:42),[6 6]);
    if nargout > 1
        x = s(1:6);
    end

end



function [mu, P] = linCov(mu0,P0,t0,t)
%LINCOV uncertainty propagation via LinCov approach
%
% INPUT:
%    mu0          mean value at epoch0
%    P0           covariance matrix at epoch0
%    t0           initial epoch
%    t            time at which to propagate
%
%
% OUTPUT:
%    mu           mean value at epoch t
%    P            covariance at epoch t
%
% CONTRIBUTORS:
%   Davide Zamblera - Politecnico di Milano - Space eng.
%
    opts = odeset('RelTol',1e-12,'AbsTol',1e-12);
    [phi, x] = STM(mu0,t0,t,opts);

    mu = x;
    P = phi*P0*phi';
    P = (P + P')/2;

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

function [mu, P] = unscentedTransform(mu0, P0, t0, t)
%LINCOV uncertainty propagation via Unscented Transform
%
% INPUT:
%    mu0          mean value at epoch0
%    P0           covariance matrix at epoch0
%    t0           initial epoch
%    t            time at which to propagate
%
%
% OUTPUT:
%    mu           mean value at epoch t
%    P            covariance at epoch t
%
% CONTRIBUTORS:
%   Davide Zamblera - Politecnico di Milano - Space eng.
%
    opts = odeset('RelTol',1e-12,'AbsTol',1e-12);
    GM = cspice_bodvrd('EARTH','GM',1);
    [x,wm,wc] = sigma_points(mu0,P0); %every sigma point is a column
    nsig = size(x,2);
    xprg = zeros(size(x));
    for i = 1:nsig
        [~,y] = ode113(@TBP,[t0 t],x(:,i),opts,GM);
        xprg(:,i) = y(end,:)';
    end

    %compute mean and covariance
    mu = zeros(size(mu0));
    P = zeros(size(P0));
    for i = 1:nsig
        mu = mu + wm(i)*xprg(:,i);
    end

    for i = 1:nsig
        P = P + wc(i)*(xprg(:,i) - mu)*(xprg(:,i) - mu)';
    end
end


function [mu, P, pop0, pop] = monteCarlo(mu0, P0, t0, t, n_sample)
%MONTECARLO uncertainty propagation via monte carlo simulation
%
% INPUT:
%    mu0          mean value at epoch0
%    P0           covariance matrix at epoch0
%    t0           initial epoch
%    t            time at which to propagate
%    n_sample     number of sample to propagate
%
%
% OUTPUT:
%    mu           mean value at epoch t
%    P            covariance at epoch t
%    pop0         states of the initial population of points
%    pop          states of the population at epoch t
%
% CONTRIBUTORS:
%   Davide Zamblera - Politecnico di Milano - Space eng.
%
    opts = odeset('RelTol',1e-12,'AbsTol',1e-12);
    GM = cspice_bodvrd('EARTH','GM',1);
    if nargin < 5
        n_sample = 100;
    end
    % generate initial population
    pop0 = mvnrnd(mu0,P0,n_sample);
    pop = zeros(size(pop0));
    for i = 1:n_sample
        [~,y] = ode113(@TBP,[t0 t],pop0(i,:),opts,GM);
        pop(i,:) = y(end,:);
    end

    mu = mean(pop)';
    P = cov(pop);

end


function h = ellipse_visual(mu,P,conf,varargin)
% Display a representation of a confidence ellipse for gaussian distributions
%
% INPUT:
%    mu           mean value
%    P            covariance matrix
%    conv         confidence interval as percentage
%
%
% OUTPUT:
%	 none
%
% CONTRIBUTORS:
%   Davide Zamblera - Politecnico di Milano - Space eng.
%

    %generate base circle

    alfa = 0:0.01:2*pi;
    points = [cos(alfa); sin(alfa)];


    %extract statistic data
    [V,D] = eigs(P);

    if nargin < 3
        conf = 95;
    end
    vec_conf = [99.5 99 97.5 95 90];
    K = [10.597   9.21   7.378   5.991   4.605];
    K = K(vec_conf == conf);
    if isempty(K)
        error('Selected confidence interval unknown.')
    end
    a = sqrt(K*D(1,1));
    b = sqrt(K*D(2,2));
    Mstretch = [a 0; 0 b];
    v = V(:,1); %eigenvector related to highest eigenvalue
    theta = atan2(v(2),v(1));
    Mrotate = [cos(theta), -sin(theta); sin(theta) cos(theta)];
    points = Mrotate*Mstretch*points;
    h = plot(points(1,:)+mu(1),points(2,:)+mu(2),varargin{:});
end