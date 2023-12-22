% Spacecraft Guidance and Navigation (2020/2021)
% Assignment # 1
% Author: Davide Zamblera

%% Section 1: Understand dynamics and Order of Magnitudes of input
clearvars; close all; clc;

%load kernels
cspice_furnsh('metakernel1.tm');

%data initialization and change of units
%%% SET OF UNITS [AU,YEARS,KG]
date_lnc = '2022-08-03-12:45:20.000 UTC';
ti = cspice_str2et(date_lnc);
xert = units(cspice_spkezr('EARTH',ti,'ECLIPJ2000','NONE','SUN'));

k = 1e-3; % scale to take [m]-->[km]
Nthr = input('Choose the number of thrusters: ');
const.ti = ti;
const.m0 = 1500; %Kg
const.Tmax = units(k*Nthr*150e-3,[1 -2 1]); %N
const.Isp = units(3000,[0 1 0]); %s
const.mu = units(cspice_bodvrd('SUN','GM',1),[3 -2 0]);
const.g0 = units(k*9.80665,[1 -2 0]);
const.x0 = xert;


date_frs = '01-Aug-2022';
date_end = '01-Dec-2024';

% Validate Azimuth and Elevation

% azimuth([-1 0 0])
% azimuth([0.5 0.5 0])
% elevation([0.5 0.5 1])
% elevation([0.5 0.5 -1])

% Begin an infinite cycle to explore dynamics by user input
while 1

l0 = input('Enter vector of lambda0: ');
if strcmp(l0,'quit')
    break;
end
if isrow(l0)
    l0 = l0';
end
date_arr = input('Enter arrival date: ','s');
tf = cspice_str2et(date_arr);
lm0 = input('Enter value of lambdaM: ');
init_state = [const.x0; const.m0; l0; lm0];

[~,y] = plot_results(l0,lm0,tf,const);

hml(init_state,const)
hml(y(end,:)',const)
pause;
close all;
end

%clear kernels
cspice_kclear;

%% Section 2: Search for case with 4 thruster

clearvars; close all; clc;

%load kernels
cspice_furnsh('metakernel1.tm');

%data initialization and change of units
%%% SET OF UNITS [AU,YEARS,KG]
date_lnc = '2022-08-03-12:45:20.000 UTC';
options = optimoptions(@fsolve,'Display','iter');
ti = cspice_str2et(date_lnc);
xert = units(cspice_spkezr('EARTH',ti,'ECLIPJ2000','NONE','SUN'));

k = 1e-3; % scale to take [m]-->[km]
Nthr = 4;
const.ti = ti;
const.m0 = 1500; %Kg
const.Tmax = units(k*Nthr*150e-3,[1 -2 1]); %N
const.Isp = units(3000,[0 1 0]); %s
const.mu = units(cspice_bodvrd('SUN','GM',1),[3 -2 0]);
const.g0 = units(k*9.80665,[1 -2 0]);
const.x0 = xert;


date_arr_1 = '2022-02-04-12:45:20.000 UTC';
date_arr_2 = '2023-06-04-12:45:20.000 UTC';
tf_1 = cspice_str2et(date_arr_1);
tf_2 = cspice_str2et(date_arr_2);


%Iteratively construct random guess and solve shooting problem
rng(0,'twister');

l0 = zeros(6,1);


while 1

l0(1) = randrange(-10,10)*1e-3;
l0(2) = randrange(-10,10)*1e-3;
l0(3) = randrange(0,0.05);
l0(4) = randrange(-10,-10)*1e-3;
l0(5) = randrange(-10,-10)*1e-3;
l0(6) = randrange(0,0.05)*1e-10;
tf = tf_1 + (tf_2 - tf_1)*rand;
lm0 = randrange(0,0.05);



z0 = [l0; lm0; tf];
[z,fval,exitflag,output] = fsolve(@(z) low_thrust_zero(z,const),z0,options);


if exitflag > 0
    fprintf('****************** Solution found ********************\n');
    l0s = z(1:6);
    lm0s = z(7);
    tfs = z(8);
    tfs_str = cspice_timout(tfs,'YYYY-MON-DD-HR:MN:SC.####::UTC');
    fprintf('l0 = [');
    fprintf('%g ', l0s);
    fprintf(']\n');
    fprintf('lm0 = %f\n',lm0s);
    fprintf('tf = %s\n',tfs_str);
    init_state = [const.x0; const.m0; l0s; lm0s];
    [~,y] = plot_results(l0s,lm0s,tfs,const);
    H1 = hml(init_state,const);
    H2 = hml(y(end,:)',const);
    fprintf('The Hamiltonian at the initial state is: %f\n',H1)
    fprintf('The Hamiltonian at the final state is: %f\n',H2)
    break;
end
end

%clear kernels
cspice_kclear;

%%%%%%%%%% SOLUTION %%%%%%%%%%%%%

% l0 = [-4.04035632136414	2.91114172227394	0.200557648631166	-0.826551898588282	-0.385114835106420	-0.0106874342906355 ]
% lm0 = 0.000476969004661305
% tf = 2023-05-12-08:42:14.726 UTC


%% Section 3: Search for case with 3 thruster

clearvars; close all; clc;

%load kernels
cspice_furnsh('metakernel1.tm');


%data initialization and change of units
%%% SET OF UNITS [AU,YEARS,KG]
date_lnc = '2022-08-03-12:45:20.000 UTC';
options = optimoptions(@fsolve,'Display','iter');
ti = cspice_str2et(date_lnc);
xert = units(cspice_spkezr('EARTH',ti,'ECLIPJ2000','NONE','SUN'));

k = 1e-3; % scale to take [m]-->[km]
Nthr = 3;
const.ti = ti;
const.m0 = 1500; %Kg
const.Tmax = units(k*Nthr*150e-3,[1 -2 1]); %N
const.Isp = units(3000,[0 1 0]); %s
const.mu = units(cspice_bodvrd('SUN','GM',1),[3 -2 0]);
const.g0 = units(k*9.80665,[1 -2 0]);
const.x0 = xert;


date_arr_1 = '2022-02-04-12:45:20.000 UTC';
date_arr_2 = '2023-06-04-12:45:20.000 UTC';
tf_1 = cspice_str2et(date_arr_1);
tf_2 = cspice_str2et(date_arr_2);

%Iteratively construct random guess and solve shooting problem
rng(0,'twister');

l0 = zeros(6,1);


while 1

l0(1) = randrange(-10,10)*1e-3;
l0(2) = randrange(-10,10)*1e-3;
l0(3) = randrange(0,0.05);
l0(4) = randrange(-10,-10)*1e-3;
l0(5) = randrange(-10,-10)*1e-3;
l0(6) = randrange(0,0.05)*1e-10;
tf = tf_1 + (tf_2 - tf_1)*rand;
lm0 = randrange(0,0.05);


%init_state = [x0; const.m0; l0; lm0];


z0 = [l0; lm0; tf];
[z,fval,exitflag,output] = fsolve(@(z) low_thrust_zero(z,const),z0,options);


if exitflag > 0
    fprintf('****************** Solution found ********************\n');
    l0s = z(1:6);
    lm0s = z(7);
    tfs = z(8);
    tfs_str = cspice_timout(tfs,'YYYY-MON-DD-HR:MN:SC.####::UTC');
    fprintf('l0 = [');
    fprintf('%g ', l0s);
    fprintf(']\n');
    fprintf('lm0 = %f\n',lm0s);
    fprintf('tf = %s\n',tfs_str);
    [~,y] = plot_results(l0s,lm0s,tfs,const);
    break;
end
end

%clear kernels
cspice_kclear;

%%%%%%%%%% SOLUTION %%%%%%%%%%%%%

% l0 = [-11.4406359406688	9.46747474260996	0.156422799663520	-2.37511149621976	-1.25128436452581	0.0182548021769254 ]
% lm0 = 0.00118609139151830
% tf = 2023-10-16-03:55:14.519 UTC

%% End of File Functions

function f = TBP(~,x,mu)
%TBP ode rhs for the restricted two body problem
%
% PROTOTYPE:
%     f = TBP(~,x,mu)
%
% INPUT:
%    t[1]            indipendent variable time [s]
%    x[6x1]          state x=[r,v] [km, km/s]
%    mu[1]           var description [km^3/s^2]
%
% OUTPUT:
%    f[6x1]          right hand side function [-]
%
% CONTRIBUTORS:
%   Davide Zamblera - Politecnico di Milano - Space eng.
%
 
r = x(1:3);
v = x(4:6);
r_3 = norm(r)^3;

f(1:3,1) = v;
f(4:6,1) = -mu*r/r_3;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = ode_tbp_agm(~,y,const)
%ODE_TBP_AGM ode of the 2 body problem augmented with costate and
%continuous thrust
%
% PROTOTYPE:
%     f = ode_tbp_agm(~,y,const)
%
% INPUT:
%    t[1]            indipendent variable time [s]
%    y[6x1]          state y=[r,v,lambdar,lambdav] [km, km/s, - , - ]
%    const[str]      struct that contains all constants needed by the model
%
% OUTPUT:
%    f[6x1]          right hand side function [-]
%
% CONTRIBUTORS:
%   Davide Zamblera - Politecnico di Milano - Space eng.
%

    mu = const.mu;
    Tmax = const.Tmax;
    Isp = const.Isp;
    g0 = const.g0;
    r = y(1:3);
    rnrm = norm(r);
    v = y(4:6);
    m = y(7);
    lr = y(8:10);
    lv = y(11:13);
    lm = y(14);
    lvnrm = norm(lv);
    
    swtc = (-lvnrm*Isp*g0/m - lm);
    if swtc > 0
        u = 0;
    elseif swtc < 0
        u = 1;
    else
        warning('ode_time_opt->ulaw: Switch value is zero, behaviour is undetermined');
        u = inf;
    end
    

    f = zeros(14,1);
    f(1:3) = v;
    
    f(4:6) = -mu/rnrm^3*r - u*Tmax/m*lv/lvnrm;
    f(7) = -u*Tmax/(Isp*g0);
    f(8:10) = -3*mu/rnrm^5*(r'*lv)*r + mu/rnrm^3*lv;
    f(11:13) = -lr;
    f(14) = -u*lvnrm*Tmax/m^2;
end


function y = low_thrust_prg(x0,tspan,const,opts)
%LOW_THRUST_PRG propagator of state with TBPAGM model
%
% PROTOTYPE:
%     y = low_thrust_prg(x0,t0,t,const,opts)
%
% INPUT:
%    x0[dim]           initial state [unit]
%    tspan[2]          integration time span [unit]
%    const[str]        struct that contains all constants needed by the model
%    opts[dim]         options of the integrator scheme
%
% OUTPUT:
%    y[dim]             state at final time [unit]
%
% CONTRIBUTORS:
%   Davide Zamblera - Politecnico di Milano - Space eng.
%
    if nargin < 5
        opts = odeset('RelTol',1e-12,'AbsTol',1e-12);
    end

    [~,y] = ode113(@(t,y) ode_tbp_agm(t,y,const),tspan,x0,opts);
    y = y(end,:)';
end

function [tt,y] = low_thrust_traj(x0,tspan,const,opts)
%LOW_THRUST_TRAJ return trajectory by integration of x0 with TBP_AGM model
%
% PROTOTYPE:
%     [tt,y] = low_thrust_traj(x0,tspan,const,opts)
%
% INPUT:
%    x0[dim]           initial state [unit]
%    tspan[2]          integration time span [unit]
%    const[str]        struct that contains all constants needed by the model
%    opts[dim]         options of the integrator scheme
%
% OUTPUT:
%    tt[dim]            times of integration
%    y[dim]             state at the times of integration [unit]
%
% CONTRIBUTORS:
%   Davide Zamblera - Politecnico di Milano - Space eng.
%

    if nargin < 5
        opts = odeset('RelTol',1e-12,'AbsTol',1e-12);
    end

    [tt,y] = ode113(@(t,y) ode_tbp_agm(t,y,const),tspan,x0,opts);
end



function fcon = low_thrust_zero(z0,const)
%LOW_THRUST_ZERO nonlinear equations to solve optimal control
%
% PROTOTYPE:
%     fcon = low_thrust_zero(z0,const)
%
% INPUT:
%    z0[dim]           state of shooting problem z0 = [l0, lm0, tf] [unit]
%    const[str]        struct that contains all constants needed by the model
%
% OUTPUT:
%    fcon[dim]         nonlinear equations of the shooting problem
%
% CONTRIBUTORS:
%   Davide Zamblera - Politecnico di Milano - Space eng.
%
    mu = const.mu;

    l0 = z0(1:6);
    lm0 = z0(7);
    tf = z0(8);
    ti = const.ti;
    x0 = const.x0;
    m0 = const.m0;
    w0 = [x0; m0; l0; lm0];
    tspan = [units(ti,[0 1 0]), units(tf,[0 1 0])];
    wf = low_thrust_prg(w0,tspan,const);
    fcon = zeros(8,1);
    xf = wf(1:6);
    lf = wf(8:13);
    lmf = wf(14);
    xmars = units(cspice_spkezr('MARS',tf,'ECLIPJ2000','NONE','SUN'));
    fcon(1:6) = (xf - xmars);
    fcon(7) = lmf;
    fMars = TBP(tf,xmars,mu);
    vmars = xmars(4:6);
    amars = fMars(4:6);
    H = hml(wf,const);
    fcon(8) = (H - lf(1:3)'*vmars - lf(4:6)'*amars);

    
end

function H = hml(w,const)
%HML compute hamiltonian given the state
%
% PROTOTYPE:
%     H = hml(w,const)
%
% INPUT:
%    w[dim]           state w=[r,v,m,lambdar,lambdav,lambdam] [km, km/s,kg , - , - , - ]
%    const[str]       struct that contains all constants needed by the model
%
% OUTPUT:
%    H[dim]         hamiltonian value
%
% CONTRIBUTORS:
%   Davide Zamblera - Politecnico di Milano - Space eng.
%
    mu = const.mu;
    Tmax = const.Tmax;
    Isp = const.Isp;
    g0 = const.g0;

    r = w(1:3);
    rnrm = norm(r);
    v = w(4:6);
    m = w(7);
    lr = w(8:10);
    lv = w(11:13);
    lvnrm = norm(lv);
    lm = w(14);

    swtc = (-lvnrm*Isp*g0/m - lm);
    if swtc > 0
        u = 0;
    elseif swtc < 0
        u = 1;
    else
        error('ode_time_opt->ulaw: Switch value is zero, behaviour is undetermined');
    end
    
    H = 1 + lr'*v - mu/rnrm^3*r'*lv + Tmax/(Isp*g0)*u*swtc;

end


function tdate = xaxis_date(tt)
%XAXIS_DATE generate date axis values
%
% PROTOTYPE:
%     tdate = xaxis_date(tt)
%
% INPUT:
%    tt[dim]          time values expressed in TDB seconds
%
% OUTPUT:
%    tdate[dim]       time in datetime data type
%
% CONTRIBUTORS:
%   Davide Zamblera - Politecnico di Milano - Space eng.
%
    pictor = 'DD-Mon-YYYY HR:MN:SC';
    tstr = cspice_timout(tt',pictor);
    tstr = cellstr(tstr);
    tdate = datetime(tstr);

end


function h = elevation(v)
%ELEVATION compute elevation angle from state
%
% PROTOTYPE:
%     h = elevation(v)
%
% INPUT:
%    v[dim]       position state
%
% OUTPUT:
%    h[dim]       elevation angle [rad]
%
% CONTRIBUTORS:
%   Davide Zamblera - Politecnico di Milano - Space eng.
%
    v1 = [v(1:2) 0];
    v2 = v;
    x = cross(v1,v2);

    n = x/norm(x); %n cross v1, show + angles w.r.t v1
    if v(3) < 0
        n = -n;
    end
    c = sign(dot(x,n)) * norm(x);
    h = atan2(c,dot(v1,v2));
end

function az = azimuth(v)
%AZIMUTH compute azimuth angle from state, positive clockwise from x
%
% PROTOTYPE:
%     az = azimuth(v)
%
% INPUT:
%    v[dim]       position state
%
% OUTPUT:
%    az[dim]      azimuth angle [rad]
%
% CONTRIBUTORS:
%   Davide Zamblera - Politecnico di Milano - Space eng.
%
    v1 = [1 0 0];
    v2 = v;
    x = cross(v1,v2);
    n = [0 0 -1]; %n cross v1, show + angles w.r.t v1
    c = sign(dot(x,n)) * norm(x);
    az = atan2(c,dot(v1,v2));
end

function alfa = compute_angles(s)
%COMPUTE_ANGLES compute thrust angles from state
%
% PROTOTYPE:
%     az = azimuth(v)
%
% INPUT:
%    s[dim]       state s = [r, v, m, lambdar, lambdav, lambdam]
%
% OUTPUT:
%    alfa[dim]    vector of angles alfa = [az, h] (deg)
%
% CONTRIBUTORS:
%   Davide Zamblera - Politecnico di Milano - Space eng.
%
    lv = -s(:,11:13); %minus because of primer vector definition 
    n = size(lv,1);
    alfa = zeros(n,2);
    for i = 1:n
        alfa(i,1) = azimuth(lv(i,:));
        alfa(i,2) = elevation(lv(i,:));
    end
    alfa = unwrap(alfa,[],1); %avoid jumps
    alfa = rad2deg(alfa);
end

function u = thrust_level(s, const)
%THRUST_LEVEL compute thrust levels from state via switch function
%
% PROTOTYPE:
%     u = thrust_level(s, const)
%
% INPUT:
%    s[dim]       state s = [r, v, m, lambdar, lambdav, lambdam]
%    const[str]       struct that contains all constants needed by the model
%
% OUTPUT:
%    u[dim]    thrust level
%
% CONTRIBUTORS:
%   Davide Zamblera - Politecnico di Milano - Space eng.
%
    Isp = const.Isp;
    g0 = const.g0;

    lv = s(:,11:13);
    lm = s(:,14);
    m = s(:,7);

    n = length(m);
    u = zeros(n,1);
    for i = 1:n
        lvnrm = norm(lv(i,:));
        swtc = (-lvnrm*Isp*g0/m - lm);
        if swtc > 0
            u(i) = 0;
        else
            u(i) = 1;
        end
    end

end


function xn = units(x,n)
%UNITS convert units to the triplet [AU,years,kg]
%
% PROTOTYPE:
%     xn = units(x,n)
%
% INPUT:
%    x[dim]       variable x, if it is the state x=[r,v] do not provide n
%    n[dim]       dimensions of x [L,T,M]
%
% OUTPUT:
%    xn[dim]    variable x in the triplet [AU,years,kg]
%
% CONTRIBUTORS:
%   Davide Zamblera - Politecnico di Milano - Space eng.
%
    cnv_l = cspice_convrt(1,'KM','AU');
    cnv_t = cspice_convrt(1,'SECONDS','YEARS');
    cnv_v = cnv_l/cnv_t;
    cnv_m = 1;
    if nargin == 2
        xn = x*cnv_l^n(1)*cnv_t^n(2)*cnv_m^n(3);
    else
        cnv = [cnv_l, cnv_l, cnv_l, cnv_v, cnv_v, cnv_v]';
        xn = cnv.*x;
    end

end

function [tt,y] = plot_results(l0,lm0,tf,const)
%PLOT_RESULTS given the initial guess [l0,lm0,tf] integrate and plot results
%
% PROTOTYPE:
%     y = plot_results(l0,lm0,tf,const)
%
% INPUT:
%    l0[dim]       initial costate [lambdar, lambdav]
%    lm0[dim]      initial costate lambdam
%    tf[dim]       arrival time
%    const[str]       struct that contains all constants needed by the model
%
% OUTPUT:
%    tt[dim]       integration times
%    y[dim]        integrated state
%
% CONTRIBUTORS:
%   Davide Zamblera - Politecnico di Milano - Space eng.
%

init_state = [const.x0; const.m0; l0; lm0];
tspan = [units(const.ti,[0 1 0]), units(tf,[0 1 0])];
[tt,y] = low_thrust_traj(init_state,tspan,const);
tt = units(tt,[0 -1 0]);
xmrs = units(cspice_spkezr('MARS',tt','ECLIPJ2000','NONE','SUN'));
xmrs = xmrs';
tt = xaxis_date(tt);
err_pos = zeros(length(tt),3);
err_pos(:,1) = xmrs(:,1)-y(:,1);
err_pos(:,2) = xmrs(:,2)-y(:,2);
err_pos(:,3) = xmrs(:,3)-y(:,3);
err_pos_nrm = vecnorm(err_pos,2,2);


date_frs = '01-Aug-2022';
date_end = '01-Dec-2024';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot of the position errors
subplot(2,2,1);
plot(tt,err_pos(:,1),tt,err_pos(:,2),tt,err_pos(:,3),tt,err_pos_nrm)
grid on;
xlabel('Time [UTC]','Interpreter','latex');
ylabel('Position Errors [AU]','Interpreter','latex');
legend('$err_{x}$','$err_{y}$','$err_{z}$','$||err_{norm}||$','Interpreter','latex')
title('Positions','Interpreter','latex');
xtickformat('dd-MMM-yyyy')
xticks(datetime(date_frs) : calmonths(3) : datetime(date_end))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot of the trajectories
subplot(2,2,2);
plot3(y(:,1),y(:,2),y(:,3))
grid on;
axis equal
hold on;
plot3(xmrs(:,1),xmrs(:,2),xmrs(:,3))
plot3(0,0,0,'o','MarkerFaceColor','yellow','MarkerSize',7);
xert = const.x0;
plot3(xert(1),xert(2),xert(3),'o','MarkerFaceColor','blue','MarkerSize',4);
plot3(xmrs(end,1),xmrs(end,2),xmrs(end,3),'o','MarkerFaceColor','red','MarkerSize',3);
xlabel('x [AU]','Interpreter','latex');
ylabel('y [AU]','Interpreter','latex');
zlabel('z [AU]','Interpreter','latex');
title('ECLIPJ2000 Frame @Sun center','Interpreter','latex');
legend({'Spacecraft Trajectory','Mars Orbit','','',''},'Interpreter','latex');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot of the errors in velocity and difference in norm
subplot(2,2,3);
err_vel = zeros(length(tt),3);
err_vel(:,1) = xmrs(:,4)-y(:,4);
err_vel(:,2) = xmrs(:,5)-y(:,5);
err_vel(:,3) = xmrs(:,6)-y(:,6);
err_vel_nrm = vecnorm(y(:,4:6),2,2) - vecnorm(xmrs(:,4:6),2,2);
k = cspice_convrt(1,'AU','KM')/cspice_convrt(1,'YEARS','SECONDS');
err_vel = k*err_vel;
err_vel_nrm = k*err_vel_nrm;
plot(tt,err_vel(:,1),tt,err_vel(:,2),tt,err_vel(:,3),tt,err_vel_nrm)
grid on;
xlabel('Time [UTC]','Interpreter','latex');
ylabel('Velocity Errors [km/s]','Interpreter','latex');
legend('$\dot{x}_{err}$','$\dot{y}_{err}$','$\dot{z}_{err}$','$||v_{s/c}|| - ||v_{mars}||$','Interpreter','latex')
title('Velocities','Interpreter','latex');
xtickformat('dd-MMM-yyyy')
xticks(datetime(date_frs) : calmonths(3) : datetime(date_end))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot of mass consumption
subplot(2,2,4)
plot(tt,y(:,7))
grid on;
xlabel('Time [UTC]','Interpreter','latex');
ylabel('Mass [kg]','Interpreter','latex');
title('Mass consumption','Interpreter','latex');
xtickformat('dd-MMM-yyyy')
xticks(datetime(date_frs) : calmonths(3) : datetime(date_end))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot of lambdam
figure;
subplot(1,3,1);
plot(tt,y(:,14))
grid on;
xlabel('Time [UTC]','Interpreter','latex');
ylabel('$\lambda_m [-]$','Interpreter','latex');
xtickformat('dd-MMM-yyyy')
xticks(datetime(date_frs) : calmonths(3) : datetime(date_end))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot of the thrust angles
theta = compute_angles(y);
subplot(1,3,2);
plot(tt,theta)
grid on;
xlabel('Time [UTC]','Interpreter','latex');
ylabel('Angles [deg]','Interpreter','latex');
title('Thrust Angles [deg]','Interpreter','latex');
legend({'Azimuth','Elevation'},'Interpreter','latex')
xtickformat('dd-MMM-yyyy')
xticks(datetime(date_frs) : calmonths(3) : datetime(date_end))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot of the thrust level
u = thrust_level(y, const);
subplot(1,3,3);
plot(tt,u)
grid on;
xlabel('Time [UTC]','Interpreter','latex');
ylabel('u [-]','Interpreter','latex');
title('Thrust level','Interpreter','latex');
xtickformat('dd-MMM-yyyy')
xticks(datetime(date_frs) : calmonths(3) : datetime(date_end))

err_pos_fin = units(err_pos_nrm(end),[-1 0 0]);
fprintf('Final error in position norm [km]: %.11f\n',err_pos_fin);
k = cspice_convrt(1,'AU','M')/cspice_convrt(1,'YEARS','SECONDS'); % [AU/year] -> [m/s]
err_vel_fin = k*err_vel_nrm(end);
fprintf('Final error in velocity norm [m/s]: %.11f\n',err_vel_fin);
tof = (tf - const.ti)*cspice_convrt(1,'SECONDS','DAYS');
fprintf('Transfer time [days]: %f\n',tof);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot of the hamiltonian
hml_v = zeros(length(tt),1);
for j = 1:length(tt)
    hml_v(j) = hml(y(j,:)',const);

end
figure;
plot(tt,hml_v);
grid on;
xlabel('Time [UTC]','Interpreter','latex');
ylabel('Hamiltonian [-]','Interpreter','latex');
xtickformat('dd-MMM-yyyy')
title('Hamiltonian Evolution in time','Interpreter','latex');
xticks(datetime(date_frs) : calmonths(3) : datetime(date_end))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot of the thrust direction in the x-y plane
lambdav = y(:,11:13);
lvnrm = vecnorm(lambdav,2,2);
alfa = -lambdav./lvnrm;
k = 0.5; %scale factor
alfapl = k*alfa(:,1:2);

figure;
plot(y(:,1),y(:,2),'b');
grid on;
hold on;

plot(xmrs(:,1),xmrs(:,2),'r');

plot3(0,0,0,'o','MarkerFaceColor','yellow','MarkerSize',7);
xert = const.x0;
plot3(xert(1),xert(2),xert(3),'o','MarkerFaceColor','blue','MarkerSize',4);
plot3(xmrs(end,1),xmrs(end,2),xmrs(end,3),'o','MarkerFaceColor','red','MarkerSize',3);
xlabel('x [AU]','Interpreter','latex');
ylabel('y [AU]','Interpreter','latex');
zlabel('z [AU]','Interpreter','latex');
title('ECLIPJ2000 Frame @Sun center','Interpreter','latex');

spc = 5; %space between arrows
h = quiver(y(1:spc:end,1),y(1:spc:end,2),alfapl(1:spc:end,1),alfapl(1:spc:end,2),'k');
set(h,'MaxHeadSize',0.05,'AutoScaleFactor',1);
legend({'Spacecraft Trajectory','Mars Orbit','','','','Thrust Directions'},'Interpreter','latex');
axis equal;

end

function n = randrange(lv,uv)
    n = lv + (uv - lv)*rand;
end