% Spacecraft Guidance and Navigation (2020/2021)
% Assignment # 2
% Author: Davide Zamblera

%% Section 1
clearvars; close all; clc;
set(0, 'defaulttextinterpreter','latex');  
set(0, 'defaultAxesTickLabelInterpreter','latex');  
set(0, 'defaultLegendInterpreter','latex');

addpath('sgp4');
cspice_furnsh('assignment02.tm');
rng('default');
%data definition
t_utc = '2022-11-11T19:08:49.824';
et_init = cspice_str2et(t_utc);

r0 = [6054.30795817484, -3072.03883303992, -133.115352431876]';
v0 = [4.64750094824087, 9.18608475681236, -0.62056520749034]';
y0 = [r0; v0];
mu0 = y0;

GM = cspice_bodvrd('EARTH','GM',1);
str_t0 = '2022-11-12T04:30:00.000';
str_tf = '2022-11-14T16:30:00.000';
et_t0 = cspice_str2et(str_t0);
et_tf = cspice_str2et(str_tf);

opts = odeset('RelTol',1e-12,'AbsTol',1e-12);
t = et_t0:60:et_tf;
[~,x_eci] = TBP_traj(y0,[et_init t],GM,opts);
x_eci = x_eci(2:end,:);

figure;
plotEarth()
hold on;
plot3(x_eci(:,1),x_eci(:,2),x_eci(:,3))
axis equal;
x_eci = x_eci';
R=42241.08;
alfa = 0:0.01:3*pi;
x = R*cos(alfa);
y = R*sin(alfa);
z = zeros(size(y));
plot3(x,y,z)


%compute visibility windows

spacecraftName = 'Ariane 5 US';
stationNames = {'KOUROU','PERTH'};
stationMinEl = [10 5];
n_stat = length(stationNames);
i_visibility = cell(n_stat,1);
windows = cell(n_stat,1);
t_begin_window = cell(n_stat,1);
t_end_window = cell(n_stat,1);
h_min = [10, 5]; %minimum visible elevation
f_handle = figure;
styles = ['r.';'b.'];
for i = 1:n_stat
    [Az, h, rho] = antenna_pointing(stationNames{i}, t, x_eci);
    figure
    plot(t/cspice_spd(), unwrap(Az)*cspice_dpr(), 'DisplayName', spacecraftName)
    grid on;
    title(['Topocentric Frame @',stationNames{i}(1),lower(stationNames{i}(2:end))])
    xlabel('Epoch [MJD2000]')
    ylabel('Azimuth [deg]')
    legend(spacecraftName)
    
    figure
    plot(t/cspice_spd(), h*cspice_dpr(),'DisplayName', spacecraftName)
    grid on;
    title(['Topocentric Frame @',stationNames{i}(1),lower(stationNames{i}(2:end))])
    xlabel('Epoch [MJD2000]')
    ylabel('Elevation [deg]')
    legend(spacecraftName)
    ylim([stationMinEl(i) max(ylim)])
    
    % Plot passes (azimuth and elevation)
    figure(f_handle)
    hold on;
    i_visibility{i} = h > stationMinEl(i)*cspice_rpd();
    windows{i} = t(i_visibility{i});
    plot(Az(i_visibility{i})*cspice_dpr(), h(i_visibility{i})*cspice_dpr(),styles(i,:),'DisplayName', spacecraftName)
    grid on;
    axis([-180,180,0, 90])
    title('Visibility passes')
    xlabel('Azimuth [deg]')
    ylabel('Elevation [deg]')

    [t_begin_window{i},t_end_window{i}] = compute_visibility(h*cspice_dpr,t,h_min(i));
end
figure(f_handle)
legend({'Kourou passes', 'Perth passes'});
% Simulate measurements

tlestr1 = '1 87654U 22110B   22316.00967942  .00000002  00000-0  32024-3 0  9990';
tlestr2 = '2 87654   3.6309 137.1541 8138191 196.1632  96.6141  1.26411866   834';
sigma_angles = 100e-3*cspice_rpd; %[rad]
sigma_range = 0.01; %[km]
mu = zeros(3,1);
R = diag([sigma_angles^2 sigma_angles^2 sigma_range^2]);

measure = cell(n_stat,1);
measure_corrupt = cell(n_stat,1);
for i = 1:n_stat
    [reci, veci] = tle2eci_state(tlestr1,tlestr2, windows{i});
    x_eci = [reci; veci];
    [Az, h, rho] = antenna_pointing(stationNames{i}, windows{i}, x_eci);
    measure{i} = [Az; h; rho];
    noise = mvnrnd(mu,R,length(Az))';
    measure_corrupt{i} = measure{i} + noise;
end


% 3) Batch filters


% 3a) batch filter: solve least squares

tref = et_t0;
sigma_meas = [sigma_angles sigma_angles sigma_range];
[r_ref, v_ref] = tle2eci_state(tlestr1,tlestr2, tref);
xref = [r_ref; v_ref];

% Plot reference orbit
figure;
plotEarth()
hold on;
[~,xv] = ode113(@TBP,[et_t0 et_tf],xref,opts,GM);
plot3(xv(:,1),xv(:,2),xv(:,3))
axis equal;


opt = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt', 'Display', 'iter');
t = windows{2};
meas_real = measure_corrupt{2};
fun = @(xest) BatchPerthKep(xest, t, tref, meas_real, sigma_meas);
[xest, resnorm, residual, exitflag, ~, ~, jac] = lsqnonlin(fun, xref, [], [], opt);
% figures here


% plot perth estimation
Pl1 = post_process_batch(xest, residual, exitflag, t, meas_real, sigma_meas, resnorm, jac)

% 3b) batch filter: use all measurements

tref = et_t0;
sigma_meas = [sigma_angles sigma_angles sigma_range];
[r_ref, v_ref] = tle2eci_state(tlestr1,tlestr2, tref);

opt = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt', 'Display', 'iter');
xref = [r_ref; v_ref];
fun = @(xest) BatchAllKep(xest, windows, tref, measure_corrupt, sigma_meas);
[xest, resnorm, residual, exitflag, ~, ~, jac] = lsqnonlin(fun, xref, [], [], opt);

% plot kourou+perth estimation
Pl2 = post_process_batch(xest, residual, exitflag, windows, measure_corrupt, sigma_meas, resnorm, jac)

% 3c) batch filter: use all measurements and J2 model

tref = et_t0;
sigma_meas = [sigma_angles sigma_angles sigma_range];
[r_ref, v_ref] = tle2eci_state(tlestr1,tlestr2, tref);

opt = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt', 'Display', 'iter');
t = windows{2};
xref = [r_ref; v_ref];
fun = @(xest) BatchAllJ2(xest, windows, tref, measure_corrupt, sigma_meas);
[xest, resnorm, residual, exitflag, ~, ~, jac] = lsqnonlin(fun, xref, [], [], opt);

% plot perth estimation
Pl3 = post_process_batch(xest, residual, exitflag, windows, measure_corrupt, sigma_meas, resnorm, jac)


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

function f = TBP_J2(t,s,mu)
%TBP_J2 two body problem with J2 perturbation rhs
%
% PROTOTYPE:
%     f = TBP_J2(t,s,mu)
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
    Re = cspice_bodvrd('EARTH','RADII',3);
    Re = Re(1);
    J2 =  0.0010826269;

    r = s(1:3);
    v = s(4:6);
    rnrm = norm(r);

    f = zeros(6,1);
    f(1:3) = v;
    f(4:6) = -mu*r/rnrm^3;

    rotm = cspice_pxform('J2000','ITRF93',t);
    r_ecef = rotm*r;
    aJ2_ecef = 3/2*mu*J2*r_ecef/rnrm^3*(Re/rnrm)^2.*(5*(r_ecef(3)/rnrm)^2 - [1; 1; 3]);
    aJ2 = rotm'*aJ2_ecef;
    f(4:6) = f(4:6) + aJ2;
end


function [tt,y] = TBP_traj(x0,tspan,mu,opts)
%TBP_traj compute the trajectory of the spacecraft, two body problem
%
% PROTOTYPE:
%     [tt,y] = TBP_traj(x0,tspan,mu,opts)
%
% INPUT:
%    x0[dim]           initial state
%    tspan[dim]        interval of integration [s]
%    mu[dim]           system gravitational constant [km^3/s^2]
%    opts[struct]      options of integration [-]
%
% OUTPUT:
%    tt[dim]           time [s]
%    y[dim]            state for every time of integration [-]
%
% CONTRIBUTORS:
%   Davide Zamblera - Politecnico di Milano - Space eng.
%
    if nargin < 4
        opts = [];
    end
    [tt,y] = ode113(@TBP,tspan,x0,opts,mu);
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

function [Az, h, rho, rho_dot] = antenna_pointing(stationName, t, x_eci)
%ANTENNA_POINTING from a state in eci compute measurements of radar antenna
%
% PROTOTYPE:
%     [Az, h, rho, rho_dot] = antenna_pointing(stationName, t, x_eci)
%
% INPUT:
%    stationName[dim]   name of observing station
%    t[dim]             time of observation [s]
%    x_eci[dim]         state in ECI frame
%
% OUTPUT:
%    Az[dim]           azimuth angle (topocentric NW non-inertial) [rad]
%    h[dim]            elevation [rad]
%    rho[dim]          range [km]
%    rho_dot[dim]      range rate [km/s]
%
% CONTRIBUTORS:
%   Davide Zamblera - Politecnico di Milano - Space eng.
%
    % Define station name
    topoFrame = [stationName, '_TOPO'];
    
    % Transformation from ECI to topocentric frame
    ROT_ECI2TOPO = cspice_sxform('J2000', topoFrame, t);
    
    % Compute spacecraft position in ECI
    rv_station_eci = cspice_spkezr(stationName, t, 'J2000', 'NONE', 'EARTH');
    if isequal(size(rv_station_eci),flip(size(x_eci)))
        x_eci = x_eci'; %deal with data in transposed format (each row is a state)
    end
    
    % Compute station-satellite vector in ECI
    rv_station_sat_eci = x_eci - rv_station_eci;
    rv_station_sat_eci = reshape(rv_station_sat_eci,[size(rv_station_sat_eci,1), 1, size(rv_station_sat_eci,2)]);
    
    % Convert state into topocentric frame
    rv_station_sat_topo = pagemtimes(ROT_ECI2TOPO,rv_station_sat_eci);
    rv_station_sat_topo = reshape(rv_station_sat_topo,[size(rv_station_sat_topo,1), size(rv_station_sat_topo,3)]);
    
    % Compute range, azimuth and elevation using cspice_xfmsta
    rll_station_sat = cspice_xfmsta(rv_station_sat_topo,'RECTANGULAR','LATITUDINAL','EARTH');
    
    rho   = rll_station_sat(1,:);   % [km]
    Az = rll_station_sat(2,:);   % [rad]
    h = rll_station_sat(3,:); % [rad]
    rho_dot = rll_station_sat(4,:);   % [km/s]
end

function [ti,te] = compute_visibility(h,t,h_min)
%COMPUTE_VISIBILITY compute epochs where visibility of the s/c change
%
% PROTOTYPE:
%     [ti,te] = compute_visibility(h,t,h_min)
%
% INPUT:
%    h[dim]             elevations for various epochs
%    t[dim]             time of observation [s]
%    h_min[dim]         minimum visible elevation
%
% OUTPUT:
%    ti[dim]            vector of visibility windows start time
%    te[dim]            vector of visibility windows end time
%
% CONTRIBUTORS:
%   Davide Zamblera - Politecnico di Milano - Space eng.
%
    ti = [];
    te = [];
    for i = 1:length(h)-1
        if h(i)<h_min && h(i+1)>h_min
            ti = [ti, t(i+1)/cspice_spd];
        end

        if h(i)>h_min && h(i+1)<h_min
            te = [te, t(i)/cspice_spd];
        end
    end

end

function [reci, veci] = tle2eci_state(tlestr1,tlestr2, et_vec)
%TLE2ECI_STATE from a two line element generate the state at the specified epoch
%
% PROTOTYPE:
%     [reci, veci] = tle2eci_state(tlestr1,tlestr2, et_vec)
%
% INPUT:
%    tlestr1[dim]       first line of TLE
%    tlestr2[dim]       second line of TLE
%    et_vec[dim]        epochs at which state is requested [s]
%
% OUTPUT:
%    reci[dim]          position state in ECI-J2000 frame (each column is an observation)
%    veci[dim]          velocity state in ECI-J2000 frame (each column is an observation)
%
% CONTRIBUTORS:
%   Davide Zamblera - Politecnico di Milano - Space eng.
%
    typerun    = 'u';  % user-provided inputs to SGP4 Matlab function
    opsmode    = 'a';  % afspc approach ('air force space command')
    whichconst =  72;  % WGS72 constants (radius, gravitational parameter)
    arcsec2rad = pi / (180*3600);

    reci = zeros(3,length(et_vec));
    veci = zeros(3,length(et_vec));
    satrec = twoline2rv(tlestr1, tlestr2, typerun,'e', opsmode, whichconst);
    % In satrec we find epoch as (JD = jdsatepoch + jdsatepochf)
    
    % Get TLE epoch
    [year,mon,day,hr,min,sec] = invjday(satrec.jdsatepoch, satrec.jdsatepochf);
    sat_epoch_str = sprintf('%d-%02d-%02dT%02d:%02d:%02.6f', [year,mon,day,hr,min,sec]); %gen string
    sat_epoch_et = cspice_str2et(sat_epoch_str); %from string to TDB/ET

    ateme = [0;0;0];
    %59895  0.185352  0.195074 -0.0160849 -0.0002340 -0.113613 -0.007056  0.000263 -0.000057  37
    ddpsi = -0.113613*arcsec2rad; %  [rad]
    ddeps = -0.007056*arcsec2rad; %  [rad]
    % ATT. OBTAINED BY CELESTRACK EOP DATA - CONSTANT FOR SIMPLICITY
    for i = 1:length(et_vec)
        tepoch = (et_vec(i) - sat_epoch_et)/60.0;
        [~,rteme,vteme] = sgp4(satrec,  tepoch);
        ttt = cspice_unitim(et_vec(i), 'ET', 'TDT')/cspice_jyear()/100; % Epoch we want to perform our conversion TEME to ECI

        [r, v, ~] = teme2eci(rteme, vteme, ateme, ttt, ddpsi, ddeps);
        reci(:,i) = r;
        veci(:,i) = v;
    end

end


function residual = BatchPerthKep(x, t, tref, meas_real, sigma_meas)
%BATCHPERTHKEP compute residuals for lsqnonline, Case: Perth-Keplerian motion
%
% PROTOTYPE:
%     residual = BatchPerthKep(x, t, tref, meas_real, sigma_meas)
%
% INPUT:
%    x[dim]             parameters to estimate, (initial state)
%    t[dim]             vector of measurement epochs
%    tref[dim]          epoch of initial state
%    meas_real[dim]     real measurements (each column is an observation) [s]
%    sigma_meas[dim]    covariance of the measurement instrument
%
% OUTPUT:
%    residual[dim]      vectors of residual = meas_sim - meas_real [3*n_obs x 1]
%
% CONTRIBUTORS:
%   Davide Zamblera - Politecnico di Milano - Space eng.
%
    residual = []; % Initialize output variable
    % Propagate x to the epochs of the measurements
    opts = odeset('RelTol',1e-12,'AbsTol',1e-12);
    GM = cspice_bodvrd('EARTH','GM',1);
    [~,x_prop] = ode113(@TBP,[tref t],x,opts,GM);
    x_prop = x_prop(2:end,:);
    % Compute predicted measurements
    [Az, h, rho] = antenna_pointing('PERTH', t, x_prop);
    %meas_pred = [Az; h; rho];
    % Compute the residual of the measurements and append it to the output
    W_m = diag(1./sigma_meas);
    diff_meas = zeros(size(meas_real));
    diff_meas(1,:) = angdiff(Az, meas_real(1,:));
    diff_meas(2,:) = angdiff(h, meas_real(2,:));
    diff_meas(3,:) = rho - meas_real(3,:);
    diff_meas_weighted = W_m * diff_meas;
    residual = [residual; diff_meas_weighted(:)];
end 

function residual = BatchAllKep(x, windows, tref, measure, sigma_meas)
%BATCHALLKEP compute residuals for lsqnonline, Case: Kourou+Perth-Keplerian motion
%
% PROTOTYPE:
%     residual = BatchAllKep(x, windows, tref, measure, sigma_meas)
%
% INPUT:
%    x[dim]             parameters to estimate, (initial state)
%    t[dim]             vector of measurement epochs
%    tref[dim]          epoch of initial state
%    meas_real[dim]     real measurements (each column is an observation) [s]
%    sigma_meas[dim]    covariance of the measurement instrument
%
% OUTPUT:
%    residual[dim]      vectors of residual = meas_sim - meas_real [3*n_obs x 1]
%
% CONTRIBUTORS:
%   Davide Zamblera - Politecnico di Milano - Space eng.
%
    residual = []; % Initialize output variable
    % Propagate x to the epochs of the measurements
    opts = odeset('RelTol',1e-12,'AbsTol',1e-12);
    GM = cspice_bodvrd('EARTH','GM',1);
    stationNames = {'KOUROU','PERTH'};
    for i = 1:length(windows)
        t = windows{i};
        meas_real = measure{i};
        [~,x_prop] = ode113(@TBP,[tref t],x,opts,GM);
        x_prop = x_prop(2:end,:);
        % Compute predicted measurements
        [Az, h, rho] = antenna_pointing(stationNames{i}, t, x_prop);
        %meas_pred = [Az; h; rho];
        % Compute the residual of the measurements and append it to the output
        W_m = diag(1./sigma_meas);
        diff_meas = zeros(size(meas_real));
        diff_meas(1,:) = angdiff(Az, meas_real(1,:));
        diff_meas(2,:) = angdiff(h, meas_real(2,:));
        diff_meas(3,:) = rho - meas_real(3,:);
        diff_meas_weighted = W_m * diff_meas;
        residual = [residual; diff_meas_weighted(:)];
    end
end 


function residual = BatchAllJ2(x, windows, tref, measure, sigma_meas)
%BATCHALLJ2 compute residuals for lsqnonline, Case: Kourou+Perth-main body gravity + J2
%
% PROTOTYPE:
%     residual = BatchPerthKep(x, t, tref, meas_real, sigma_meas)
%
% INPUT:
%    x[dim]             parameters to estimate, (initial state)
%    t[dim]             vector of measurement epochs
%    tref[dim]          epoch of initial state
%    meas_real[dim]     real measurements (each column is an observation) [s]
%    sigma_meas[dim]    covariance of the measurement instrument
%
% OUTPUT:
%    residual[dim]      vectors of residual = meas_sim - meas_real [3*n_obs x 1]
%
% CONTRIBUTORS:
%   Davide Zamblera - Politecnico di Milano - Space eng.
%
    residual = []; % Initialize output variable
    % Propagate x to the epochs of the measurements
    opts = odeset('RelTol',1e-12,'AbsTol',1e-12);
    GM = cspice_bodvrd('EARTH','GM',1);
    stationNames = {'KOUROU','PERTH'};
    for i = 1:length(windows)
        t = windows{i};
        meas_real = measure{i};
        [~,x_prop] = ode113(@TBP_J2,[tref t],x,opts,GM);
        x_prop = x_prop(2:end,:);
        % Compute predicted measurements
        [Az, h, rho] = antenna_pointing(stationNames{i}, t, x_prop);
        %meas_pred = [Az; h; rho];
        % Compute the residual of the measurements and append it to the output
        W_m = diag(1./sigma_meas);
        diff_meas = zeros(size(meas_real));
        diff_meas(1,:) = angdiff(Az, meas_real(1,:));
        diff_meas(2,:) = angdiff(h, meas_real(2,:));
        diff_meas(3,:) = rho - meas_real(3,:);
        diff_meas_weighted = W_m * diff_meas;
        residual = [residual; diff_meas_weighted(:)];
    end
end


function P_ls = post_process_batch(xest, residual, exitflag, t, meas, sigma_meas, resnorm, jac)
%POST_PROCESS_BATCH generate relevant graphs and compute covariance of
%initial state
%
% PROTOTYPE:
%     P_ls = post_process_batch(xest, residual, exitflag, t, meas, sigma_meas, resnorm, jac)
%
% INPUT:
%    xest[dim]          estimated initial state
%    residual[dim]      vector of residuals from batch filter
%    exitflag[dim]      result flag of batch filter estimation
%    t[dim]             vector of measurement epochs
%    meas[dim]          real measurements (each column is an observation) [s]
%    sigma_meas[dim]    covariance of the measurement instrument
%    resnorm            norm of residual computed by lsqnonlin
%    jac                jacobian computed by lsqnonlin
%
% OUTPUT:
%    P_ls[dim]          covariance of estimated initial state
%
% CONTRIBUTORS:
%   Davide Zamblera - Politecnico di Milano - Space eng.
%
    fprintf('****************************\n')
    fprintf('Results of Batch Filter: Least Squares Minimum Variance\n')
    fprintf('Initial state estimation: x = [')
    fprintf('%g, ',xest)
    fprintf(']\n')
    fprintf('exitflag = %d\n',exitflag)

    % Compute resulting covariance
    Jac = full(jac);
    P_ls = resnorm / (length(residual)-length(xest)) .* inv(Jac'*Jac);
    sigmar = sqrt(trace(P_ls(1:3,1:3)));
    sigmav = sqrt(trace(P_ls(4:6,4:6)));
    fprintf('Trace of upper left matrix sigma_r = %f\n',sigmar);
    fprintf('Trace of lower right matrix sigma_r = %f\n',sigmav);

    figure;
    n_single_meas = length(residual)/3;
    plot(1:n_single_meas,residual(1:3:end),'.')
    hold on;
    plot(1:n_single_meas,residual(2:3:end),'.')
    plot(1:n_single_meas,residual(3:3:end),'.')
    legend('Az','h','rho')

    figure;
    if iscell(t)
        for i = 1:length(t)
            ti = t{i};
            size_meas = size(meas{i});
            n = numel(meas{i}); %get number of el in first
            res_i = residual(1:n);
            residual = residual(n+1:end);
            res_i = reshape(res_i,size_meas);
            res_i = diag(sigma_meas)*res_i;

            % plot of angular difference of prediction and measure in radians 
            hold on;
            plot(ti/cspice_spd(),vecnorm(res_i(1:2,:))*cspice_dpr(),'.');
            

        end
        legend({'Kourou','Perth'})


    else
        residual = reshape(residual, size(meas));
        residual = diag(sigma_meas)*residual;

        % plot of angular difference of prediction and measure in radians 
        plot(t/cspice_spd(),vecnorm(residual(1:2,:))*cspice_dpr(),'.');
        legend('Perth')

    end

    grid on;
    
    title('Angular displacement of predicted w.r.t. measured')
    xlabel('Time of observation [MJD2000]')
    ylabel('$||err_{angle}||$ [deg]')


end