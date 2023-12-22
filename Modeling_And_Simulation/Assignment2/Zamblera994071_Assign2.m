% Modeling and Simulation of Aerospace Systems (2021/2022)
% Assignment # 2
% Author: Davide Zamblera

%% Ex 1
clearvars; close all; clc;

sampl = readmatrix('samples.txt');
t = sampl(2:end,1);
acc = sampl(2:end,2:3);

data.t = t;
xdata = t;
ydata = acc;

data.T = 0.1; %Nm
data.J1 = 0.2; %kgm
data.J2 = 0.1; %kgm
data.lin_damp = 1; %change this to show behaviour with nonlinear damper
k = 0.5;
b = 5;
% k = 3.257516;
% b = 2.927782;
param = [k b];

opts = [];
[t,yout] = ode45(@ode_mech,t,[0 0 0 0]',opts,data,param);
figure;
grid on;
plot(t,yout(:,1),t,yout(:,2));
xlabel('t [s]','Interpreter','latex');
ylabel('Angles [rad]','Interpreter','latex');
legend({'$\theta_1$','$\theta_2$'},'Interpreter','latex')

figure;
grid on;
plot(t,yout(:,3),t,yout(:,4));
xlabel('t [s]','Interpreter','latex');
ylabel('Angular velocities [rad]','Interpreter','latex');
legend({'$\dot{\theta}_1$','$\dot{\theta}_2$'},'Interpreter','latex')

%model identification

opts = optimoptions(@lsqcurvefit,'Display','iter');
[x,resnorm,residual,exitflag,output] = lsqcurvefit(@(param,xdata) paramObj(param,xdata,data),param,xdata,ydata,[],[],opts);
acc_sim = paramObj(x,xdata,data);
figure;
subplot(1,2,1)
plot(t,acc)
grid on;
xlabel('t [s]','Interpreter','latex');
ylabel('$Angular Accelerations [rad/s^2]$','Interpreter','latex');
legend({'$\ddot{\theta}_1$','$\ddot{\theta}_2$'},'Interpreter','latex')
title('Experimental data points')

subplot(1,2,2);
plot(t,acc_sim)
grid on;
xlabel('t [s]','Interpreter','latex');
ylabel('$Angular Accelerations [rad/s^2]$','Interpreter','latex');
legend({'$\ddot{\theta}_1$','$\ddot{\theta}_2$'},'Interpreter','latex')
title('Simulated behaviour')



figure;
plot(t,abs(acc - acc_sim))
title('Differences between simulated and real system','Interpreter','latex')
xlabel('$t [s]$','Interpreter','latex');
ylabel('$|\mathbf{a} - \mathbf{a}_{sim}|$','Interpreter','latex');
grid on;

fprintf('Estimation of elastic spring coefficent: %f [Nm]\n',x(1));
fprintf('Estimation of viscous damping coefficent: %f [kg m s^-1]\n',x(2));

% Non Linear Damper
data.lin_damp = 0;
opts = optimoptions(@lsqcurvefit,'Display','iter');
[x,resnorm,residual,exitflag,output] = lsqcurvefit(@(param,xdata) paramObj(param,xdata,data),param,xdata,ydata,[],[],opts);
acc_sim = paramObj(x,xdata,data);
figure;
subplot(1,2,1)
plot(t,acc)
grid on;
xlabel('t [s]','Interpreter','latex');
ylabel('$Angular Accelerations [rad/s^2]$','Interpreter','latex');
legend({'$\ddot{\theta}_1$','$\ddot{\theta}_2$'},'Interpreter','latex')
title('Experimental data points')

subplot(1,2,2);
plot(t,acc_sim)
grid on;
xlabel('t [s]','Interpreter','latex');
ylabel('$Angular Accelerations [rad/s^2]$','Interpreter','latex');
legend({'$\ddot{\theta}_1$','$\ddot{\theta}_2$'},'Interpreter','latex')
title('Simulated behaviour')



figure;
plot(t,abs(acc - acc_sim))
title('Differences between simulated and real system','Interpreter','latex')
xlabel('$t [s]$','Interpreter','latex');
ylabel('$|\mathbf{a} - \mathbf{a}_{sim}|$','Interpreter','latex');
grid on;

fprintf('Estimation of elastic spring coefficent: %f [Nm]\n',x(1));
fprintf('Estimation of viscous damping coefficent: %f [kg m]\n',x(2));

%% Ex 2
clearvars; close all; clc;

data.R1 = 100;
data.R2k = 10;
data.L = 10;
data.C = 1e-3;
data.vgen = 0;

Vc0 = 1;
y0 = [Vc0 0]';
tspan = [0 1.5];
opts = [];
[t,yout] = ode45(@ode_elec,tspan,y0,opts,data);

figure;
plot(t,yout(:,1),'k')
grid on;
xlabel('t [s]','Interpreter','latex');
ylabel('$V_c \; [V]$','Interpreter','latex');
title('Free dynamics','Interpreter','latex')

data.vgen = 1;
data.f = 5;
tspan = [0 5];
[t,yout] = ode45(@ode_elec,tspan,y0,opts,data);

figure;
plot(t,yout(:,1),'k')
grid on;
xlabel('t [s]','Interpreter','latex');
ylabel('$V_c \; [V]$','Interpreter','latex');
title('Dynamics with Oscillating Source','Interpreter','latex')


%% Ex 3

clearvars; close all; clc;

opts = odeset('RelTol',1e-12,'AbsTol',1e-12);
n_res = 5;
n_el = 2;


A = 1;
l = [0.001 0.024 0 0.012 0.001];
ll_m1 = cumsum([l(1)/2 (l(1) + l(2))/2  (l(2) + l(3))/2 (l(3) + l(4))/2 (l(4) + l(5))/2]);

k = [50 100 1 0.7 150];
rho = [8000 2500];
c = [500 1500];

% compute resistances
R = zeros(1,n_res);
for i = 1:n_res
    R(i) = l(i)/(k(i)*A);
end
R(3) = 0.002; %interface has directly thermal resistance
Rnet = zeros(6,1);
Rnet(1) = R(1)/2;
Rnet(2) = (R(1) + R(2))/2;
Rnet(3) = (R(2) + R(3))/2;
Rnet(4) = (R(3) + R(4))/2;
Rnet(5) = (R(4) + R(5))/2;
Rnet(6) = R(5)/2;

%compute capacitances
C = zeros(1,n_el);
j = [2 4];
for i = 1:n_el
    M = rho(i)*A*l(j(i));
    C(i) = M*c(i);
end

data.R = R;
data.C = C;
data.Rnet = Rnet;

tspan = [0 60];
Tin = 20 ;
y0 = Tin*ones(2,1);
[t,yout] = ode45(@ode_therm,tspan,y0,opts,data);
T = nodes_response(t,yout,data);
figure;
plot(t,T);
grid on;
xlabel('t [s]','Interpreter','latex');
ylabel('T [$^{\circ}$C]','Interpreter','latex');
legend({'$T_1$','$T_2$','$T_{int}$','$T_3$','$T_4$'},'Interpreter','latex')
title('Lumped Model: 5 Points','Interpreter','latex');

[~,idx] = min(abs(t-20));
t_start_m1 = T(idx,:);
[~,idx] = min(abs(t-40));
t_middle_m1 = T(idx,:);
t_fin_m1 = T(end,:);

% CASE WITH 2 NODES for CENTRAL LAYERS

Rnet = zeros(6,1);
Rnet(1) = R(1)/2;
Rnet(2) = R(1)/2 + R(2)/3;
Rnet(3) = R(2)/3;
Rnet(4) = R(2)/3 + R(3)/2;
Rnet(5) = R(3)/2 + R(4)/3;
Rnet(6) = R(4)/3;
Rnet(7) = R(4)/3 + R(5)/2;
Rnet(8) = R(5)/2;

Cnet = zeros(4,1);
Cnet(1) = C(1)/2;
Cnet(2) = C(1)/2;
Cnet(3) = C(2)/2;
Cnet(4) = C(2)/2;

data.R = R;
data.C = C;
data.Rnet = Rnet;
data.Cnet = Cnet;

tspan = [0 60];
Tin = 20 ;
y0 = Tin*ones(4,1);
[t,yout] = ode45(@ode_therm_multi,tspan,y0,opts,data);
T = nodes_response_multi(t,yout,data);
figure;
plot(t,T);
grid on;
xlabel('t [s]','Interpreter','latex');
ylabel('T [$^{\circ}$C]','Interpreter','latex');
legend({'$T_{lin}$','$T_1$','$T_2$','$T_{int}$','$T_3$','$T_4$','$T_{out}$'},'Interpreter','latex')
title('Lumped Model: 7 Points','Interpreter','latex');

% compute temperatures of all the layers for various time instants
ll_m2 = cumsum([l(1)/2, (l(1)/2 + l(2)/3),  l(2)/3, l(2)/3 + l(3)/2, l(3)/2 + l(4)/3, l(4)/3, l(4)/3 + l(5)/2]);
[~,idx] = min(abs(t-20));
t_start_m2 = T(idx,:);
[~,idx] = min(abs(t-40));
t_middle_m2 = T(idx,:);
t_fin_m2 = T(end,:);

% Plot of temperatures evolution
f1 = figure;
firstax = axes (f1); 
plot(ll_m1*1e3,t_start_m1,'ro-',ll_m1*1e3,t_middle_m1,'ro--',ll_m1*1e3,t_fin_m1,'ro-.')
hold on;
plot(ll_m2*1e3,t_start_m2,'ko-',ll_m2*1e3,t_middle_m2,'ko--',ll_m2*1e3,t_fin_m2,'ko-.')
grid on;
legend({'$T_{5-lumps}$','','','$T_{7-lumps}$'},'Interpreter','latex')
secondax = copyobj(firstax, gcf);
delete( get(secondax, 'Children'))

H1 = plot(ll_m2*1e3,t_start_m2, '-', 'Color', [0 0 0]/255, 'Parent', secondax);
H2 = plot(ll_m2*1e3,t_middle_m2, '--', 'Color', [0 0 0]/255, 'Parent', secondax);
H3 = plot(ll_m2*1e3,t_fin_m2, '-.', 'Color', [0 0 0]/255, 'Parent', secondax);
set(secondax, 'Color', 'none', 'XTick', [], 'YAxisLocation', 'right', 'Box', 'Off') 
legend ([H1 H2 H3], {'t = 20s', 't = 40s', 't = 60s'}, 'Location', 'northeast','Interpreter','latex','Color',[1 1 1]);
hold off
set(secondax, 'Color', 'none', 'XTick', [], 'YAxisLocation', 'right', 'Box', 'Off', 'Visible', 'off') %make it transparent
xlabel('x [mm]','Interpreter','latex');
ylabel('T [$^{\circ}$C]','Interpreter','latex');
title('Model Comparison: Temperatures evolution','Interpreter','latex');
%% Ex 4

clearvars; close all; clc;

% Data definition
Km = 20;
R = 200;
data.v0 = 2;
data.om = 5;
data.beta = 0.2;
data.L = 2e-3;
data.J1 = 0.5;
data.J2 = 0.3;
data.b = 0.1;
data.k = 0.5;


v0 = data.v0;
om = data.om;
beta = data.beta;
L = data.L;
J1 = data.J1;
J2 = data.J2;
b = data.b;
k = data.k;

% compute system eigenvalues
A = [-R/L     0      0     -Km/L     0;
      0       0      0      1        0;
      0       0      0      0        1;
      Km/J1  -k/J1   k/J1  -b/J1     b/J1;
      0       k/J2  -k/J2   b/J2    -b/J2];

e = eig(A);
rank(A)

% selection of a stiff integrator and computation of response

tspan = [0 30];
y0 = zeros(5,1);
param = [Km, R];
[t,yout] = ode15s(@ode_elec_mech,tspan,y0,[],data,param);

figure;
plot(t,yout(:,1));
grid on;
xlabel('t [s]','Interpreter','latex');
ylabel('i [A]','Interpreter','latex');
title('Current response','Interpreter','latex');

figure;
plot(t,yout(:,2:3));
grid on;
xlabel('t [s]','Interpreter','latex');
ylabel('Angles [rad]','Interpreter','latex');
legend({'$\theta_1$','$\theta_2$'},'Interpreter','latex')
title('Angular displacements','Interpreter','latex');

figure;
plot(t,yout(:,4:5));
grid on;
xlabel('t [s]','Interpreter','latex');
ylabel('Angular Velocities [rad/s]','Interpreter','latex');
legend({'$\Omega_1$','$\Omega_2$'},'Interpreter','latex')
title('Angular velocities','Interpreter','latex');

% System identification and associated plots

data_exp = readmatrix('Profile.txt');
xdata = data_exp(:,1);
ydata = data_exp(:,2);


param = [Km R]; %first guess of parameters
opts = optimoptions(@lsqcurvefit,'Display','iter');
[x,resnorm,residual,exitflag,output] = lsqcurvefit(@(param,xdata) paramObjOmega(param,xdata,data),param,xdata,ydata,[],[],opts);


[t,yout] = ode15s(@ode_elec_mech,xdata,y0,[],data,x);
om_sim = yout(:,5);
figure;
subplot(1,2,1)
plot(t,ydata)
grid on;
xlabel('t [s]','Interpreter','latex');
ylabel('$\Omega_1 [rad/s]$','Interpreter','latex');
title('Experimental data points')

subplot(1,2,2);
plot(t,om_sim)
grid on;
xlabel('t [s]','Interpreter','latex');
ylabel('$\Omega_1 [rad/s]$','Interpreter','latex');
title('Simulated behaviour')

figure;
plot(t,abs(ydata - om_sim))
title('Differences between simulated and real system','Interpreter','latex')
xlabel('$t [s]$','Interpreter','latex');
ylabel('$|\mathbf{\Omega_1} - \mathbf{\Omega_1}_{sim}|$','Interpreter','latex');
grid on;

fprintf('Estimation of coefficent Km: %f [-]\n',x(1));
fprintf('Estimation of electric resistance: %f [Ohm]\n',x(2));



%% Ex 5
clearvars; close all; clc;
%data definition


data.pump.N = 9; %number of pistons
data.pump.Dp = 0.7e-2; %diameter of pistons
data.pump.dp = 1.5e-2; %distance between opposite pistons/shafts
data.pump.n = 4000/60; %rotation speed
data.pump.Ap = pi*data.pump.Dp^2/4;

data.pilot.lc = 10e-2; %control lever length
data.pilot.th_max = 20; %max control swashplate angle
data.pilot.c = data.pilot.lc*tand(data.pilot.th_max); %max height of control lever
data.pilot.kp = 2.5; %pilot pipe head loss
data.pilot.F0 = 5; %pre loaded force
data.pilot.dk = 1e-3; %diameter of the pipe
data.pilot.Dk = 1e-2; %diameter of the piston
data.pilot.rk = 1; %friction coefficent
data.pilot.mk = 2; %equivalent mass
data.pilot.Ak = pi*data.pilot.Dk^2/4;

data.check_valve.kcv = 2; %coefficent pressure drop check valve
data.check_valve.D = 20e-3; %diameter valve
data.check_valve.A = pi*data.check_valve.D^2/4; %Area valve

data.fluid.rho = 1000; %Density
data.fluid.cw = 4186; %Specific heat water

data.network.tube_data = zeros(6,5);
data.network.tube_data(1,:) = [0.5 1.5 2.7 2.5 1]; %Lengths
data.network.tube_data(2,:) = 0.02*ones(1,5); %Diameters
data.network.tube_data(3,:) = pi*data.network.tube_data(2,:).^2/4; %Areas
data.network.tube_data(4,:) = data.network.tube_data(1,:).*data.network.tube_data(3,:); %Volumes
data.network.tube_data(5,:) = [0.032 0.032 0.040 0.028 0.032]; %Friction factor
data.network.tube_data(6,:) = data.network.tube_data(5,:).*data.network.tube_data(1,:)./data.network.tube_data(2,:);
%V = sum(data.network.tube_data(4,:));
V0 = 0.025;

data.distributor.kd = 15; %pressure drop coefficent distributor
data.distributor.r0 = 5e-3; %radius of the valve
data.distributor.DT = 2; %time to open valve completely

data.filter.kf = 35;
data.filter.D = 20e-3; %diameter filter
data.filter.A = pi*data.filter.D^2/4; %Area filter
data.filter.kl = 0.025; %leaking coefficent

data.heat_exc.K = 1000; %[W]
data.heat_exc.T0 = 350-273.15; %[°C]
data.heat_exc.kT = 20;
data.heat_exc.om = 5;
% compute resistances
n_res = 3;
l = [1 2.5 1]*1e-2;
k = [395 310 125];
Ae = 1000*1e-4;
h = 20;
R = zeros(1,n_res);
for i = 1:n_res
    R(i) = l(i)/(k(i)*Ae);
end
data.heat_exc.Rconv = 1/(h*Ae);
data.heat_exc.Rnet = zeros(1,4);
data.heat_exc.Rnet(1) = R(1)/2;
data.heat_exc.Rnet(2) = (R(1) + R(2))/2;
data.heat_exc.Rnet(3) = (R(2) + R(3))/2;
data.heat_exc.Rnet(4) = R(3)/2;
rho = 8620;
c2 = 100;
data.heat_exc.C = c2*rho*l(2)*Ae;

% Compute capacitance of the fluid

n_tubes = 5;
% total_volume = sum(data.network.tube_data(4,:));
% Mw = data.fluid.rho*total_volume;
D_pipe = 20e-3;
L_pipe = 0.5;
V_pipe = pi*D_pipe^2*L_pipe/4;
data.fluid.Cw = data.fluid.cw*data.fluid.rho*V_pipe;

% design h
pnom = 5*1.01325*1e5; %nominal pressure pump
h = design_pump_control(pnom,data);
fprintf('Spring coefficent of pump control piston: %f [N/m]\n',h)
data.pilot.h = h;

tspan = [0 25];
y0 = [0 0 V0 320-273.15 320-273.15];
[t,yout] = ode45(@ode_hyd,tspan,y0,[],data);

figure(1);
plot(t,yout(:,1)*1e3)
hold on;
%plot([0 t(end)],[data.pilot.c data.pilot.c]*1e3,'--');
title('Position of the pilot piston','Interpreter','latex')
xlabel('t [s]','Interpreter','latex');
ylabel('Position [mm]','Interpreter','latex');
grid on;

figure(1);
curves = gca().Children;
axes('position',[.65 .2 .25 .25])
box on
plot(curves(1).XData,curves(1).YData);
axis tight


figure;
plot(t,yout(:,2)*1e3)
title('Velocity of the pilot piston','Interpreter','latex')
xlabel('t [s]','Interpreter','latex');
ylabel('Velocity [mm/s]','Interpreter','latex');
grid on;


figure;
plot(t,yout(:,3))
title('Volume of the tank','Interpreter','latex')
xlabel('t [s]','Interpreter','latex');
ylabel('Volume [$m^3$]','Interpreter','latex');
grid on;

figure;
plot(t,yout(:,4))
title('Temperature of layer 2 lump','Interpreter','latex')
xlabel('t [s]','Interpreter','latex');
ylabel('Temperature [$^{\circ}$C]','Interpreter','latex');
grid on;

figure;
plot(t,yout(:,5))
title('Temperature of the fluid','Interpreter','latex')
xlabel('t [s]','Interpreter','latex');
ylabel('Temperature [$^{\circ}$C]','Interpreter','latex');
grid on;

parout = zeros(length(t),11);
for i = 1:length(t)
    [~, out] = ode_hyd(t(i),yout(i,:)',data);
    parout(i,:) = out;
end 

figure;
plot(t,parout(:,1:9))
legend({'$p_1$','$p_2$','$p_3$','$p_4$','$p_5$','$p_6$','$p_7$','$p_8$','$p_9$'},'Interpreter','latex')
title('Pressures in the network','Interpreter','latex')
xlabel('t [s]','Interpreter','latex');
ylabel('Pressure [Pa]','Interpreter','latex');
grid on;

figure;
plot(t,parout(:,10))
title('Pump Volumetric Flow rate','Interpreter','latex')
xlabel('t [s]','Interpreter','latex');
ylabel('$Q_{pump}$ [$m^3$/s]','Interpreter','latex');
grid on;


figure;
plot(t,parout(:,11))
title('Heat flow rate to fluid','Interpreter','latex')
xlabel('t [s]','Interpreter','latex');
ylabel('$Q_h$ [W]','Interpreter','latex');
grid on;

%% End of File Functions

function f = ode_mech(t,s,data,param)
%ODE_MECH rhs dynamics of exercise 1 mechanical system
%
% PROTOTYPE:
%   f = ode_mech(t,s,data,param)
%
% INPUT:
%    t[dim]            time [s]
%    s[dim]            state of the system [th1,th2,phi1,phi2] [rad,rad/s]
%    data[dim]         struct containing informations about system characteristics
%                       -data.lin_damp choose between linear or nonlinear damper
%    param[dim]        variable parameters of the system [k,b] [Nm, Nms]
%
% OUTPUT:
%	 f[dim] 	       rhs of the system [unit]
%
% CONTRIBUTORS:
%   Davide Zamblera - Politecnico di Milano - Space eng.
%

    k = param(1);
    b = param(2);

    J1 = data.J1;
    J2 = data.J2;
    T = data.T;

    th1 = s(1,:);
    th2 = s(2,:);
    phi1 = s(3,:);
    phi2 = s(4,:);

    f = zeros(4,length(t));
    f(1,:) = phi1;
    f(2,:) = phi2;
    f(3,:) = - k/J1*th1 + k/J1*th2;
    
    if data.lin_damp
        f(4,:) = - b/J2*phi2 + k/J2*th1 - k/J2*th2 + T/J2;
    else
        f(4,:) = - b/J2*sign(phi2).*phi2.^2 + k/J2*th1 - k/J2*th2 + T/J2;
    end

end

function acc = paramObj(param,xdata,data)
%PARAMOBJ function needed by lsqcurvefit
%
% PROTOTYPE:
%   F = ode_elec(t,s,data)
%
% INPUT:
%    param[dim]        variable parameters of the system [k,b] [Nm, Nms]
%    xdata[dim]        time vector of experiment result [s]
%    data[dim]         output matrix of experiment result [rad/s^2]
%
% OUTPUT:
%	 F[dim] 	       rhs of the system [unit]
%
% CONTRIBUTORS:
%   Davide Zamblera - Politecnico di Milano - Space eng.
%
    t = xdata;
    opts = [];
    [t,yout] = ode45(@ode_mech,t,[0 0 0 0]',opts,data,param);

    acc = ode_mech(t',yout',data,param)';
    acc = acc(:,3:4);

end

function F = ode_elec(t,s,data)
%ODE_ELEC rhs dynamics of exercise 2 electrical system
%
% PROTOTYPE:
%   F = ode_elec(t,s,data)
%
% INPUT:
%    t[dim]            time [s]
%    s[dim]            state of the system [Vc,dVc] [V,V/s]
%    data[dim]         struct containing informations about system characteristics
%
% OUTPUT:
%	 F[dim] 	       rhs of the system [unit]
%
% CONTRIBUTORS:
%   Davide Zamblera - Politecnico di Milano - Space eng.
%
    R1 = data.R1;
    R2k = data.R2k;
    L = data.L;
    C = data.C;
    vgen = data.vgen;
    

    Vc = s(1);
    dVc = s(2);
    R2 = R2k*C*dVc;
    K = L*C*(2*R2/R1 + 1); %2 due to R=R2k*i

    F = zeros(2,1);

    F(1) = dVc;

    if vgen
        f = data.f;
        V = sin(2*pi*f*t)*atan(t);
        dV = sin(2*pi*f*t)/(t^2 + 1) + 2*f*pi*atan(t)*cos(2*pi*f*t);
        F(2) = -(L/R1 + R2*C)/K * dVc - 1/K*Vc + 1/K*(V + L/R1*dV); 
    else
        F(2) = -(L/R1 + R2*C)/K * dVc - 1/K*Vc;
    end


end


function f = ode_therm(t,s,data)
%ODE_THERM rhs dynamics of exercise 3 thermal system
%
% PROTOTYPE:
%   f = ode_elec(t,s,data)
%
% INPUT:
%    t[dim]            time [s]
%    s[dim]            state of the system [T1,T2] [°C]
%    data[dim]         struct containing informations about system characteristics
%
% OUTPUT:
%	 f[dim] 	       rhs of the system [unit]
%
% CONTRIBUTORS:
%   Davide Zamblera - Politecnico di Milano - Space eng.
%
    T1 = s(1);
    T2 = s(2);

    if t <= 1
        T_ramp1 = 20 ;
        T_ramp2 = 1000 ;
        dt = 1;
        Ti = T_ramp1 + (T_ramp2 - T_ramp1)/dt*t;
    else
        Ti = 1000 ;
    end
    To = 20 ;
    C = data.C;
    Rnet = data.Rnet;

    f = zeros(2,1);
    f(1) = (Ti - T1)/((Rnet(1) + Rnet(2))*C(1)) - (T1 - T2)/((Rnet(3) + Rnet(4))*C(1));
    f(2) = (T1 - T2)/((Rnet(3) + Rnet(4))*C(2)) - (T2 - To)/((Rnet(5) + Rnet(6))*C(2));

end

function T = nodes_response(t,yout,data)
%NODES_RESPONSE_MULTI identify temperatures which are not in the state
%
% PROTOTYPE:
%   T = nodes_response_multi(t,yout,data)
%
% INPUT:
%    t[dim]            time [s]
%    yout[dim]         output of thermal model system dynamics [°C]
%    data[dim]         struct containing informations about system characteristics
%
% OUTPUT:
%    T[dim]            vector of temperatures [°C]
%
% CONTRIBUTORS:
%   Davide Zamblera - Politecnico di Milano - Space eng.
%
    T1 = yout(:,1);
    T2 = yout(:,2);

    T_ramp1 = 20 ;
    T_ramp2 = 1000 ;
    dt = 1;
    
    n = length(t);
    Ti = zeros(n,1);
    for i = 1:n
        if t(i) <= 1
            Ti(i) = T_ramp1 + (T_ramp2 - T_ramp1)/dt*t(i);
        else
            Ti(i) = T_ramp2;
        end
    end
    To = (20 )*ones(n,1);

    Rnet = data.Rnet;
    Q1 = (Ti - T1)/(Rnet(1) + Rnet(2));
    Tlin = Ti - Rnet(1)*Q1; %compute inner lining temperatures
    Q2 = (T1 - T2)/(Rnet(3) + Rnet(4));
    Tint = T1 - Rnet(3)*Q2; %compute interface temperatures
    Q3 = (T2 - To)/(Rnet(5) + Rnet(6));
    Twall = T2 - Rnet(5)*Q3; %compute outer wall temperatures

    T = [Tlin, T1, Tint, T2, Twall];
end


function f = ode_therm_multi(t,s,data)
%ODE_THERM_MULTI rhs dynamics of exercise 3 thermal system, with multi-nodal layers
%
% PROTOTYPE:
%   f = ode_therm_multi(t,s,data)
%
% INPUT:
%    t[dim]            time [s]
%    s[dim]            state of the system [T1,T2,T3,T4] [°C]
%    data[dim]         struct containing informations about system characteristics
%
% OUTPUT:
%	 f[dim] 	       rhs of the system [unit]
%
% CONTRIBUTORS:
%   Davide Zamblera - Politecnico di Milano - Space eng.
%
    T1 = s(1);
    T2 = s(2);
    T3 = s(3);
    T4 = s(4);

    if t <= 1
        T_ramp1 = 20 ;
        T_ramp2 = 1000 ;
        dt = 1;
        Ti = T_ramp1 + (T_ramp2 - T_ramp1)/dt*t;
    else
        Ti = 1000 ;
    end
    To = 20 ;
    Cnet = data.Cnet;
    Rnet = data.Rnet;

    f = zeros(2,1);
    f(1) = (Ti - T1)/((Rnet(1) + Rnet(2))*Cnet(1)) - (T1 - T2)/(Rnet(3)*Cnet(1));
    f(2) = (T1 - T2)/(Rnet(3)*Cnet(2)) - (T2 - T3)/((Rnet(4) + Rnet(5))*Cnet(2));
    f(3) = (T2 - T3)/((Rnet(4) + Rnet(5))*Cnet(3)) - (T3 - T4)/(Rnet(6)*Cnet(3));
    f(4) = (T3 - T4)/(Rnet(6)*Cnet(4)) - (T4 - To)/((Rnet(7) + Rnet(8))*Cnet(4));

end

function T = nodes_response_multi(t,yout,data)
%NODES_RESPONSE_MULTI identify temperatures which are not in the state for
%the multiple node model
%
% PROTOTYPE:
%   T = nodes_response_multi(t,yout,data)
%
% INPUT:
%    t[dim]            time [s]
%    yout[dim]         output of thermal model system dynamics [°C]
%    data[dim]         struct containing informations about system characteristics
%
% OUTPUT:
%    T[dim]            vector of temperatures [°C]
%
% CONTRIBUTORS:
%   Davide Zamblera - Politecnico di Milano - Space eng.
%
    T1 = yout(:,1);
    T2 = yout(:,2);
    T3 = yout(:,3);
    T4 = yout(:,4);

    T_ramp1 = 20 ;
    T_ramp2 = 1000 ;
    dt = 1;
    
    n = length(t);
    Ti = zeros(n,1);
    for i = 1:n
        if t(i) <= 1
            Ti(i) = T_ramp1 + (T_ramp2 - T_ramp1)/dt*t(i);
        else
            Ti(i) = T_ramp2;
        end
    end
    To = (20 )*ones(n,1);

    Rnet = data.Rnet;
    Q1 = (Ti - T1)/(Rnet(1) + Rnet(2));
    Tlin = Ti - Rnet(1)*Q1; %compute inner lining temperatures
    Q2 = (T2 - T3)/(Rnet(4) + Rnet(5));
    Tint = T2 - Rnet(4)*Q2; %compute interface temperatures
    Q3 = (T4 - To)/(Rnet(7) + Rnet(8));
    Twall = T4 - Rnet(7)*Q3; %compute outer wall temperatures

    T = [Tlin, T1, T2, Tint, T3, T4, Twall];
end

function f = ode_elec_mech(t,s,data,param)
%ODE_ELEC_MECH rhs dynamics of exercise 4 multi domain electrical mechanical system
%
% PROTOTYPE:
%   f = ode_elec_mech(t,s,data,param)
%
% INPUT:
%    t[dim]            time [s]
%    s[dim]            state of the system [i,th1,th2,Om1,Om2] [A,rad,rad/s]
%    data[dim]         struct containing informations about system characteristics
%
% OUTPUT:
%	 f[dim] 	       rhs of the system [unit]
%
% CONTRIBUTORS:
%   Davide Zamblera - Politecnico di Milano - Space eng.
%
    
    Km = param(1);
    R = param(2);

    v0 = data.v0;
    om = data.om;
    beta = data.beta;
    L = data.L;
    J1 = data.J1;
    J2 = data.J2;
    b = data.b;
    k = data.k;


    i = s(1);
    th1 = s(2);
    th2 = s(3);
    Om1 = s(4);
    Om2 = s(5);

    v = v0*cos(om*t)*exp(-beta*t);

    f = zeros(5,1);
    f(1) = -R/L*i -Km/L*Om1 + v/L;
    f(2) = Om1;
    f(3) = Om2;
    f(4) = Km/J1*i - k/J1*th1 + k/J1*th2 - b/J1*Om1 + b/J1*Om2;
    f(5) = k/J2*th1 - k/J2*th2 + b/J2*Om1 - b/J2*Om2;

end


function om = paramObjOmega(param,xdata,data)
%PARAMOBJ function needed by lsqcurvefit
%
% PROTOTYPE:
%   F = ode_elec(t,s,data)
%
% INPUT:
%    param[dim]        variable parameters of the system [Km,R] [-, Ohm]
%    xdata[dim]        time vector of experiment result [s]
%    data[dim]         output matrix of experiment result [rad/s]
%
% OUTPUT:
%	 F[dim] 	       rhs of the system [unit]
%
% CONTRIBUTORS:
%   Davide Zamblera - Politecnico di Milano - Space eng.
%
    t = xdata;
    opts = [];
    y0 = zeros(5,1);
    [~,yout] = ode15s(@ode_elec_mech,t,y0,opts,data,param);

    om = yout(:,5);

end

% Exercise 5 functions

function Q1 = pump(t,x,data)
%PUMP subfunction of ODE_HYD
%
% INPUT:
%    x[dim]            position of the piston [m]
%    data[dim]         struct containing informations about system characteristics
%
% OUTPUT:
%    Q1[dim]           volumetric flow rate of the pump [m^3/s]
%
% CONTRIBUTORS:
%   Davide Zamblera - Politecnico di Milano - Space eng.
%
    N = data.pump.N; %number of pistons
    dp = data.pump.dp; %distance between opposite pistons
    lc = data.pilot.lc; %control lever length
    n = data.pump.n; %rotation speed
    c = data.pilot.c; %max height of control lever
    Ap = data.pump.Ap; %Area of a single piston
    s = dp/lc*(c - x); %stroke

    % pump start up
    %tr = 2;
    %if t < tr
    %  n = n*t/tr;
    %end
    
    Q1 = n*N*Ap*s;

end

function p1 = check_valve(p2,Q,data)
%CHECK_VALVE subfunction of ODE_HYD
%
% INPUT:
%    p2[dim]           pressure at the exit point [Pa]
%    Q[dim]            volume flow rate through the element [m^3/s]
%    data[dim]         struct containing informations about system characteristics
%
% OUTPUT:
%    p1[dim]           pressure at the entrance point [Pa]
%
% CONTRIBUTORS:
%   Davide Zamblera - Politecnico di Milano - Space eng.
%
    rho = data.fluid.rho;
    kcv = data.check_valve.kcv;
    A = data.check_valve.A;
    p1 = p2 + 1/2*kcv*rho*Q*abs(Q)/A^2;
end

function p1 = pdrop(p2,Q,label,data)
%PDROP subfunction of ODE_HYD
%
% INPUT:
%    p2[dim]           pressure at the exit point [Pa]
%    Q[dim]            volume flow rate through the element [m^3/s]
%    label[dim]        tube label
%    data[dim]         struct containing informations about system characteristics
%
% OUTPUT:
%    p1[dim]           pressure at the entrance point [Pa]
%
% CONTRIBUTORS:
%   Davide Zamblera - Politecnico di Milano - Space eng.
%
    rho = data.fluid.rho;
    p1 = p2 - 1/2*tube('k',label,data)*rho*abs(Q)*Q/tube('A',label,data)^2;
end

function val = tube(datatype,label,data)
%TUBE subfunction of ODE_HYD: check tube database and output requested data
%
% INPUT:
%    datatype[dim]     variable that user is searching 
%    ['L','D','A','V','f','k'] = [Length,Diameter,Area,Volume,friction factor, pressure drop constant]
%    label[dim]        tube label, att. direction counts [k12 = -k21]
%    data[dim]         struct containing informations about system characteristics
%
% OUTPUT:
%    val[dim]          output requested value [-]
%
% CONTRIBUTORS:
%   Davide Zamblera - Politecnico di Milano - Space eng.
%
    datatype_list = {'L','D','A','V','f','k'};
    label_list = {'T1','34','56','78','9T'};
    inv_label_list = {'1T','43','65','87','T9'};
    m = strcmp(datatype_list,datatype);
    flag = 1;
    n = strcmp(label_list,label);
    if all(~n)
        n = strcmp(inv_label_list,label);
        flag = -1;
    end
    tube_data = data.network.tube_data;
    tube_data(6,:) = flag*tube_data(6,:);
    val = tube_data(m,n);

end


function [xdot,vdot] = pump_control(x,v,p,data)
%PUMP_CONTROL subfunction of ODE_HYD
%
% INPUT:
%    x[dim]            position of the pilot piston [m]
%    v[dim]            position of the pilot piston [m/s]
%    p[dim]            pressure at pump exit point [Pa]
%    data[dim]         struct containing informations about system characteristics
%
% OUTPUT:
%    xdot[dim]         derivative of the position of the pilot piston [m/s]
%    vdot[dim]         derivative of the velocity of the pilot piston [m/s^2]
%
% CONTRIBUTORS:
%   Davide Zamblera - Politecnico di Milano - Space eng.
%
    rho = data.fluid.rho;

    kp = data.pilot.kp; %pilot pipe head loss
    F0 = data.pilot.F0; %pre loaded force
    dk = data.pilot.dk; %diameter of the pipe
    Dk = data.pilot.Dk; %diameter of the piston
    rk = data.pilot.rk; %friction coefficent
    h = data.pilot.h; %spring constant
    mk = data.pilot.mk; %equivalent mass
    Ak = data.pilot.Ak;

    xdot = v;
    vk = v*(Dk/dk)^2;
    pk = p - 1/2*kp*rho*vk*abs(vk);
    vdot = 1/mk*(pk*Ak - F0 - h*x - rk*v);
end

function h = design_pump_control(pnom,data)
%DESIGN_PUMP_CONTROL compute the spring elastic coefficent associated with
%an equilibrium at pnom
%
% INPUT:
%    pnom[dim]         nominal pressure to achieve at equilibrium [Pa]
%    data[dim]         struct containing informations about system characteristics
%
% OUTPUT:
%    h[dim]            spring elastic coefficent [N/m]
%
% CONTRIBUTORS:
%   Davide Zamblera - Politecnico di Milano - Space eng.
%
    F0 = data.pilot.F0; %pre loaded force
    Ak = data.pilot.Ak;
    c = data.pilot.c; %max height of control lever
    h = (pnom*Ak - F0)/c;
end

function p1 = distributor(t,p2,Q,data)
%FILTER subfunction of ODE_HYD
%
% INPUT:
%    t[dim]            time [s]
%    p2[dim]           pressure at the exit point [Pa]
%    Q[dim]            volume flow rate through the element [m^3/s]
%    data[dim]         struct containing informations about system characteristics
%
% OUTPUT:
%    p1[dim]           pressure at the entrance point [Pa]
%
% CONTRIBUTORS:
%   Davide Zamblera - Politecnico di Milano - Space eng.
%
    rho = data.fluid.rho;
    kd = data.distributor.kd; %pressure drop coefficent distributor
    r0 = data.distributor.r0; %radius of the valve
    DT = data.distributor.DT;
    if t <= DT
        u = t/DT;
        %alfa = pi + pi*t/DT;
    else
        u = 1;
        %alfa = 2*pi;
    end
    alfa = pi + 2*acos(1 - abs(u));
    A = r0^2/2*(alfa - sin(alfa)); %open area of the valve
    p1 = p2 + 1/2*kd*rho*Q*abs(Q)/A^2;
end

function p1 = filter(p2,Q,data)
%FILTER subfunction of ODE_HYD
%
% INPUT:
%    p2[dim]           pressure at the exit point [Pa]
%    Q[dim]            volume flow rate through the element [m^3/s]
%    data[dim]         struct containing informations about system characteristics
%
% OUTPUT:
%    p1[dim]           pressure at the entrance point [Pa]
%
% CONTRIBUTORS:
%   Davide Zamblera - Politecnico di Milano - Space eng.
%
    rho = data.fluid.rho;
    kf = data.filter.kf; %pressure drop coefficent filter 
    A = data.filter.A;
    p1 = p2 + 1/2*kf*rho*Q*abs(Q)/A^2;
end

function [Q2,Qleak] = filter_leak(Q1,data)
%FILTER_LEAK subfunction of ODE_HYD
%
% INPUT:
%    Q1[dim]          volumetric flow entering the filter [m^3/s]
%    data[dim]        struct containing informations about system characteristics
%
% OUTPUT:
%    Q2[dim]          volumetric flow exiting the filter [m^3/s]
%    Qleak[dim]       volumetric flow of the leak [m^3/s]
%
% CONTRIBUTORS:
%   Davide Zamblera - Politecnico di Milano - Space eng.
%
    kl = data.filter.kl; %leaking coefficent
    Qleak = kl*Q1;
    Q2 = Q1 - Qleak;
end

function Vdot = tank(Qin, Qout)
%TANK subfunction of ODE_HYD
%
% INPUT:
%    Qin[dim]          volumetric flow entering the tank [m^3/s]
%    Qout[dim]         volumetric flow exiting the tank [m^3/s]
%
% OUTPUT:
%    Vdot[dim]         derivative of the volume of fluid in tank [m^3/s]
%
% CONTRIBUTORS:
%   Davide Zamblera - Politecnico di Milano - Space eng.
%
    Vdot = -Qout + Qin;
end

function p1 = rayleigh(p2,Qdot,data)
%RAYLEIGH subfunction of ODE_HYD
%
% INPUT:
%    p2[dim]           pressure at the exit point [Pa]
%    Qdot[dim]         heat flux from the heat exchanger [W]
%    data[dim]         struct containing informations about system characteristics
%
% OUTPUT:
%    p1[dim]           pressure at the entrance point [Pa]
%
% CONTRIBUTORS:
%   Davide Zamblera - Politecnico di Milano - Space eng.
%
    K = data.heat_exc.K;
    p1 = p2/(exp(Qdot/K));
end

function [Tdot, Tfluid_dot, Qh] = heat_exc(t,T2,Tfluid,data)
%HEAT_EXC subfunction of ODE_HYD
%
% INPUT:
%    t[dim]            time [s]
%    T2[dim]           temperature of the lumped second layer [°C]
%    Tfluid[dim]       temperature of the fluid [°C]
%    data[dim]         struct containing informations about system characteristics
%
% OUTPUT:
%    Tdot[dim]         temperature T2 derivative [°C/s]
%    Tfluid_dot[dim]   fluid temperature derivative [°C/s]
%    Qh[dim]           heat flux exiting from heat exchanger [W]
%
% CONTRIBUTORS:
%   Davide Zamblera - Politecnico di Milano - Space eng.
%
    T0 = data.heat_exc.T0;
    kT = data.heat_exc.kT;
    om = data.heat_exc.om;
    T = T0 + kT*cos(om*t);
    
    Rconv = data.heat_exc.Rconv;
    Rnet = data.heat_exc.Rnet;
    C = data.heat_exc.C;
    Cw = data.fluid.Cw;
    Qin = 1/(Rnet(1) + Rnet(2))*(T - T2);
    Qh = 1/(Rnet(3) + Rnet(4) + Rconv)*(T2 - Tfluid);
    Tdot = (Qin - Qh)/C;
    Tfluid_dot = Qh/Cw;

end


function [f, parout] = ode_hyd(t,s,data)
%ODE_HYD rhs dynamics of exercise 5 thermal-hydraulic system
%
% PROTOTYPE:
%   [f, parout] = ode_hyd(t,s,data)
%
% INPUT:
%    t[dim]            time [s]
%    s[dim]            state of the system [m,m/s,m^3,°C]
%    data[dim]         struct containing informations about system characteristics
%
% OUTPUT:
%	 f[dim] 	       rhs of the system [unit]
%
% CONTRIBUTORS:
%   Davide Zamblera - Politecnico di Milano - Space eng.
%
    %state
    x = s(1);
    v = s(2);
    %Vt = s(3);
    T2 = s(4);
    Tfluid = s(5);
    c = data.pilot.c; %max height of control lever
    if x < 0
        x = 0;
    end
    if x <= 0 && v < 0
        x = 0;
        v = 0;
    end
    if x > c
        x = c;
    end
    if x >= c && v > 0
        x = c;
        v = 0;
    end


    pt = 0.1e6;
    Q1 = pump(t,x,data);
    %conservation of mass
    Q2 = Q1;
    Q3 = Q2;
    Q4 = Q3;
    Q5 = Q4;
    Q6 = Q5;
    
    %compute leak
    [Q7,Qleak] = filter_leak(Q6,data);
    Q8 = Q7;
    Q9 = Q8;

    Qin = Q9 + Qleak;
    Qout = Q1;
    Vdot = tank(Qin,Qout);

    %compute pressures
    p1 = pdrop(pt,Q1,'T1',data);
    p9 = pdrop(pt,Q9,'T9',data);

    %compute heat flow Tdot
    [T2dot, Tfluid_dot, Qh] = heat_exc(t,T2,Tfluid,data);

    %continue temperature computation
    p8 = rayleigh(p9,Qh,data);
    p7 = pdrop(p8,Q8,'87',data);
    p6 = filter(p7,Q6,data);
    p5 = pdrop(p6,Q6,'65',data);
    p4 = distributor(t,p5,Q5,data);
    p3 = pdrop(p4,Q4,'43',data);
    p2 = check_valve(p3,Q3,data);

    [xdot,vdot] = pump_control(x,v,p2,data);
    % pilot piston physical constraints
    if (x == c && vdot > 0) || (x == 0 && vdot < 0)
        xdot = 0;
        vdot = 0;
    end

    

    f = zeros(4,1);
    f(1) = xdot;
    f(2) = vdot;
    f(3) = Vdot;
    f(4) = T2dot;
    f(5) = Tfluid_dot;
    parout = [p1 p2 p3 p4 p5 p6 p7 p8 p9 Q1 Qh];
end