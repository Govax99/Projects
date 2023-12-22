% Spacecraft Guidance and Navigation (2020/2021)
% Assignment # 1
% Author: Davide Zamblera

%% Section 1
clearvars; close all; clc;

%constants
mu = 1.21506683e-2;
const.mu = mu;

ms = 3.28900541e5; %Scaled mass of the Sun
const.ms = ms;

rho = 3.88811143e2; %Scaled Sun–(Earth + Moon) distance
const.rho = rho;

oms = -9.25195985e-1; %Scaled angular velocity of the Sun
const.oms = oms;

k = 1e3; % from [km]->[m]
Re = k*6378; %m Mean Earth radius
const.Re = Re;

Rm = k*1738; %m Mean Moon radius
const.Rm = Rm;

hi = k*167; %m Altitude of departure orbit
hf = k*100; %m Altitude of arrival orbit

DU = 3.84405000e8; %m Distance unit, from [m]->[adimensional]

ri = (Re + hi)/DU;
rf = (Rm + hf)/DU;
const.ri = ri;
const.rf = rf;

%data initialization
alfa = 1.5*pi;
beta = 1.41;
delta = 7;
ti = 0;
prm = [alfa, beta, delta, ti]';

s = param2state(prm,const);
x0 = s(1:4);
tspan = [s(5) s(6)];

yMoon = [1-mu 0 0 0];
yEarth = [-mu 0 0 0];

%integration of the state, obtain trajectory
[tt,y] = PBRFBP_traj(x0,tspan,const);

%plot results: rotating frame - origin baricenter earth moon
figure;
plot(y(:,1),y(:,2),'k')

hold on;
plot(yEarth(1),yEarth(2),'ok','MarkerFaceColor','red')
plot(yMoon(1),yMoon(2),'ok','MarkerFaceColor','blue')
axis equal;
grid on;
xlabel('x [-]','Interpreter','latex');
ylabel('y [-]','Interpreter','latex');
zlabel('z [-]','Interpreter','latex');
title('Rotating Earth-Moon frame @Earth-Moon baricenter, adimensional',Interpreter="latex")
legend({'Guess Trajectory','Earth','Moon'},'Interpreter','latex');
axis(axis+sign(axis)*0.2)




%transformation of state in new frame
YEarth = [0 0 0 0];
YMoon = rot2in_ctr_earth(tt(end),yMoon,mu);
Y = rot2in_ctr_earth(tt,y,mu);

%plot results: inertial frame - origin earth center
figure;
plot(Y(:,1),Y(:,2),'k')
hold on;
plot(YEarth(1),YEarth(2),'ok','MarkerFaceColor','red')
plot(YMoon(1),YMoon(2),'ok','MarkerFaceColor','blue')
axis equal;
grid on;
xlabel('x [-]','Interpreter','latex');
ylabel('y [-]','Interpreter','latex');
zlabel('z [-]','Interpreter','latex');
title('Inertial Earth fixed frame @Earth Center, adimensional',Interpreter="latex")
legend({'Guess Trajectory','Earth','Moon at $t_f=7$'},'Interpreter','latex');
axis(axis+sign(axis)*0.2)

%% Section 2
clearvars; close all; clc;

%constants
mu = 1.21506683e-2;
const.mu = mu;

ms = 3.28900541e5; %Scaled mass of the Sun
const.ms = ms;

rho = 3.88811143e2; %Scaled Sun–(Earth + Moon) distance
const.rho = rho;

oms = -9.25195985e-1; %Scaled angular velocity of the Sun
const.oms = oms;

k = 1e3; % from [km]->[m]
Re = k*6378; %m Mean Earth radius
const.Re = Re;

Rm = k*1738; %m Mean Moon radius
const.Rm = Rm;

hi = k*167; %m Altitude of departure orbit
hf = k*100; %m Altitude of arrival orbit

DU = 3.84405000e8; %m Distance unit, from [m]->[adimensional]

ri = (Re + hi)/DU;
rf = (Rm + hf)/DU;
const.ri = ri;
const.rf = rf;

%data initialization
alfa = 1.5*pi;
beta = 1.41;
delta = 7;
ti = 0;
prm = [alfa, beta, delta, ti]';

s = param2state(prm,const);

%solve simple shooting without gradients/jacobians
opts = optimoptions(@fmincon,'Display','iter','Algorithm','active-set','MaxFunctionEvaluations',2000);
tic;
[yopt,fval_smp] = simple_shoot(s,const,opts);
toc

%plot initial guess
x0 = s(1:4);
tspan = [s(5) s(6)];
[~,y] = PBRFBP_traj(x0,tspan,const);
yMoon = [1-mu 0 0 0];
yEarth = [-mu 0 0 0];
figure;
plot(y(:,1),y(:,2),'r') %plot initial guess - comparison
hold on;

%plot results: rotating frame - origin baricenter earth moon
[~,y] = PBRFBP_traj(yopt(1:4),yopt(5:6),const);
plot(y(:,1),y(:,2),'k')


%plot astronomical objects
plot(yEarth(1),yEarth(2),'ok','MarkerFaceColor','red')
plot(yMoon(1),yMoon(2),'ok','MarkerFaceColor','blue')
legend({'Guess Trajectory','Solution with Gradients','Earth','Moon'},'Interpreter','latex');

%finish the plot
axis equal;
grid on;
legend({'Guess Trajectory','Solution','Earth','Moon'},'Interpreter','latex');
title('Rotating Earth-Moon frame @Earth-Moon baricenter, adimensional',Interpreter="latex")
xlabel('x [-]','Interpreter','latex');
ylabel('y [-]','Interpreter','latex');
zlabel('z [-]','Interpreter','latex');
axis(axis+sign(axis)*0.2)




%solve simple shooting with gradients/jacobians
opts = optimoptions(@fmincon,'Display','iter-detailed','Algorithm','active-set',...
    'MaxFunctionEvaluations',2000,'SpecifyConstraintGradient',true,'SpecifyObjectiveGradient',true);
tic;
[yopt_grd,fval_grd_smp] = simple_shoot(s,const,opts);
toc

%plot initial guess
x0 = s(1:4);
tspan = [s(5) s(6)];
[~,y] = PBRFBP_traj(x0,tspan,const);
yMoon = [1-mu 0 0 0];
yEarth = [-mu 0 0 0];
figure;
plot(y(:,1),y(:,2),'r') %plot initial guess - comparison
hold on;

%plot results: rotating frame - origin baricenter earth moon
[~,y] = PBRFBP_traj(yopt_grd(1:4),yopt_grd(5:6),const);
plot(y(:,1),y(:,2),'k')


%plot astronomical objects
plot(yEarth(1),yEarth(2),'ok','MarkerFaceColor','red')
plot(yMoon(1),yMoon(2),'ok','MarkerFaceColor','blue')
legend({'Guess Trajectory','Solution with Gradients','Earth','Moon'},'Interpreter','latex');

%finish the plot
axis equal;
grid on;
legend({'Guess Trajectory','Solution with gradients','Earth','Moon'},'Interpreter','latex');
title('Rotating Earth-Moon frame @Earth-Moon baricenter, adimensional',Interpreter="latex")
xlabel('x [-]','Interpreter','latex');
ylabel('y [-]','Interpreter','latex');
zlabel('z [-]','Interpreter','latex');
axis(axis+sign(axis)*0.2)
% you can implement multiple optimization rounds 

%% Section 3
clearvars; close all; clc;

%constants
mu = 1.21506683e-2;
const.mu = mu;

ms = 3.28900541e5; %Scaled mass of the Sun
const.ms = ms;

rho = 3.88811143e2; %Scaled Sun–(Earth + Moon) distance
const.rho = rho;

oms = -9.25195985e-1; %Scaled angular velocity of the Sun
const.oms = oms;

k = 1e3; % from [km]->[m]
Re = k*6378; %m Mean Earth radius
const.Re = Re;

Rm = k*1738; %m Mean Moon radius
const.Rm = Rm;

hi = k*167; %m Altitude of departure orbit
hf = k*100; %m Altitude of arrival orbit

DU = 3.84405000e8; %m Distance unit, from [m]->[adimensional]
const.DU = DU;

ri = (Re + hi)/DU;
rf = (Rm + hf)/DU;
const.ri = ri;
const.rf = rf;

%data initialization
alfa = 1.5*pi;
beta = 1.41;
delta = 7;
ti = 0;
prm = [alfa, beta, delta, ti]';

s = param2state(prm,const);
x0 = s(1:4);

%%%%%%%%%%%%%% N = 4 %%%%%%%%%%%%%%%
%prepare the initial state
N = 4;
t1 = ti;
tN = ti+delta;

tspan = [t1 tN];
[~,x] = sample_trajectory(tspan,x0,N,const);


%solve multiple shooting with gradients/jacobians
opts = optimoptions(@fmincon,'Display','iter-detailed','Algorithm','active-set',...
    'SpecifyConstraintGradient',true,'SpecifyObjectiveGradient',true);
tic;
[yopt_mlt,~] = multi_shoot(x,const,opts);
toc

%plot the results
figure;
t1 = yopt_mlt(end-1);
tN = yopt_mlt(end);
t = zeros(N,1);
hold on;
axis equal;
grid on;
xlabel('x [-]','Interpreter','latex');
ylabel('y [-]','Interpreter','latex');
zlabel('z [-]','Interpreter','latex');

for i = 1:N-1
    t(i) = t1 + (i-1)/(N-1)*(tN-t1);
    t(i+1) = t1 + i/(N-1)*(tN-t1);
    [~,y] = PBRFBP_traj(yopt_mlt(4*(i-1)+1:4*i),[t(i) t(i+1)],const);
    plot(y(:,1),y(:,2),'k')
    hold on;
    plot(yopt_mlt(4*(i-1)+1),yopt_mlt(4*(i-1)+2),'o','MarkerFaceColor','red')
end
plot(yopt_mlt(end-5),yopt_mlt(end-4),'o','MarkerFaceColor','red')
axis(axis+sign(axis)*0.2)
title('Rotating Earth-Moon frame @Earth-Moon baricenter, adimensional',Interpreter="latex")

%%%%%%%%%%%%%% N = 10 %%%%%%%%%%%%%%%
%prepare the initial state
N = 10;
t1 = ti;
tN = ti+delta;

tspan = [t1 tN];
[~,x] = sample_trajectory(tspan,x0,N,const);


%solve multiple shooting with gradients/jacobians
tic;
[yopt_mlt,~] = multi_shoot(x,const,opts);
toc

%plot the results
figure;
t1 = yopt_mlt(end-1);
tN = yopt_mlt(end);
t = zeros(N,1);
hold on;
axis equal;
grid on;
xlabel('x [-]','Interpreter','latex');
ylabel('y [-]','Interpreter','latex');
zlabel('z [-]','Interpreter','latex');

for i = 1:N-1
    t(i) = t1 + (i-1)/(N-1)*(tN-t1);
    t(i+1) = t1 + i/(N-1)*(tN-t1);
    [~,y] = PBRFBP_traj(yopt_mlt(4*(i-1)+1:4*i),[t(i) t(i+1)],const);
    plot(y(:,1),y(:,2),'k')
    hold on;
    plot(yopt_mlt(4*(i-1)+1),yopt_mlt(4*(i-1)+2),'o','MarkerFaceColor','red')
end
plot(yopt_mlt(end-5),yopt_mlt(end-4),'o','MarkerFaceColor','red')
axis(axis+sign(axis)*0.2)
title('Rotating Earth-Moon frame @Earth-Moon baricenter, adimensional',Interpreter="latex")

%%%%%%%%%%%%%% N = 20 %%%%%%%%%%%%%%%
%prepare the initial state
N = 20;
t1 = ti;
tN = ti+delta;

tspan = [t1 tN];
[~,x] = sample_trajectory(tspan,x0,N,const);


%solve multiple shooting with gradients/jacobians
tic;
[yopt_mlt,fval_mlt] = multi_shoot(x,const,opts);
toc

%plot the results
figure;
t1 = yopt_mlt(end-1);
tN = yopt_mlt(end);
t = zeros(N,1);
hold on;
axis equal;
grid on;
xlabel('x [-]','Interpreter','latex');
ylabel('y [-]','Interpreter','latex');
zlabel('z [-]','Interpreter','latex');

for i = 1:N-1
    t(i) = t1 + (i-1)/(N-1)*(tN-t1);
    t(i+1) = t1 + i/(N-1)*(tN-t1);
    [~,y] = PBRFBP_traj(yopt_mlt(4*(i-1)+1:4*i),[t(i) t(i+1)],const);
    plot(y(:,1),y(:,2),'k')
    hold on;
    plot(yopt_mlt(4*(i-1)+1),yopt_mlt(4*(i-1)+2),'o','MarkerFaceColor','red')
end
plot(yopt_mlt(end-5),yopt_mlt(end-4),'o','MarkerFaceColor','red')
axis(axis+sign(axis)*0.2)
title('Rotating Earth-Moon frame @Earth-Moon baricenter, adimensional',Interpreter="latex")

%% Optimization of Function calls
clearvars; close all; clc;

%constants
mu = 1.21506683e-2;
const.mu = mu;

ms = 3.28900541e5; %Scaled mass of the Sun
const.ms = ms;

rho = 3.88811143e2; %Scaled Sun–(Earth + Moon) distance
const.rho = rho;

oms = -9.25195985e-1; %Scaled angular velocity of the Sun
const.oms = oms;

k = 1e3; % from [km]->[m]
Re = k*6378; %m Mean Earth radius
const.Re = Re;

Rm = k*1738; %m Mean Moon radius
const.Rm = Rm;

hi = k*167; %m Altitude of departure orbit
hf = k*100; %m Altitude of arrival orbit

DU = 3.84405000e8; %m Distance unit, from [m]->[adimensional]

ri = (Re + hi)/DU;
rf = (Rm + hf)/DU;
const.ri = ri;
const.rf = rf;

%data initialization
alfa = 1.5*pi;
beta = 1.41;
delta = 7;
ti = 0;
prm = [alfa, beta, delta, ti]';

s = param2state(prm,const);

%solve simple shooting without gradients/jacobians
opts = optimoptions(@fmincon,'Display','iter','Algorithm','active-set','MaxFunctionEvaluations',2000);
tic;
[yopt,fval] = shoot_optimized(s,const,opts);
toc

%plot initial guess
x0 = s(1:4);
tspan = [s(5) s(6)];
[~,y] = PBRFBP_traj(x0,tspan,const);
yMoon = [1-mu 0 0 0];
yEarth = [-mu 0 0 0];
figure;
plot(y(:,1),y(:,2),'r') %plot initial guess - comparison
hold on;

%plot results: rotating frame - origin baricenter earth moon
[~,y] = PBRFBP_traj(yopt(1:4),yopt(5:6),const);
plot(y(:,1),y(:,2),'k')


%plot astronomical objects
plot(yEarth(1),yEarth(2),'ok','MarkerFaceColor','red')
plot(yMoon(1),yMoon(2),'ok','MarkerFaceColor','blue')
legend({'Guess Trajectory','Solution with Gradients','Earth','Moon'},'Interpreter','latex');

%finish the plot
axis equal;
grid on;
legend({'Guess Trajectory','Solution','Earth','Moon'},'Interpreter','latex');
title('Rotating Earth-Moon frame @Earth-Moon baricenter, adimensional',Interpreter="latex")
xlabel('x [-]','Interpreter','latex');
ylabel('y [-]','Interpreter','latex');
zlabel('z [-]','Interpreter','latex');
axis(axis+sign(axis)*0.2)




%solve simple shooting with gradients/jacobians
opts = optimoptions(@fmincon,'Display','iter-detailed','Algorithm','active-set',...
    'MaxFunctionEvaluations',2000,'SpecifyConstraintGradient',true,'SpecifyObjectiveGradient',true);
tic;
[yopt_grd,fval_grd] = shoot_optimized(s,const,opts);
toc

%plot initial guess
x0 = s(1:4);
tspan = [s(5) s(6)];
[~,y] = PBRFBP_traj(x0,tspan,const);
yMoon = [1-mu 0 0 0];
yEarth = [-mu 0 0 0];
figure;
plot(y(:,1),y(:,2),'r') %plot initial guess - comparison
hold on;

%plot results: rotating frame - origin baricenter earth moon
[~,y] = PBRFBP_traj(yopt_grd(1:4),yopt_grd(5:6),const);
plot(y(:,1),y(:,2),'k')


%plot astronomical objects
plot(yEarth(1),yEarth(2),'ok','MarkerFaceColor','red')
plot(yMoon(1),yMoon(2),'ok','MarkerFaceColor','blue')
legend({'Guess Trajectory','Solution with Gradients','Earth','Moon'},'Interpreter','latex');

%finish the plot
axis equal;
grid on;
legend({'Guess Trajectory','Solution with gradients','Earth','Moon'},'Interpreter','latex');
title('Rotating Earth-Moon frame @Earth-Moon baricenter, adimensional',Interpreter="latex")
xlabel('x [-]','Interpreter','latex');
ylabel('y [-]','Interpreter','latex');
zlabel('z [-]','Interpreter','latex');
axis(axis+sign(axis)*0.2)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PERFORMANCE DATA
%                             WTO GRADIENT   | WITH GRADIENT
%PAR. COMP.                   4.700190/8       0.759241/8
%PAR. COMP. + SHARED STATE    16.788716/8      0.306691/8
% SHARED STATE                2.794537         0.316812
% NORMAL                      5.482328         0.705622 

%% End of File Functions
function s = param2state(prm,const)
%PARAM2STATE transformation from parameters to position, velocity state
%
% PROTOTYPE:
%   s = param2state(prm,const)
%
% INPUT:
%    prm[dim]          parameters alfa, beta, gamma, delta [unit]
%    const[str]        struct that contains all constants needed by the model
%
% OUTPUT:
%    s[dim]            state position and velocity [unit]
%
% CONTRIBUTORS:
%   Davide Zamblera - Politecnico di Milano - Space eng.
%

    r0 = const.ri;
    mu = const.mu;
    alfa = prm(1);
    beta = prm(2);
    delta = prm(3);
    ti = prm(4);
    v0 = beta*sqrt((1-mu)/r0);
    s = zeros(6,1);
    s(1) = r0*cos(alfa)-mu;
    s(2) = r0*sin(alfa);
    s(3) = -(v0-r0)*sin(alfa);
    s(4) = (v0-r0)*cos(alfa);
    s(5) = ti;
    s(6) = ti+delta;
end

function [x,y,xdot,ydot,ti,tf] = ext_st(s)
    x = s(1);
    y = s(2);
    xdot = s(3);
    ydot = s(4);
    if nargout > 4
        ti = s(5);
        tf = s(6);
    end

end

function [F,GF] = singleShootObj(s,const)
%SINGLESHOOTOBJ implementation of the shooting problem objective function
%
% PROTOTYPE:
%     [F,GF] = singleShootObj(s,const)
%
% INPUT:
%    s[dim]           variables of the optimization problem [unit]
%    const[str]        struct that contains all constants needed by the model
%
% OUTPUT:
%    F[dim]           objective function [unit]
%    GF[dim]          gradient of objective function [unit]
%
% CONTRIBUTORS:
%   Davide Zamblera - Politecnico di Milano - Space eng.
%

    mu = const.mu;
    sgn1 = 1;
    sgn2 = 1;
    
    [xi,yi,xdoti,ydoti,ti,tf] = ext_st(s);
    si = s(1:4);
    tspan = [ti tf];
    if nargout > 1
        [PHI, sf] = STMvrt(si,ti,tf,const);
    else
        sf = PBRFBP_prg(si,tspan,const);
    end
    [xf,yf,xdotf,ydotf] = ext_st(sf);
    
    ri = const.ri;
    rf = const.rf;
    Dvi = sqrt((xdoti-yi)^2 + (ydoti + xi + mu)^2) - sqrt((1-mu)/ri); %abs
    if Dvi < 0
        sgn1 = -1;
        Dvi = sgn1*Dvi;
    end
    Dvf = sqrt((xdotf-yf)^2 + (ydotf + xf + mu - 1)^2) - sqrt(mu/rf);
    if Dvf < 0
        sgn2 = -1;
        Dvf = sgn2*Dvf;
    end
    F = Dvi + Dvf;
    if nargout > 1
        GF = zeros(6,1);
      
        fi = PBRFBP(ti,si,const);
        ff = PBRFBP(tf,sf,const);

        dDV1dx1 = sgn1*1/sqrt((xdoti-yi)^2+(ydoti+xi+mu)^2)*[ydoti+xi+mu, yi-xdoti, xdoti-yi, ydoti+xi+mu]';
        dDV2dx2 = sgn2*1/sqrt((xdotf-yf)^2 + (ydotf+xf+mu-1)^2)*[ydotf+xf+mu-1, yf-xdotf, xdotf-yf ydotf+xf+mu-1]';
        dDV2dx1 = PHI'*dDV2dx2;
        dDV2dt1 = -dDV2dx2'*(PHI*fi);
        dDV2dt2 = dDV2dx2'*ff;
        GF(1:4) = dDV1dx1 + dDV2dx1;
        GF(5) = dDV2dt1;
        GF(6) = dDV2dt2;
    end
end

function [c,ceq,GC,GCeq] = singleShootCon(s,const)
%SINGLESHOOTCON implementation of the shooting problem constraint
%
% PROTOTYPE:
%     [c,ceq,GC,GCeq] = singleShootCon(y)
%
% INPUT:
%    s[dim]           variables of the optimization problem [unit]
%    const[str]        struct that contains all constants needed by the model
%
% OUTPUT:
%    c[dim]           inequality constraints [unit]
%    ceq[dim]         equality constraints [unit]
%    GC[dim]          jacobian of inequality constraints [unit]
%    GCeq[dim]        jacobian of equality constraints [unit]
%
% CONTRIBUTORS:
%   Davide Zamblera - Politecnico di Milano - Space eng.
%
     
    c = [];
    GC = [];
    GCeq = zeros(4,6);
    mu = const.mu;
    
    si = s(1:4);
    [xi,yi,xdoti,ydoti,ti,tf] = ext_st(s);
    tspan = [ti tf];

    if nargout > 2
        [PHI,sf] = STMvrt(si,ti,tf,const);
        [xf,yf,xdotf,ydotf] = ext_st(sf);
        fi = PBRFBP(ti,si,const);
        ff = PBRFBP(tf,sf,const);
        cidxi = [2*mu+2*xi 2*yi 0 0 ; xdoti ydoti mu+xi yi];
        cfdxf = [2*(mu+xf-1) 2*yf 0 0; xdotf ydotf mu+xf-1 yf];
        GCeq(1:2,1:4) = cidxi;
        GCeq(3:4,1:4) = cfdxf*PHI;
        GCeq(3:4,5) = - cfdxf*PHI*fi;
        GCeq(3:4,6) = cfdxf*ff;
        GCeq = GCeq';

    else
        sf = PBRFBP_prg(si,tspan,const);
        [xf,yf,xdotf,ydotf] = ext_st(sf);
    end
    ri = const.ri;
    rf = const.rf;
    ceq = zeros(4,1);
    ceq(1) = (xi + mu)^2 + yi^2 - ri^2;
    ceq(2) = (xi + mu)*(xdoti - yi) + yi*(ydoti + xi + mu);
    ceq(3) = (xf + mu - 1)^2 + yf^2 - rf^2;
    ceq(4) = (xf + mu - 1)*(xdotf - yf) + yf*(ydotf + xf + mu - 1);
end

function f = PBRFBP(t,s,const)
%PBRFBP planar bicircular restricted four body problem 
%
% PROTOTYPE:
%     f = PBRFBP(t,s,const)
%
% INPUT:
%    t[dim]            time [s]
%    s[dim]            state [km; km/s]
%    const[str]        struct that contains all constants needed by the model
%
% OUTPUT:
%    f[dim]           rhs of the system [-]
%
% CONTRIBUTORS:
%   Davide Zamblera - Politecnico di Milano - Space eng.
%
% VERSIONS
% 2022-10-04: First version
%

    mu = const.mu;
    rho = const.rho;
    oms = const.oms;
    ms = const.ms;
    [x,y,vx,vy] = ext_st(s);
    f = zeros(4,1);
    f(1) = vx;
    f(2) = vy;
    %r1 = ((x + mu)^2 + y^2)^(1/2);
    %r2 = ((x + mu - 1)^2 + y^2)^(1/2);

    %OM3 = 1/2*(x^2+y^2) + (1-mu)/r1 + mu/r2 + 1/2*mu*(1-mu);
    %r3 = ((x-rho*cos(oms*t))^2 + (y - rho*sin(oms*t))^2)^(1/2);
    %OM4 = OM3 + ms/r3 - ms/rho^2*(x*cos(oms*t)+y*sin(oms*t));
    OM4dx = x - (mu*(2*mu + 2*x - 2))/(2*((mu + x - 1)^2 + y^2)^(3/2)) - (ms*cos(oms*t))/rho^2 + ((2*mu + 2*x)*(mu - 1))/(2*((mu + x)^2 + y^2)^(3/2)) - (ms*(2*x - 2*rho*cos(oms*t)))/(2*((x - rho*cos(oms*t))^2 + (y - rho*sin(oms*t))^2)^(3/2));
    OM4dy = y - (ms*sin(oms*t))/rho^2 - (mu*y)/((mu + x - 1)^2 + y^2)^(3/2) - (ms*(2*y - 2*rho*sin(oms*t)))/(2*((x - rho*cos(oms*t))^2 + (y - rho*sin(oms*t))^2)^(3/2)) + (y*(mu - 1))/((mu + x)^2 + y^2)^(3/2);
    f(3) = 2*vy + OM4dx;
    f(4) = -2*vx + OM4dy;
    

end

function y = PBRFBP_prg(x0,tspan,const,opts)
%CRTBP_prg propagator of state with PBRFBP model
%
% PROTOTYPE:
%     y = PBRFBP_prg(x0,tspan,const,opts)
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
    if nargin < 4
        opts = odeset('RelTol',1e-12,'AbsTol',1e-12);
    end
    [~,yn] = ode113(@(t,x) PBRFBP(t,x,const),tspan,x0,opts);
    y = yn(end,:)';
end

function [tt,y] = PBRFBP_traj(x0,tspan,const,opts)
%PBRFBP_prg return trajectory by integration of x0 with CRTBP model
%
% PROTOTYPE:
%     [tt,y] = PBRFBP_traj(x0,tspan,const,opts)
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
    if nargin < 4
        opts = odeset('RelTol',1e-12,'AbsTol',1e-12);
    end
    [tt,y] = ode113(@(t,x) PBRFBP(t,x,const),tspan,x0,opts);
end

function SS = rot2in_ctr_earth(tt,ss,mu)
%ROT2IN_CTR_EARTH transform state from rotating frame to inertial earth centered
%
% PROTOTYPE:
%   SS = rot2in_ctr_earth(tt,ss,mu)
%
% INPUT:
%    tt[dim]           time of the states [unit]
%    ss[dim]           states in rotating baricenter centered frame[unit]
%    mu[dim]           system gravitational constant [km^3/s^2]
%
% OUTPUT:
%	 SS[dim] 	       states in inertial earth centered frame[unit]
%
% CONTRIBUTORS:
%   Davide Zamblera - Politecnico di Milano - Space eng.
%
    SS = zeros(size(ss));
    for i = 1:length(tt)
        t = tt(i);
        x = ss(i,1);
        y = ss(i,2);
        xdot = ss(i,3);
        ydot = ss(i,4);
        SS(i,1) = (x + mu)*cos(t) - y*sin(t);
        SS(i,2) = (x + mu)*sin(t) + y*cos(t);
        SS(i,3) = (xdot - y)*cos(t) - (ydot + x + mu)*sin(t);
        SS(i,4) = (xdot - y)*sin(t) + (ydot + x + mu)*cos(t);
    end
end

function [yopt,fval] = simple_shoot(x0,const,opts)
%SIMPLE_SHOOT shooting problem, optimization via fmincon
%
% PROTOTYPE:
%     [yopt,fval] = simple_shoot(x0,const,opts)
%
% INPUT:
%    x0[dim]           initial guess of optimization problem [unit]
%    const[str]        struct that contains all constants needed by the model
%    opts[dim]         option of fmincon
%
% OUTPUT:
%    yopt[dim]        state that minimizes objective function [unit]
%    fval[dim]         minimum of objective function [unit]
%
% CONTRIBUTORS:
%   Davide Zamblera - Politecnico di Milano - Space eng.
%
    if nargin < 3
        opts = [];
    end
    A = zeros(1,6);
    A(5) = 1;
    A(6) = -1;
    b = 0;

    lb(1:6,1) = -inf;
    lb(1) = -inf;
    lb(2) = -inf;
    lb(5) = 0;
    lb(6) = 0;

    ub(1:6,1) = inf;
    ub(1) = inf;
    ub(2) = inf;
    ub(5) = 23;
    ub(6) = 23;
    [yopt,fval] = fmincon(@(s) singleShootObj(s,const),x0,A,b,[],[],lb,ub,@(s) singleShootCon(s,const),opts);
end

function [PHI, x] = STMvrt(x0,t0,t,const,opts)
%STM return state transition matrix computed with variational approach
%
% PROTOTYPE:
%     [PHI, x] = STM(x0,t0,t,mu,om,opts)
%
% INPUT:
%    x0[dim]           initial state [unit]
%    t0[dim]           initial time [unit]
%    t[dim]            final time [unit]
%    mu[dim]           system gravitational constant [km^3/s^2]
%    om[dim]           system angular velocity [rad/s]
%    opts[dim]         options of the integrator scheme
%
% OUTPUT:
%    x[dim]            state at final time
%    STM[dim]          State Transition Matrix from t0 to t [unit]
%
% CONTRIBUTORS:
%   Davide Zamblera - Politecnico di Milano - Space eng.
%
    if nargin < 5
        opts = odeset('RelTol',1e-12,'AbsTol',1e-12);
    end
    y0 = [x0; reshape(eye(4),[16 1])];
    tspan = [t0 t];
    [~,yout] = ode113(@(t,x) stmPBRFBP(t,x,const),tspan,y0,opts);
    yout = yout(end,:);
    PHI = reshape(yout(5:end),[4 4]);
    x = yout(1:4)';
end

function [PHI, x] = STM(x0,t0,t,const,opts)
%STM return state transition matrix computed with finite difference
%
% PROTOTYPE:
%     [PHI, x] = STM(x0,t0,t,mu,om,opts)
%
% INPUT:
%    x0[dim]           initial state [unit]
%    t0[dim]           initial time [unit]
%    t[dim]            final time [unit]
%    mu[dim]           system gravitational constant [km^3/s^2]
%    om[dim]           system angular velocity [rad/s]
%    opts[dim]         options of the integrator scheme
%
% OUTPUT:
%    x[dim]            state at final time
%    STM[dim]          State Transition Matrix from t0 to t [unit]
%
% CONTRIBUTORS:
%   Davide Zamblera - Politecnico di Milano - Space eng.
%
    if nargin < 5
        opts = odeset('RelTol',1e-12,'AbsTol',1e-12);
    end
    ns = 4;
    e = zeros(ns);
    PHI = zeros(ns);
    tspan = [t0 t];
    phi0 = PBRFBP_prg(x0,tspan,const,opts);
    for k = 1:ns
        e(k,k) = sqrt(eps)*max(1, abs(x0(k)));
        phix = PBRFBP_prg(x0 + e(:,k),tspan,const,opts);
        PHI(:,k) = (phix - phi0)/e(k,k);

    end
    x = phi0;
end


function f = stmPBRFBP(t,s,const)
%STMPBRFBP planar bicircular restricted four body problem augmented by the
%stm variational equations
%
% PROTOTYPE:
%     f = stmPBRFBP(t,s,const)
%
% INPUT:
%    t[dim]            time [s]
%    s[dim]            state [km; km/s]
%    const[str]        struct that contains all constants needed by the model
%
% OUTPUT:
%    f[dim]           rhs of the system [-]
%
% CONTRIBUTORS:
%   Davide Zamblera - Politecnico di Milano - Space eng.
%
% VERSIONS
% 2022-10-04: First version
%

    mu = const.mu;
    rho = const.rho;
    oms = const.oms;
    ms = const.ms;
    [x,y,vx,vy] = ext_st(s);
    PHI = reshape(s(5:20),[4 4]);
    f = zeros(4,1);
    f(1) = vx;
    f(2) = vy;
    %r1 = ((x + mu)^2 + y^2)^(1/2);
    %r2 = ((x + mu - 1)^2 + y^2)^(1/2);

    %OM3 = 1/2*(x^2+y^2) + (1-mu)/r1 + mu/r2 + 1/2*mu*(1-mu);
    %r3 = ((x-rho*cos(oms*t))^2 + (y - rho*sin(oms*t))^2)^(1/2);
    %OM4 = OM3 + ms/r3 - ms/rho^2*(x*cos(oms*t)+y*sin(oms*t));
    OM4dx = x - (mu*(2*mu + 2*x - 2))/(2*((mu + x - 1)^2 + y^2)^(3/2)) - (ms*cos(oms*t))/rho^2 + ((2*mu + 2*x)*(mu - 1))/(2*((mu + x)^2 + y^2)^(3/2)) - (ms*(2*x - 2*rho*cos(oms*t)))/(2*((x - rho*cos(oms*t))^2 + (y - rho*sin(oms*t))^2)^(3/2));
    OM4dy = y - (ms*sin(oms*t))/rho^2 - (mu*y)/((mu + x - 1)^2 + y^2)^(3/2) - (ms*(2*y - 2*rho*sin(oms*t)))/(2*((x - rho*cos(oms*t))^2 + (y - rho*sin(oms*t))^2)^(3/2)) + (y*(mu - 1))/((mu + x)^2 + y^2)^(3/2);
    f(3) = 2*vy + OM4dx;
    f(4) = -2*vx + OM4dy;
    A = zeros(4);
    A(1:2,3:4) = eye(2);
    A(3,4) = 2;
    A(4,3) = -2;
    A(3,1) = (mu - 1)/((mu + x)^2 + y^2)^(3/2) - mu/((mu + x - 1)^2 + y^2)^(3/2) - ms/((x - rho*cos(oms*t))^2 + (y - rho*sin(oms*t))^2)^(3/2) + (3*ms*(2*x - 2*rho*cos(oms*t))^2)/(4*((x - rho*cos(oms*t))^2 + (y - rho*sin(oms*t))^2)^(5/2)) + (3*mu*(2*mu + 2*x - 2)^2)/(4*((mu + x - 1)^2 + y^2)^(5/2)) - (3*(2*mu + 2*x)^2*(mu - 1))/(4*((mu + x)^2 + y^2)^(5/2)) + 1;
    A(3,2) = (3*ms*(2*x - 2*rho*cos(oms*t))*(2*y - 2*rho*sin(oms*t)))/(4*((x - rho*cos(oms*t))^2 + (y - rho*sin(oms*t))^2)^(5/2)) + (3*mu*y*(2*mu + 2*x - 2))/(2*((mu + x - 1)^2 + y^2)^(5/2)) - (3*y*(2*mu + 2*x)*(mu - 1))/(2*((mu + x)^2 + y^2)^(5/2));
    A(4,1) = (3*ms*(2*x - 2*rho*cos(oms*t))*(2*y - 2*rho*sin(oms*t)))/(4*((x - rho*cos(oms*t))^2 + (y - rho*sin(oms*t))^2)^(5/2)) + (3*mu*y*(2*mu + 2*x - 2))/(2*((mu + x - 1)^2 + y^2)^(5/2)) - (3*y*(2*mu + 2*x)*(mu - 1))/(2*((mu + x)^2 + y^2)^(5/2));
    A(4,2) = (mu - 1)/((mu + x)^2 + y^2)^(3/2) - mu/((mu + x - 1)^2 + y^2)^(3/2) - ms/((x - rho*cos(oms*t))^2 + (y - rho*sin(oms*t))^2)^(3/2) + (3*ms*(2*y - 2*rho*sin(oms*t))^2)/(4*((x - rho*cos(oms*t))^2 + (y - rho*sin(oms*t))^2)^(5/2)) - (3*y^2*(mu - 1))/((mu + x)^2 + y^2)^(5/2) + (3*mu*y^2)/((mu + x - 1)^2 + y^2)^(5/2) + 1;
    f(5:20,1) = reshape(A*PHI,[16 1]);
end

function [yopt,fval] = multi_shoot(x0,const,opts)
%MULTI_SHOOT  multiple shooting problem, optimization via fmincon
%
% PROTOTYPE:
%     [yopt,fval] = multi_shoot(x0,const,opts)
%
% INPUT:
%    x0[dim]           initial guess of optimization problem [unit]
%    const[str]        struct that contains all constants needed by the model
%    opts[dim]         option of fmincon
%
% OUTPUT:
%    yopt[dim]        state that minimizes objective function [unit]
%    fval[dim]         minimum of objective function [unit]
%
% CONTRIBUTORS:
%   Davide Zamblera - Politecnico di Milano - Space eng.
%
    if nargin < 3
        opts = [];
    end
    
    m = length(x0);
    const.m = m;
    N = (m-2)/4;
    const.N = N;
    A = zeros(1,m);
    A(m-1) = 1;
    A(m) = -1;
    b = 0;

    lb = zeros(m,1);
    lb(1:end-2) = -inf;
    lb(end-1:end) = 0;
    ub = zeros(m,1);
    ub(1:end-2) = inf;
    ub(end-1) = 23;
    ub(end) = 23;
    [yopt,fval] = fmincon(@(s) multiShootObj(s,const),x0,A,b,[],[],lb,ub,@(s) multiShootCon(s,const),opts);

end


function [F,GF] = multiShootObj(s,const)
%MULTISHOOTOBJ implementation of the multiple shooting problem objective function
%
% PROTOTYPE:
%     [F,GF] = multiShootObj(s,const)
%
% INPUT:
%    y[dim]           variables of the optimization problem [unit]
%    const[str]        struct that contains all constants needed by the model
%
% OUTPUT:
%    F[dim]           objective function [unit]
%    GF[dim]          gradient of objective function [unit]
%
% CONTRIBUTORS:
%   Davide Zamblera - Politecnico di Milano - Space eng.
%
    sgn1 = 1;
    sgn2 = 1;
    mu = const.mu;
    N = const.N;
    si = s(1:4);
    [xi,yi,xdoti,ydoti] = ext_st(si);
    
    sf = s(end-5:end-2);
    [xf,yf,xdotf,ydotf] = ext_st(sf);
    
    
    ri = const.ri;
    rf = const.rf;
    Dvi = sqrt((xdoti-yi)^2 + (ydoti + xi + mu)^2) - sqrt((1-mu)/ri);
    if Dvi < 0
        sgn1 = -1;
        Dvi = sgn1*Dvi;
    end
    Dvf = sqrt((xdotf-yf)^2 + (ydotf + xf + mu - 1)^2) - sqrt(mu/rf);
    if Dvf < 0
        sgn2 = -1;
        Dvf = sgn2*Dvf;
    end
    F = Dvi + Dvf;
    if nargout > 1
        
        GF = zeros(4*N+2,1);

        dDV1dx1 = sgn1*1/sqrt((xdoti-yi)^2+(ydoti+xi+mu)^2)*[ydoti+xi+mu, yi-xdoti, xdoti-yi, ydoti+xi+mu]';
        dDV2dx2 = sgn2*1/sqrt((xdotf-yf)^2 + (ydotf+xf+mu-1)^2)*[ydotf+xf+mu-1, yf-xdotf, xdotf-yf ydotf+xf+mu-1]';
        GF(1:4) = dDV1dx1;
        GF(end-5:end-2) = dDV2dx2;
    end
end



function [c,ceq,GC,GCeq] = multiShootCon(s,const)
%MULTISHOOTCON implementation of the multiple shooting problem constraints
%
% PROTOTYPE:
%     [c,ceq,GC,GCeq] = multiShootCon(s,const)
%
% INPUT:
%    y[dim]           variables of the optimization problem [unit]
%    const[str]        struct that contains all constants needed by the model
%
% OUTPUT:
%    c[dim]           inequality constraints [unit]
%    ceq[dim]         equality constraints [unit]
%    GC[dim]          jacobian of inequality constraints [unit]
%    GCeq[dim]        jacobian of equality constraints [unit]
%
% CONTRIBUTORS:
%   Davide Zamblera - Politecnico di Milano - Space eng.
%
     
    mu = const.mu;
    N = const.N;
    Re = const.Re;
    Rm = const.Rm;
    DU = const.DU;
    
    si = s(1:4);
    [xi,yi,xdoti,ydoti] = ext_st(si);
    
    sf = s(end-5:end);
    [xf,yf,xdotf,ydotf,t1,tN] = ext_st(sf);
    
    ri = const.ri;
    rf = const.rf;
    
    ceq = zeros(4*N,1);
    c = zeros(2*N,1);
    
    
    if nargout > 2
        GCeq = spalloc(4*(N-1)+4,4*N+2,20*(N-1)+2*4*(N-1)+12);
        GC = spalloc(2*N,4*N+2,4*N);  
    end
    
    t = zeros(N,1);
    for j = 1:N
        t(j) = t1 + (j-1)/(N-1)*(tN-t1);
    end
    
    for j = 1:N-1
        idx = 1+4*(j-1):4*j;
        xj = s(idx);
        xjp1 = s(4*j+1:4*(j+1));
        
        if nargout > 2
            [PHI,sprg] = STMvrt(xj,t(j),t(j+1),const);
            ceq(idx) = sprg - xjp1;
            c(2*(j-1)+1) = (Re/DU)^2 - (xj(1)+mu)^2 - xj(2)^2;
            c(2*j) = (Rm/DU)^2 - (xj(1)+mu-1)^2 - xj(2)^2;
            GCeq(idx,idx) = PHI;
    
            f = PBRFBP(t(j),xj,const);
            fp1 = PBRFBP(t(j+1),sprg,const);
            GCeq(idx,end-1) = -PHI*f*((N-j)/(N-1)) + fp1*((N-j-1)/(N-1));
            GCeq(idx,end) = -PHI*f*((j-1)/(N-1)) + fp1*(j/(N-1));
        else
            sprg = PBRFBP_prg(xj,[t(j) t(j+1)],const);
            ceq(idx) = sprg - xjp1;
            c(2*(j-1)+1) = (Re/DU)^2 - (xj(1)+mu)^2 - xj(2)^2;
            c(2*j) = (Rm/DU)^2 - (xj(1)+mu-1)^2 - xj(2)^2;
        end
    end
    ceq(end-3) = (xi + mu)^2 + yi^2 - ri^2;
    ceq(end-2) = (xi + mu)*(xdoti - yi) + yi*(ydoti + xi + mu);
    ceq(end-1) = (xf + mu - 1)^2 + yf^2 - rf^2;
    ceq(end) = (xf + mu - 1)*(xdotf - yf) + yf*(ydotf + xf + mu - 1);
    
    
    c(end-1) = (Re/DU)^2 - (xf+mu)^2 - yf^2;
    c(end) = (Rm/DU)^2 - (xf+mu-1)^2 - yf^2;
    
    
    if nargout > 2
        I = eye(4);
        for j = 2:N
            GCeq(1+4*(j-2):4*(j-1),1+4*(j-1):4*j) = -I;
        end
        GCeq(end-3:end-2,1:4) = [2*(xi+mu) 2*yi 0 0; xdoti ydoti (xi+mu) yi];
        GCeq(end-1:end,end-5:end-2) = [2*(xf+mu-1) 2*yf 0 0; xdotf ydotf (xf+mu-1) yf];
        GCeq = GCeq';

        for j = 1:N
            xj = s(1+4*(j-1));
            yj = s(2+4*(j-1));
            GC(2*(j-1)+1:2*j,4*(j-1)+1:4*j-2) = [-2*(xj+mu), -2*yj; -2*(xj+mu-1) -2*yj];
        end
        GC = GC';
    end
    

end

function [tsample,xsample] = sample_trajectory(tspan,x0,N,const)
%SAMPLE_TRAJECTORY extract state of N points in tspan
%
% PROTOTYPE:
%     [c,ceq,GC,GCeq] = multiShootCon(s,const)
%
% INPUT:
%    tspan[dim]           interval of integration times  [unit]
%    x0[dim]              initial state  [unit]
%    N[dim]               number of points to sample [unit]
%    const[str]        struct that contains all constants needed by the model
%
% OUTPUT:
%    tsample[dim]           times of sampled points [unit]
%    xsample[dim]         vector of sampled state joined together [unit]
%
% CONTRIBUTORS:
%   Davide Zamblera - Politecnico di Milano - Space eng.
%
    t1 = tspan(1);
    tN = tspan(2);
    xsample = zeros(4*N+2,1);
    xsample(1:4) = x0;
    tsample = zeros(N,1);
    xsample(end-1) = t1;
    xsample(end) = tN;
    for j = 1:N
        tsample(j) = t1 + (j-1)/(N-1)*(tN-t1);
    end
    
    for j = 2:N
        xj = xsample(4*(j-2)+1:4*(j-1));
        xj = PBRFBP_prg(xj,[tsample(j-1) tsample(j)],const);
        xsample(4*(j-1)+1:4*j) = xj;
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%EXTRA FUNCTIONS - CODE OPTIMIZATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [sf, fi, ff, PHI] = computeall(x,const)
    
    si = x(1:4);
    ti = x(5);
    tf = x(6);
    tspan = [ti tf];
    if nargout > 1
        fi = PBRFBP(ti,si,const);
        
        [PHI, sf] = STMvrt(si,ti,tf,const);
        ff = PBRFBP(tf,sf,const);
    else
        sf = PBRFBP_prg(si,tspan,const);
    end
end


function [yopt,fval] = shoot_optimized(x0,const,opts)

if nargin < 3 % No options supplied
    opts = [];
end

xLast = []; % Last place computeall was called
sf = [];
PHI = [];
fi = [];
ff = [];

A = zeros(1,6);
A(5) = 1;
A(6) = -1;
b = 0;

lb(1:6,1) = -inf;
lb(1) = -inf;
lb(2) = -inf;
lb(5) = 0;
lb(6) = 0;

ub(1:6,1) = inf;
ub(1) = inf;
ub(2) = inf;
ub(5) = 23;
ub(6) = 23;

% Call fmincon
fun = @(s) objfun(s,const,opts);
cfun = @(s) constr(s,const,opts);
[yopt,fval] = fmincon(fun,x0,A,b,[],[],lb,ub,cfun,opts);

    function [F,GF] = objfun(s,const,opts)
        if ~isequal(s,xLast) % Check if computation is necessary
            if opts.SpecifyObjectiveGradient
                [sf, fi, ff, PHI] = computeall(s,const);
            else
                sf = computeall(s,const);
            end
            xLast = s;
        end
        % Now compute objective function
        mu = const.mu;
        sgn1 = 1;
        sgn2 = 1;
        
        [xi,yi,xdoti,ydoti] = ext_st(s);
        
        [xf,yf,xdotf,ydotf] = ext_st(sf);
        
        ri = const.ri;
        rf = const.rf;
        Dvi = sqrt((xdoti-yi)^2 + (ydoti + xi + mu)^2) - sqrt((1-mu)/ri); %abs
        if Dvi < 0
            sgn1 = -1;
            Dvi = sgn1*Dvi;
        end
        Dvf = sqrt((xdotf-yf)^2 + (ydotf + xf + mu - 1)^2) - sqrt(mu/rf);
        if Dvf < 0
            sgn2 = -1;
            Dvf = sgn2*Dvf;
        end
        F = Dvi + Dvf;
        if opts.SpecifyObjectiveGradient
            GF = zeros(6,1);
          
    
            dDV1dx1 = sgn1*1/sqrt((xdoti-yi)^2+(ydoti+xi+mu)^2)*[ydoti+xi+mu, yi-xdoti, xdoti-yi, ydoti+xi+mu]';
            dDV2dx2 = sgn2*1/sqrt((xdotf-yf)^2 + (ydotf+xf+mu-1)^2)*[ydotf+xf+mu-1, yf-xdotf, xdotf-yf ydotf+xf+mu-1]';
            dDV2dx1 = PHI'*dDV2dx2;
            dDV2dt1 = -dDV2dx2'*(PHI*fi);
            dDV2dt2 = dDV2dx2'*ff;
            GF(1:4) = dDV1dx1 + dDV2dx1;
            GF(5) = dDV2dt1;
            GF(6) = dDV2dt2;
        end
    end

    function [c,ceq,GC,GCeq] = constr(s,const,opts)
        if ~isequal(s,xLast) % Check if computation is necessary
            if opts.SpecifyConstraintGradient
                [sf, fi, ff, PHI] = computeall(s,const);
            else
                sf = computeall(s,const);
            end
            xLast = s;
        end
        % Now compute constraint function
        c = [];
        GC = [];
        GCeq = [];
        mu = const.mu;
        
        [xi,yi,xdoti,ydoti] = ext_st(s);
    
        if opts.SpecifyConstraintGradient
            GCeq = zeros(4,6);
            [xf,yf,xdotf,ydotf] = ext_st(sf);
            cidxi = [2*mu+2*xi 2*yi 0 0 ; xdoti ydoti mu+xi yi];
            cfdxf = [2*(mu+xf-1) 2*yf 0 0; xdotf ydotf mu+xf-1 yf];
            GCeq(1:2,1:4) = cidxi;
            GCeq(3:4,1:4) = cfdxf*PHI;
            %GCeq(1:2,5) = cidxi*fi;
            GCeq(3:4,5) = - cfdxf*PHI*fi;
            GCeq(3:4,6) = cfdxf*ff;
            GCeq = GCeq';
    
        else
            [xf,yf,xdotf,ydotf] = ext_st(sf);
        end
        ri = const.ri;
        rf = const.rf;
        ceq = zeros(4,1);
        ceq(1) = (xi + mu)^2 + yi^2 - ri^2;
        ceq(2) = (xi + mu)*(xdoti - yi) + yi*(ydoti + xi + mu);
        ceq(3) = (xf + mu - 1)^2 + yf^2 - rf^2;
        ceq(4) = (xf + mu - 1)*(xdotf - yf) + yf*(ydotf + xf + mu - 1);
    end

end