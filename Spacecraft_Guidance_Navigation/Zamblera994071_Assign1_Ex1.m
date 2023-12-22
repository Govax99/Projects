% Spacecraft Guidance and Navigation (2020/2021)
% Assignment # 1
% Author: Davide Zamblera

%% Section 1
clearvars; close all; clc;

%data initialization
mu = 3.0359e-6;

f =  @(x) (x - (1-mu)*(x+mu)./abs(x+mu).^(3) - mu*(x+mu-1)./abs(x+mu-1).^(3));

%find L2 point
xL2 = fbisection(f,[1.001,1.25]);
fprintf('The L2 point x coordinate is: %.10f\n',xL2);
fprintf('Value of the f(x) at L2: %g\n',f(xL2));


%plot the function, search of the zero
figure;

xx = linspace(1.001,1.25,1000);
subplot(1,3,1);
plot(xx,f(xx),'k')
hold on;
plot(xL2,0,'o')
grid on;
xlabel('x [-]',Interpreter="latex")
ylabel('f(x) [-]',Interpreter="latex")


xx = linspace(0.8,0.999,1000);
subplot(1,3,2);
plot(xx,f(xx),'k')
grid on;
xlabel('x [-]',Interpreter="latex")
ylabel('f(x) [-]',Interpreter="latex")


xx = linspace(-2,-0.21,1000);
subplot(1,3,3);
plot(xx,f(xx),'k')
grid on;
xlabel('x [-]',Interpreter="latex")
ylabel('f(x) [-]',Interpreter="latex")


%plot to show the L2 point on the function
figure;
xx = linspace(1.001,1.1,1000);
plot(xx,f(xx),'k')
hold on;
plot(xL2,0,'xr')
grid on;
xlabel('x [-]',Interpreter="latex")
ylabel('f(x) [-]',Interpreter="latex")

%% Section 2
clearvars; close all; clc;

%data initialization
mu = 3.0359e-6;
xL2 = 1.0100701875;
t0 = 0;
t = 2;
x0 = 1.008296144180133;
y0 = 0;
z0 = 0.001214294450297;
vx = 0;
vy = 0.010020975499502;
vz = 0;

x = [x0 y0 z0 vx vy vz]'; %initial state

om = [0 0 1]';
tspan = [t0, t];
opts = odeset('RelTol',2.5e-14,'AbsTol',2.5e-14,'Events',@myEventsFcn); %stop at y = 0

%plot initial guess
[~,y] = CRTBP_traj(x,tspan,mu,om,opts);
figure;
plot3(xL2,0,0,'o','MarkerFaceColor','blue');
hold on;
plot3(y(:,1),y(:,2),y(:,3),'r')
grid on;


%find the initial state for periodic orbit
tol = 1e-12;

xper = find_periodic(x,tol);
fprintf('The periodic orbit has initial state: [');
fprintf('%g, ',xper(1:end-1));
fprintf('%d]\n',xper(end))
opts = odeset('RelTol',2.5e-14,'AbsTol',2.5e-14);

%plot the results for some periods
t = 4;
tspan = [t0 t];
[tt,y] = CRTBP_traj(xper,tspan,mu,om,opts);

plot3(y(:,1),y(:,2),y(:,3),'k')
view(45,30)
axis equal;
grid on;
xlabel('X [-]',Interpreter="latex")
ylabel('Y [-]',Interpreter="latex")
zlabel('Z [-]',Interpreter="latex")
title('Rotating Sun-Earth frame @Sun-Earth baricenter, adimensional',Interpreter="latex")
legend('L2 Lagrange Point','Initial guess propagation','Periodic Halo Orbit',Interpreter="latex",Location="northwest")
%% Section 3 
clearvars; close all; clc;
tic;

%data initialization
mu = 3.0359e-6;
t0 = 0;
x0 = 1.008296144180133;
y0 = 0;
z0 = 0.001214294450297;
vx = 0;
vy = 0.010020975499502;
vz = 0;
z0f = 0.0046;
n = 10; %number of orbits in the family
zz0 = linspace(z0,z0f,n);
xx = zeros(6,n);
x = [x0 y0 z0 vx vy vz]';
t = 4;
om = [0 0 1]';
tspan = [t0 t];
tol = 1e-12;
opts = odeset('RelTol',2.5e-14,'AbsTol',2.5e-14,'Events',@myEventsFcn);


%numerical continuation: solve for a certain z0 and plot results
figure;
colormap(turbo(n))
colcustom = turbo(n);
for i = 1:length(zz0)
    x(3) = zz0(i);
    x = find_periodic(x,tol);
    [tt,y] = CRTBP_traj(x,tspan,mu,om);
    plot3(y(:,1),y(:,2),y(:,3),'Color', colcustom(i,:))
    hold on;
    grid on;
end

%graphical enhancement of the plot
view(45,30)
axis equal;
c = colorbar; 
caxis([0.001214294450297  0.0046]);
c.Label.String = '$z_0 [-]$';
c.Label.Interpreter = 'Latex';
c.Label.FontSize = 14;
c.Label.Rotation = 0;
c.Label.Position(1) = 3;

xlabel('X [-]',Interpreter="latex")
ylabel('Y [-]',Interpreter="latex")
zlabel('Z [-]',Interpreter="latex")
title('Rotating Sun-Earth frame @Sun-Earth baricenter, adimensional',Interpreter="latex")
toc;
%% End of File Functions


function [x,k] = fbisection(fun,x0)
%FBISECTION find zero of fun in x0 range via bisection method
%
% PROTOTYPE:
%     x = fbisection(fun,x0)
%
% INPUT:
%    fun[dim]           function handle [unit]
%    x0[2]              search range [unit]
%
% OUTPUT:
%    x[dim]             zero of the function [unit]
%    k[dim]             number of function evaluations [unit]
%
% CONTRIBUTORS:
%   Davide Zamblera - Politecnico di Milano - Space eng.
%

tol = 1e-10;

x1 = x0(1);
x2 = x0(2);
f1 = fun(x1);
f2 = fun(x2);
k = 2;
xt = x1;
xp = x2;
if f1*f2 > 0
    error('Initial search span must contain a sign variation.')
end


while abs(xt-xp)/abs(xt) >= tol
    xp = xt;
    xt = (x1+x2)/2;
    ft = fun(xt);
    k = k+1;
    if ft*f2 < 0
        x1 = xt;
    elseif ft*f2 > 0
        x2 = xt;
    else
        break;
    end

end
x = xt;
end





% model of the dynamics
function f = CRTBP(~,xx,mu,om)
%CRTPB circular restricted 3 body problem in 3 dimensions
%
% PROTOTYPE:
%     f = CRTBP(~,xx,mu)
%
% INPUT:
%    t[dim]            time [s]
%    xx[dim]           state [km; km/s]
%    mu[dim]           system gravitational constant [km^3/s^2]
%    om[dim]           system angular velocity [rad/s]
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

r = xx(1:3);
v = xx(4:6);
r1 = r + [mu 0 0]';
r2 = r + [mu-1 0 0]';
r1_3 = norm(r1)^3;
r2_3 = norm(r2)^3;
f(1:3,1) = v;
f(4:6,1) = -2*cross(om,v) - cross(om,cross(om,r)) - (1-mu)/r1_3*r1 - mu/r2_3*r2;
end


function f = stmCRTBP(~,s,mu,om)
%stmCRTPB extended circular restricted 3 body problem with STM variational integration
%
% PROTOTYPE:
%     f = stmCRTBP(~,s,mu,om)
%
% INPUT:
%    t[dim]            time [s]
%    s[dim]            state [km; km/s]
%    mu[dim]           system gravitational constant [km^3/s^2]
%    om[dim]           system angular velocity [rad/s]
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

f = zeros(42,1);
r = s(1:3);
x = r(1);
y = r(2);
z = r(3);
v = s(4:6);
omz = om(3);
PHI = reshape(s(7:42),[6 6]);
r1 = r + [mu 0 0]';
r2 = r + [mu-1 0 0]';
r1_3 = norm(r1)^3;
r2_3 = norm(r2)^3;
f(1:3,1) = v;
f(4:6,1) = -2*cross(om,v) - cross(om,cross(om,r)) - (1-mu)/r1_3*r1 - mu/r2_3*r2;
A = zeros(6);
A(1:3,4:6) = eye(3);
A(4,5) = 2*om(3);
A(5,4) = -2*om(3);
%made with matlab symbolic toolbox
A(4:6,1) = [(mu - 1)/((mu + x)^2 + y^2 + z^2)^(3/2) + omz^2 - mu/((mu + x - 1)^2 + y^2 + z^2)^(3/2) + (3*mu*(2*mu + 2*x - 2)*(mu + x - 1))/(2*((mu + x - 1)^2 + y^2 + z^2)^(5/2)) - (3*(2*mu + 2*x)*(mu + x)*(mu - 1))/(2*((mu + x)^2 + y^2 + z^2)^(5/2));...
            (3*mu*y*(2*mu + 2*x - 2))/(2*((mu + x - 1)^2 + y^2 + z^2)^(5/2)) - (3*y*(2*mu + 2*x)*(mu - 1))/(2*((mu + x)^2 + y^2 + z^2)^(5/2));...
            (3*mu*z*(2*mu + 2*x - 2))/(2*((mu + x - 1)^2 + y^2 + z^2)^(5/2)) - (3*z*(2*mu + 2*x)*(mu - 1))/(2*((mu + x)^2 + y^2 + z^2)^(5/2))];
A(4:6,2) = [(3*mu*y*(mu + x - 1))/((mu + x - 1)^2 + y^2 + z^2)^(5/2) - (3*y*(mu + x)*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(5/2);...
            (mu - 1)/((mu + x)^2 + y^2 + z^2)^(3/2) + omz^2 - mu/((mu + x - 1)^2 + y^2 + z^2)^(3/2) - (3*y^2*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(5/2) + (3*mu*y^2)/((mu + x - 1)^2 + y^2 + z^2)^(5/2);...
            (3*mu*y*z)/((mu + x - 1)^2 + y^2 + z^2)^(5/2) - (3*y*z*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(5/2)];
A(4:6,3) = [(3*mu*z*(mu + x - 1))/((mu + x - 1)^2 + y^2 + z^2)^(5/2) - (3*z*(mu + x)*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(5/2);...
            (3*mu*y*z)/((mu + x - 1)^2 + y^2 + z^2)^(5/2) - (3*y*z*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(5/2);...
            (mu - 1)/((mu + x)^2 + y^2 + z^2)^(3/2) - mu/((mu + x - 1)^2 + y^2 + z^2)^(3/2) - (3*z^2*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(5/2) + (3*mu*z^2)/((mu + x - 1)^2 + y^2 + z^2)^(5/2)];
dPHI = A*PHI;
f(7:42,1) = reshape(dPHI,[36 1]);
end

% state propagators and trajectory plotters
function y = CRTBP_prg(x0,tspan,mu,om,opts)
%CRTBP_prg propagator of state with CRTBP model
%
% PROTOTYPE:
%     y = CRTBP_prg(x0,tspan,mu,om,opts)
%
% INPUT:
%    x0[dim]           initial state [unit]
%    tspan[2]          integration time span [unit]
%    mu[dim]           system gravitational constant [km^3/s^2]
%    om[dim]           system angular velocity [rad/s]
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
    [~,y] = ode113(@(t,x) CRTBP(t,x,mu,om),tspan,x0,opts);
    y = y(end,:)';
end

function [tt,y] = CRTBP_traj(x0,tspan,mu,om,opts)
%CRTBP_prg return trajectory by integration of x0 with CRTBP model
%
% PROTOTYPE:
%     [tt,y] = CRTBP_traj(x0,tspan,mu,om,opts)
%
% INPUT:
%    x0[dim]           initial state [unit]
%    tspan[2]          integration time span [unit]
%    mu[dim]           system gravitational constant [km^3/s^2]
%    om[dim]           system angular velocity [rad/s]
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
    [tt,y] = ode113(@(t,x) CRTBP(t,x,mu,om),tspan,x0,opts);
end

function [value,isterminal,direction] = myEventsFcn(~,y)
    value = y(2);
    isterminal = 1;
    direction = -1;
end

%STM and state propagators
function [PHI, x] = STM(x0,t0,t,mu,om,opts)
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
    if nargin < 6
        opts = odeset('RelTol',1e-12,'AbsTol',1e-12);
    end
    ns = 6;
    e = zeros(ns);
    PHI = zeros(ns);
    tspan = [t0 t];
    phi0 = CRTBP_prg(x0,tspan,mu,om,opts);
    for k = 1:ns
        e(k,k) = sqrt(eps)*max(1, abs(x0(k)));
        phix = CRTBP_prg(x0 + e(:,k),tspan,mu,om,opts);
        PHI(:,k) = (phix - phi0)/e(k,k);

    end
    x = phi0;
end

function [PHI, x] = STMvrt(x0,t0,t,mu,om,opts)
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
    if nargin < 6
        opts = odeset('RelTol',1e-12,'AbsTol',1e-12);
    end
    y0 = [x0; reshape(eye(6),[36 1])];
    tspan = [t0 t];
    [~,yout] = ode113(@(t,x) stmCRTBP(t,x,mu,om),tspan,y0,opts);
    yout = yout(end,:);
    PHI = reshape(yout(7:end)',[6 6]);
    x = yout(1:6)';
end



function xper = find_periodic(x,tol)
%FIND_PERIODIC find a periodic orbit (identified by initial state) through STM iteration
%
% PROTOTYPE:
%     xper = find_periodic(x,tol)
%
% INPUT:
%    x[dim]           initial state guess [unit]
%    tol[dim]         tolerance of the iterative process [unit]
%
% OUTPUT:
%    xper[dim]          initial state of periodic orbit [unit]
%
% CONTRIBUTORS:
%   Davide Zamblera - Politecnico di Milano - Space eng.
%
mu = 3.0359e-6;
t0 = 0;
vfx = inf;
vfz = inf;
t = 4;
om = [0 0 1]';
opts = odeset('RelTol',1e-13,'AbsTol',1e-13,'Events',@myEventsFcn);

while abs(vfx)>tol || abs(vfz)>tol
    [PHI,y] = STM(x,t0,t,mu,om,opts);
    vfx = y(4);
    vfz = y(6);
    vold = [vfx; vfz];
    dx = - 1/(PHI(4,1)*PHI(6,5) - PHI(4,5)*PHI(6,1))*[PHI(6,5) -PHI(4,5); -PHI(6,1) PHI(4,1)]*vold;
    x(1) = x(1) + dx(1);
    x(5) = x(5) + dx(2);
end
xper = x;

end
