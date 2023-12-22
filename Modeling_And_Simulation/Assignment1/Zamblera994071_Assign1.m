% Modeling and Simulation of Aerospace Systems (2020/2021)
% Assignment # 1
% Author: Davide Zamblera

%% Ex 1
clearvars; close all; clc;

%0) data definition and search for a guess solution
f = @(x) sin(x) + x - 1;

delta = 9e-9;
fcs = 0.51097343;

xx = -10:0.00001:10;
figure;
plot(xx,f(xx),'k')
hold on;
plot([0 0],[0 f(0)],'b--');
plot([1 1],[0 f(1)],'r--');
plot(fcs,f(fcs),'o')
axis([-4 5 -6 4])
grid on;
legend({'$f(x)=sin(x) + x - 1$','','',''},Interpreter="latex",Location="best")
xlabel('$x [-]$',Interpreter='latex');
ylabel('$y [-]$',Interpreter='latex');

%from figure guess a=0, b=1
x0 = [0 1];
figure;
plot(xx,f(xx),'k')
hold on;
axis([fcs-delta fcs+delta -delta delta])
clear('xx')

%a) bisection algorithm

[x_bis, c_bis] = fbisection(f,x0);
ftime = @() fbisection(f,x0);
t_bis = timeit(ftime);
plot(x_bis,f(x_bis),'o','Color',"#0072BD",'MarkerSize',8,'LineWidth',2);
grid on;
fprintf('%+73s','Result          Time [s]          Fcn eval'); fprintf('\n');
fprintf('%-30s', 'Bisection method:'); fprintf('%15.12f  %.12f  %10.0f\n',[x_bis,t_bis,c_bis])

%b) secant algorithm

[x_sec, c_sec] = fsecant(f,x0);
ftime = @() fsecant(f,x0);
t_sec = timeit(ftime);
plot(x_sec,f(x_sec),'o','Color',"#D95319",'MarkerSize',8,'LineWidth',2);
fprintf('%-30s', 'Secant method:'); fprintf('%15.12f  %.12f  %10.0f\n',[x_sec,t_sec,c_sec])

%b) falsi algorithm

[x_fal, c_fal] = ffalsi(f,x0);
ftime = @() ffalsi(f,x0);
t_fal = timeit(ftime);
plot(x_fal,f(x_fal),'o','Color','#77AC30','MarkerSize',8,'LineWidth',2);
fprintf('%-30s', 'Falsi method:'); fprintf('%15.12f  %.12f  %10.0f\n',[x_fal,t_fal,c_fal])

legend({'','bisection','secant','falsi'},'Location','best','FontSize',10,Interpreter="latex")
xlabel('$x [-]$',Interpreter='latex');
ylabel('$y [-]$',Interpreter='latex');

%% Ex 2

clearvars; close all; clc;


%0) data definition and search for a guess solution
f = @(x) [x(1).^2 + x(2) - 5; x(2).^2 - x(1)];
f1 = @(x,y) x.^2 + y - 5;
f2 = @(x,y) y.^2 - x;
x1 = linspace(-1,3,100);
x2 = linspace(-2.5,3,100);
[X,Y] = meshgrid(x1,x2);
figure;
contour(X,Y,f1(X,Y),[0 0],'r')
hold on;
contour(X,Y,f2(X,Y),[0 0],'b')
grid on;
xlabel('$x_1 [-]$',Interpreter='latex');
ylabel('$x_2 [-]$',Interpreter='latex');
x1 = [2,1.5]';
x2 = [2,-2]';

%1a) solution of the system via analitically derived jacobian

x_jacob.p1 = fnewton(f,x1,'jacob');


x_jacob.p2 = fnewton(f,x2,'jacob');

%1b) solution of the system via finite difference jacobian

x_finite.p1 = fnewton(f,x1,'finite');


x_finite.p2 = fnewton(f,x2,'finite');


%2) Plot and Display of Results
plot(x_jacob.p1(1),x_jacob.p1(2),'ok')
plot(x_jacob.p2(1),x_jacob.p2(2),'ok')
plot(x1(1),x1(2),'x','Color','#A2142F')
plot(x1(1),x2(2),'x','Color','#A2142F')
legend({'$f_1(x)=x_1^2 + x_2 -5$','$f_2(x)=x_2^2 - x_1$','Solutions','','Initial guesses',''},Interpreter="latex",Location="best")

fprintf('%+60s','Result               Finite Difference     Analytical'); fprintf('\n');
fprintf('%-30s', 'First Point [x]:'); fprintf('%15.16f   %.16f\n',[x_finite.p1(1),x_jacob.p1(1)])
fprintf('%-30s', 'First Point [y]:'); fprintf('%15.16f   %.16f\n',[x_finite.p1(2),x_jacob.p1(2)])
fprintf('%-30s', 'Second Point [x]:'); fprintf('%15.16f   %.16f\n',[x_finite.p2(1),x_jacob.p2(1)])
fprintf('%-30s', 'Second Point [y]:'); fprintf('%15.16f  %.16f\n',[x_finite.p2(2),x_jacob.p2(2)])
fprintf('%-40s', 'Accuracy First Zero -Analytical:'); fprintf('%15.16f   %.16f\n',f([x_jacob.p1(1),x_jacob.p1(2)])')
fprintf('%-40s', 'Accuracy First Zero -Finite:'); fprintf('%15.16f   %.16f\n',f([x_finite.p1(1),x_finite.p1(2)])')
fprintf('%-40s', 'Accuracy Second Zero -Analytical:'); fprintf('%15.16f   %.16f\n',f([x_jacob.p2(1),x_jacob.p2(2)])')
fprintf('%-40s', 'Accuracy Second Zero -Finite:'); fprintf('%15.16f   %.16f\n',f([x_finite.p2(1),x_finite.p2(2)])')

%% Ex 3

clearvars; close all; clc;

%1) data definition
ode = @(t,x) - (t.^2 - 1)./(t.^2 + 1).*x;
x0 = 1;

x = @(t) exp(2*atan(t) - t);

tspan = [0 10];
hh = [0.5, 0.2, 0.05, 0.01];

err_heun = zeros(1,4);
time_heun = zeros(1,4);

%2) computation of error and time for RK2
figure;
hold on;
grid on;

for i = 1:length(hh)
    f = @() RK2(ode,tspan,x0,hh(i));
    [tout,yout] = RK2(ode,tspan,x0,hh(i));
    time_heun(i) = timeit(f);
    yexact = x(tout);
    err = abs(yout - yexact)./(abs(yexact));
    err_heun(i) = max(err);
    plot(tout,err);
end
line([2 2], [1e-8 1], 'color','r','linestyle','--');
xlabel('t [s]',Interpreter="latex")
ylabel('$\frac{||x_{true} - x_{num}||}{||x_{true}||}$','Interpreter','latex',...
    'FontSize',20)
set(gca, 'YScale', 'log')
legend({'h = 0.5','h = 0.2', 'h = 0.05', 'h = 0.01'},Interpreter="latex",Location="best")

err_rk4 = zeros(1,4);
time_rk4 = zeros(1,4);


%3) computation of error and time for RK4
figure;
hold on;
grid on;
ylabel('x')
xlabel('t')
for i = 1:length(hh)
    f = @() RK4(ode,tspan,x0,hh(i));
    [tout,yout] = RK4(ode,tspan,x0,hh(i));
    time_rk4(i) = timeit(f);
    yexact = x(tout);
    err = abs(yout - yexact)./(abs(yexact));
    err_rk4(i) = max(err);
    plot(tout,err);

end
line([2 2], [1e-14 1], 'color','r','linestyle','--');
xlabel('t [s]',Interpreter="latex")
ylabel('$\frac{||x_{true} - x_{num}||}{||x_{true}||}$','Interpreter','latex',...
    'FontSize',20)
set(gca, 'YScale', 'log')
legend({'h = 0.5','h = 0.2', 'h = 0.05', 'h = 0.01'},Interpreter="latex",Location="best")


%4) Comparison between the two methods
figure;
hold on;
grid on;
plot(hh,err_heun);
plot(hh,err_rk4);
legend('RK2','RK4',Interpreter="latex",Location="best")
xlabel('h [-]',Interpreter="latex")
ylabel('$\frac{||x_{true} - x_{num}||}{||x_{true}||}$','Interpreter','latex',...
    'FontSize',20)
set(gca, 'YScale', 'log', 'XScale', 'log')
figure;
hold on;
grid on;
plot(hh,time_heun);
plot(hh,time_rk4);
xlabel('h [-]',Interpreter="latex")
ylabel('Computational Time t [s]',Interpreter="latex")
legend('RK2','RK4',Interpreter="latex",Location="best")
set(gca, 'YScale', 'log')


%% EX 4

clearvars; close all; clc;

%1) definition of the operators
I = eye(2);
FRK2 = @(A,h) I + A*h + (A*h)^2/2;
FRK4 = @(A,h) I+A*h + 0.5*(A*h)^2 + 1/6*(A*h)^3 + 1/24*(A*h)^4;


%2) solution for alfa = pi
alfa = pi;

hpi_RK2 = alfa_stb(FRK2,alfa);
hpi_RK4 = alfa_stb(FRK4,alfa);

fprintf('h marginal at alfa = pi for RK2: %.12f\n',hpi_RK2);
fprintf('h marginal at alfa = pi for RK4: %.12f\n',hpi_RK4);

%3) plot of region of stability and eigenvalues
figure;
ode_stab(FRK2)
hold on;
hlim = [2 3];
ode_stab(FRK4,hlim)
hh = [0.5, 0.2, 0.05, 0.01];
plot(hh,zeros(1,4),'xr')
legend({'RK2','RK4','Eigenvalues, t=0'},Interpreter="latex",Location="best")


%% EX 5

clearvars; close all; clc;

%1) Definition of limits and update operators
tol = [1e-3, 1e-4, 1e-5, 1e-6];
hlim(:,:,1) = [1e-4,1.9e-2; 1e-5,1.9e-3; 1e-6,1.9e-4; 1e-7,1.9e-5];
hlim(:,:,2) = [1e-3,1; 1e-3,1e-1; 1e-3,1e-1; 1e-4,1e-2];
hlim(:,:,3) = [1e-3,1; 1e-3,1; 1e-3,1; 1e-3,1];
I = eye(2);
FRK1 = @(A,h) I + A*h;
FRK2 = @(A,h) I + A*h + (A*h)^2/2;
FRK4 = @(A,h) I + A*h + (A*h)^2/2 + (A*h)^3/6 + (A*h)^4/24;
order = [1 2 4];
solver = {FRK1,FRK2,FRK4};
ktot = zeros(3,length(tol));
hpi = zeros(3,length(tol));

%2) Plot Region of Accuracy
tic;
for i = 1:3
    figure(i);
    for j = 1:length(tol)
    
        ode_acc(solver{i},tol(j),hlim(j,:,i));
        hold on;

        [hpi(i,j),ktot(i,j)] = alfa_acc(solver{i},pi,tol(j),hlim(j,:,i),order(i));
    end
    legend({'tol = 1e-3','tol = 1e-4','tol = 1e-5','tol = 1e-6'},Interpreter="latex",Location="best")
end
T = toc;

figure(1);
curves = gca().Children;
axes('position',[.65 .2 .25 .25])
box on
plot(curves(1).XData,curves(1).YData,Color='#7E2F8E');
hold on;
plot(curves(2).XData,curves(2).YData,Color='#EDB120');
axis tight

%3) Plot of the Function Evaluations vs tolerance
figure;
tol = flip(tol);
ktot = flip(ktot,2);
plot(tol,ktot);
grid on;
xlabel('tol [-]',Interpreter="latex")
ylabel('function evaluations [-]',Interpreter="latex")
legend('RK1','RK2','RK4',Interpreter="latex",Location="best")
set(gca, 'XScale', 'log', 'YScale','log')

%% Ex 6

clearvars; close all; clc;

%1) Definition of Operators
I = eye(2);
BI2th = @(A,h,th) (I - A*h*(1-th) + (A*h*(1-th))^2/2)\(I + A*h*th + (A*h*th)^2/2);
th = [0.1 0.3 0.4 0.7 0.9];
hlim = [0.001, 11];

%2) Plot Region of Stability
figure;
for i = 1:length(th)
    BI2 = @(A,h) BI2th(A,h,th(i));
    if i > 3
        ode_stab(BI2,hlim)
    else
        ode_stab(BI2,hlim,'--')
    end
    hold on;
end
axis([-12 12 -8 8])
legend({'$\theta=0.1$','$\theta=0.3$', '$\theta=0.4$', '$\theta=0.7$', '$\theta=0.9$'},Interpreter="latex",Location="best")

%% Ex 7

clearvars; close all; clc;

%0) Prepare data
B = [-180.5, 219.5; 179.5, -220.5];
tspan = [0,5];
x0 = [1 1]';

ode = @(t,x) B*x;

h = 0.1;
tout = tspan(1):h:tspan(2);
yexact = zeros(length(x0),length(tout));


%1) Compute relative errors of methods
for i = 1:length(tout)
    yexact(:,i) = expm(B*tout(i))*x0;
end
yexact = yexact';


% a) RK4 and h=0.1
[~,yout] = RK4(ode,tspan,x0,h);
err = vecnorm(yout - yexact,2,2)./(vecnorm(yexact,2,2));
figure;
yyaxis right
plot(tout,err);
set(gca, 'YScale', 'log')
hold on;
grid on;

% a) IEX4 and h=0.1
[~,yout] = IEX4(ode,tspan,x0,h);
err = vecnorm(yout - yexact,2,2)./(vecnorm(yexact,2,2));
yyaxis left
plot(tout,err);

xlabel('t [s]',Interpreter="latex")
ylabel('$\frac{||x_{true} - x_{num}||}{||x_{true}||}$','Interpreter','latex',...
    'FontSize',20)
set(gca, 'YScale', 'log')
legend({'IEX4 h=0.1','RK4 h=0.1'},Interpreter="latex",Location="best")
set(gca,'ycolor','k')


%2) Plot eigenvalues and regions of stability
lam = eig(B);
hlam = h*lam;
I = eye(2);
FBE = @(A,h) inv(I-A*h);
FRK4 = @(A,h) I+A*h + 0.5*(A*h)^2 + 1/6*(A*h)^3 + 1/24*(A*h)^4;
FIEX4 = @(A,h) -1/6*FBE(A,h) +4*FBE(A,h/2)*FBE(A,h/2) + ...
    -27/2*FBE(A,h/3)*FBE(A,h/3)*FBE(A,h/3)+32/3*FBE(A,h/4)*FBE(A,h/4)*FBE(A,h/4)*FBE(A,h/4);
figure;

hlim = [2 3];
ode_stab(FRK4,hlim)
hold on
hlim = [0 20];
ode_stab(FIEX4,hlim,'--')
plot(hlam(1),0,'x','Color','k')
plot(hlam(2),0,'x','Color','#77AC30')
axis([-50 20 -10 10])
legend({'RK4','IEX4','$h\lambda_1$','$h\lambda_2$'},Interpreter="latex",Location='eastoutside')
%% Ex 8

clearvars; close all; clc;

%Determine accuracy of the methods
ode = @(t,x) [-5/2*(1+8*sin(t))*x(1); (1-x(1))*x(2)+x(1)];
y0 = [1 1]';
tspan = [0 3];
h = 0.1;
tt = tspan(1):h:tspan(2);
opts = odeset('RelTol',2.5e-14,'AbsTol',2.5e-14);
[tex,yex] = ode113(ode,tt,y0,opts);


figure;
yyaxis left
hold on;
xlabel('t [s]',Interpreter="latex")
ylabel('$\frac{||x_{true} - x_{num}||}{||x_{true}||}$','Interpreter','latex',...
    'FontSize',20)
set(gca,'ycolor','k')
set(gca, 'YScale', 'log')
solvers = {@AM3, @ABM3, @BDF3};
for i = 1:length(solvers)
    solver = solvers{i};
    [tout,yout] = solver(ode,tspan,y0,h);
    err = vecnorm(yout - yex,2,2)./(vecnorm(yex,2,2));
    plot(tout,err);
end
yyaxis right;
[tout,yout] = AB3(ode,tspan,y0,h);
err = vecnorm(yout - yex,2,2)./(vecnorm(yex,2,2));
plot(tout,err);
legend('AM3','ABM3','BDF3','AB3',Interpreter="latex",Location="best")
grid on;
set(gca, 'YScale', 'log')

%Eigenvalue during simulation time

lam1 = @(t) 0.1*(-5/2*(1+8*sin(t)));
lam2 = @(t) 0.1*(1-exp(-20)*exp(-5/2*t + 20*cos(t)));

%Plotting of region of stability
O = zeros(2);
I = eye(2);
FAB3 = @(A,h) [O, I, O; O O I; 5/12*A*h, -4/3*A*h, (I + 23/12*A*h)];
FAM3 = @(A,h) [O, I; -(I - 5/12*A*h)\(1/12*A*h) (I-5/12*A*h)\(I+2/3*A*h)];
FABM3 = @(A,h) [O I O; O O I; 25/144*(A*h)^2, -(1/12*A*h+5/9*(A*h)^2) (I+13/12*A*h+115/144*(A*h)^2)];
FBDF3 = @(A,h) [O I O; O O I ; 2/11*inv(I-6/11*h*A), -9/11*inv(I-6/11*h*A), 18/11*inv(I-6/11*h*A)];
figure;
ode_stab(FAB3,[0.5 5])
hold on;
ode_stab(FAM3,[0.5 8])
ode_stab(FABM3,[1 5])
ode_stab(FBDF3,[0.00001,10],'--')
tt = 0:0.1:3;
ylam = zeros(1,length(tt));
plot3(lam1(tt),ylam,tt)
plot3(lam2(tt),ylam,tt)
view(45,45)
zlabel('t [s]',Interpreter='latex')
legend({'AB3','AM3','ABM3','BDF3','$h\lambda_1$','$h\lambda_2$'},Interpreter="latex",Location="best")

%Second plot highlight time behaviour
figure;
zz = 0:-0.04:-0.16;
hmarg = [-0.545454, -6, -1.72496, 6.66667];
for i=1:3
    plot([0 hmarg(i)],[zz(i) zz(i)],'LineWidth',1.5)
    hold on;
end
plot([0 hmarg(4)],[zz(4) zz(4)],'--','LineWidth',1.5)
plot(lam1(tt),tt)
plot(lam2(tt),tt)
axis([-3 2 -0.5 3])
axis equal
grid on
xlabel('$Re(h\lambda)$',Interpreter='latex')
ylabel('t [s]',Interpreter='latex')
legend({'AB3','AM3','ABM3','BDF3','$h\lambda_1$','$h\lambda_2$'},Interpreter="latex",Location="best")
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

tol = 1e-8;

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


function [x,k] = fsecant(fun,x0)
%FSECANT find zero of fun with x0 guesses via secant method
%
% PROTOTYPE:
%     x = fsecant(fun,x0)
%
% INPUT:
%    fun[dim]           function handle [unit]
%    x0[2]              first 2 guesses [unit]
%
% OUTPUT:
%    x[dim]             zero of the function [unit]
%    k[dim]             number of function evaluations [unit]
%
% CONTRIBUTORS:
%   Davide Zamblera - Politecnico di Milano - Space eng.
%

tol = 1e-8;
x1 = x0(1);
x2 = x0(2);
f1 = fun(x1);
f2 = fun(x2);
k = 2;

while abs(x2-x1)/abs(x2) >= tol
    xt = x2 - (x2-x1)*f2/(f2-f1);
    x1 = x2;
    f1 = f2;
    x2 = xt;
    f2 = fun(x2);
    k = k+1;

end
x = x2;
end

function [x,k] = ffalsi(fun,x0)
%FFALSI find zero of fun in x0 range via falsi method
%
% PROTOTYPE:
%     x = ffalsi(fun,x0)
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

tol = 1e-8;
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
    xt = (x1*f2 - x2*f1)/(f2 - f1);
    ft = fun(xt);
    k = k+1;
    if ft*f2 < 0
        x1 = xt;
        f1 = ft;
    elseif ft*f2 > 0
        x2 = xt;
        f2 = ft;
    else
        break;
    end
end
x = xt;
end



function [x,it] = fnewton(fun,x0,mode)
%FNEWTON find zero of a system of equation via Newton's method
%
% PROTOTYPE:
%     x = fsecant(fun,x0)
%
% INPUT:
%    fun[dim]           function handle [unit]
%    x0[2]              initial guess [unit]
%    mode[chr]          if 'jacob' compute with analytical jacobian
%                       if 'finite' compute with forward finite differences
%
% OUTPUT:
%    x[dim]             zero of the system [unit]
%
% CONTRIBUTORS:
%   Davide Zamblera - Politecnico di Milano - Space eng.
%

tol = 1e-8;

if strcmp(mode,'jacob')
    dfun = @(x) [2*x(1), 1; -1, 2*x(2)];
elseif strcmp(mode,'finite')
    dfun = @(x) jacob(fun,x);
end

k = 2;
it = zeros(length(x0),50);
it(:,1) = x0;
x = x0;
xn = x + 1;
while abs(norm(xn-x)) > tol
    x = xn;
    xn = x - dfun(x)\fun(x);
    it(:,k) = xn;
    k = k+1;
end
it(:,all(it == 0))=[];
x = xn;
end

function J = jacob(F,x0)
%JACOB compute the jacobian of vector function F(x)
%
% PROTOTYPE:
%   J = jacob(F,x0)
%
% INPUT:
%    F[dim]            function  [-]
%    x0[dim]           point of evaluation [-]
%
% OUTPUT:
%	 J[dim] 	       jacobian of function F at x0 [-]
%
% CONTRIBUTORS:
%   Davide Zamblera - Politecnico di Milano - Space eng.
%


ns = length(x0);
e = zeros(ns);
J = zeros(ns);
for k=1:ns
    e(k,k) = sqrt(eps)*max(1, abs(x0(k)));
    J(:,k) = (F(x0+e(:,k)) - F(x0))/e(k,k);
end
end



function [tout,yout] = RK2(ode,tspan,y0,h)
%RK2 Heun explicit method (equal blend of FE and PC)
%
% PROTOTYPE:
%     [tout,yout] = RK2(ode,tspan,y0,h)
%
% INPUT:
%    ode[dim]           rhs right hand side [unit]
%    tspan[dim]         limits of integration [unit]
%    y0[dim]            initial state [unit]
%    h[dim]             time step [unit]
%
% OUTPUT:
%    tout[dim]           vector of times [unit]
%    yout[dim]           vector of state at multiple times [unit]
%
% CONTRIBUTORS:
%   Davide Zamblera - Politecnico di Milano - Space eng.
%
 

%check input
if isrow(y0)
    y0 = y0';
end

f0 = ode(tspan(1),y0);
if isrow(f0) && ~isscalar(f0)
    error('Function should return column array.')
end

tout = tspan(1):h:tspan(2);
n = length(tout);
yout = zeros(length(y0),n);


yout(:,1) = y0;
for i = 2:n
    yp = yout(:,i-1);

    k1 = ode(tout(i-1),yp)*h;
    k2 = ode(tout(i-1)+h,yp+k1)*h;
    yout(:,i) = yp + 1/2*(k1 + k2);
end

tout = tout';
yout = yout';
end


function [tout,yout] = RK4(ode,tspan,y0,h)
%RK4 Runge Kutta 4th order
%
% PROTOTYPE:
%     [tout,yout] = RK4(ode,tspan,y0,h)
%
% INPUT:
%    ode[dim]           rhs right hand side [unit]
%    tspan[dim]         limits of integration [unit]
%    y0[dim]            initial state [unit]
%    h[dim]             time step [unit]
%
% OUTPUT:
%    tout[dim]           vector of times [unit]
%    yout[dim]           vector of state at multiple times [unit]
%
% CONTRIBUTORS:
%   Davide Zamblera - Politecnico di Milano - Space eng.
%
 

%check input
if isrow(y0)
    y0 = y0';
end

f0 = ode(tspan(1),y0);
if isrow(f0) && ~isscalar(f0)
    error('Function should return column array.')
end

tout = tspan(1):h:tspan(2);
n = length(tout);
yout = zeros(length(y0),n);


yout(:,1) = y0;
for i = 2:n
    yp = yout(:,i-1);

    k1 = ode(tout(i-1),yp)*h;
    k2 = ode(tout(i-1)+0.5*h,yp + 1/2*k1)*h;
    k3 = ode(tout(i-1)+0.5*h,yp + 1/2*k2)*h;
    k4 = ode(tout(i-1)+h,yp + k3)*h;
    yout(:,i) = yp + 1/6*(k1 + 2*k2 + 2*k3 + k4);
end

tout = tout';
yout = yout';
end

function hmax = alfa_stb(F,alfa,hlim)
%ALFA_STB return largest h for which stability of algorithm F is
% maintained along the direction alfa
%
% PROTOTYPE:
%     hmax = alfa_stb(F,alfa,hlim)
%
% INPUT:
%    F[dim]                recursive matrix function of A,h [unit]
%    alfa[dim]             angle from plus axis where stability is tested [rad]
%    hlim[dim]             minimum and maximum h step to consider in search for hmax [unit]
%
% OUTPUT:
%    hmax[dim]             maximum step size h which guarantees stability [unit]
%
% CONTRIBUTORS:
%   Davide Zamblera - Politecnico di Milano - Space eng.
%

if nargin < 3
    hlim = [0.001 3];
end

x = cos(alfa);
A = [0 1; -1 2*x];
f = @(h) max(abs(eig(F(A,h))))-1;
hmax = fzero(f,hlim);
end




function ode_stab(F,hlim,lineopt)
%ODE_STAB display the stability region of a ode method given F(z)
%
% PROTOTYPE:
%     ode_stab(F,hlim)
%
% INPUT:
%    F[dim]                recursive matrix function of A,h [unit]
%    hlim[dim]             minimum and maximum h step to consider in search for hmax [unit]
%    lineopt[chr]          select display line style
%
% OUTPUT:
%
% CONTRIBUTORS:
%   Davide Zamblera - Politecnico di Milano - Space eng.
%

if nargin < 2
    hlim = [0.001 3];
end
alfa = 0:0.01:pi;
hh_x = zeros(size(alfa));
hh_y = hh_x;
for i = 1:length(alfa)
    try
        
        hmax = alfa_stb(F,alfa(i),hlim);
    catch
        hmax = 0;
    end
    hh_x(i) = hmax*cos(alfa(i));
    hh_y(i) = hmax*sin(alfa(i));

end

hh_x = [hh_x , flip(hh_x)];
hh_y = [hh_y , flip(-hh_y)];

if nargin < 3
    plot(hh_x,hh_y)
else
    plot(hh_x,hh_y,lineopt)
end
title('Stability Region of the Ode Solver',Interpreter='latex');
ylabel('$Im\left(h\lambda \right)$',Interpreter='latex');
xlabel('$Re\left(h\lambda \right)$',Interpreter='latex');
axis equal;
grid on;
axis([-4 2 -3 3])

end


function [hmax,k] = alfa_acc(F,alfa,tol,hlim,order)
%ALFA_ACC return largest h for which accuracy of algorithm F is
% achieved along the direction alfa
%
% PROTOTYPE:
%     [hmax,k] = alfa_acc(F,alfa,tol,hlim,order)
%
% INPUT:
%    F[dim]                update operator function of A,h [unit]
%    alfa[dim]             angle from plus axis where stability is tested [rad]
%    tol[dim]              tolerance on final infinite norm which must be achieved [unit]
%    hlim[dim]             minimum and maximum h step to consider in search for hmax [unit]
%    order[dim]            order of the method, necessary to compute k [unit] 
%
% OUTPUT:
%    hmax[dim]             maximum step size h which guarantees a certain accuracy [unit]
%    k[dim]                number of function evaluations for a method F
%
% CONTRIBUTORS:
%   Davide Zamblera - Politecnico di Milano - Space eng.
%

f = @(h) test_prb(F,h,alfa,tol);
hmax = fzero(f,hlim);

if nargout > 1
    k = ceil(order/hmax);
    
end
end

function z = test_prb(F,h,alfa,tol)
    %child function of ALFA_ACC
    tspan = [0 1];
    m = floor(1/h); %number of step must be integer
    A = [0 1; -1 2*cos(alfa)];
    x0 = [1 1]';
    x_ex = expm(A*tspan(2))*x0;
    x_num = (F(A,h))^m*x0;
    x_num = (F(A,1-h*m))*x_num; %doing a last step to arrive to tspan(2)
    z = max(abs(x_ex - x_num)) - tol;

end

function ode_acc(F,tol,hlim)
%ODE_ACC display the accuracy region of a ode method on the model problem
%
% PROTOTYPE:
%     ode_acc(F,tol,hlim)
%
% INPUT:
%    F[dim]                update operator function of A,h [unit]
%    tol[dim]              accuracy tolerance to be achieved [unit]
%    hlim[dim]             minimum and maximum h step to consider in search for hmax [unit]
%
% OUTPUT:
%
% CONTRIBUTORS:
%   Davide Zamblera - Politecnico di Milano - Space eng.
%

if nargin < 3
    hlim = [0.01, 1.3];
end



alfa = 0:0.01:2*pi;
hh_x = zeros(size(alfa));
hh_y = hh_x;
for i = 1:length(alfa)
    try
        
        hmax = alfa_acc(F,alfa(i),tol,hlim);
    catch
        hmax = 0;
    end
    hh_x(i) = hmax*cos(alfa(i));
    hh_y(i) = hmax*sin(alfa(i));

end


plot(hh_x,hh_y)
title('Accuracy Region of the Ode Solver',Interpreter='latex');
ylabel('$Im\left(h\lambda \right)$',Interpreter='latex');
xlabel('$Re\left(h\lambda \right)$',Interpreter='latex');
axis equal;
grid on;
end



function [tout,yout] = IEX4(ode,tspan,y0,h)
%IEX4 Implicit Extrapolation 4th order
%
% PROTOTYPE:
%     [tout,yout] = IEX4(ode,tspan,y0,h)
%
% INPUT:
%    ode[dim]           rhs right hand side [unit]
%    tspan[dim]         limits of integration [unit]
%    y0[dim]            initial state [unit]
%    h[dim]             time step [unit]
%
% OUTPUT:
%    tout[dim]           vector of times [unit]
%    yout[dim]           vector of state at multiple times [unit]
%
% CONTRIBUTORS:
%   Davide Zamblera - Politecnico di Milano - Space eng.
% 

%check input
if isrow(y0)
    y0 = y0';
end

f0 = ode(tspan(1),y0);
if isrow(f0) && ~isscalar(f0)
    error('Function should return column array.')
end
options = optimset('Display','off');

tout = tspan(1):h:tspan(2);
n = length(tout);
yout = zeros(length(y0),n);
pred = zeros(length(y0),4);
alfa = [-1/6 4 -27/2 32/3];

yout(:,1) = y0;
for i = 2:n
    
    for j = 1:4
        yp = yout(:,i-1);
        tp = tout(i-1);
        g = h/j;
        for k = 1:j
            tp = tp + g;
            f = @(y) y - yp - ode(tp,y)*g;
            
            yp = fsolve(f,yp,options);
        end
        pred(:,j) = yp;
    end
    for j = 1:4
        yout(:,i) = yout(:,i) + alfa(j)*pred(:,j);
    end
end

tout = tout';
yout = yout';
end


function [tout,yout] = AB3(ode,tspan,y0,h)
%AB3 Adam Bashforth 3rd order
%
% PROTOTYPE:
%     [tout,yout] = AB3(ode,tspan,y0,h)
%
% INPUT:
%    ode[dim]           rhs right hand side [unit]
%    tspan[dim]         limits of integration [unit]
%    y0[dim]            initial state [unit]
%    h[dim]             time step [unit]
%
% OUTPUT:
%    tout[dim]           vector of times [unit]
%    yout[dim]           vector of state at multiple times [unit]
%
% CONTRIBUTORS:
%   Davide Zamblera - Politecnico di Milano - Space eng.
% 

%check input
if isrow(y0)
    y0 = y0';
end

f0 = ode(tspan(1),y0);
if isrow(f0) && ~isscalar(f0)
    error('Function should return column array.')
end


tout = tspan(1):h:tspan(2);
n = length(tout);
yout = zeros(length(y0),n);


yout(:,1) = y0;

% initial build up
fmin2 = ode(tout(1),yout(:,1));
yout(:,2) = y0 + h*fmin2;
fmin1 = ode(tout(2),yout(:,2));
yout(:,3) = yout(:,2) + h*fmin1;

for i = 4:n
    f = ode(tout(i-1),yout(:,i-1));
    yout(:,i) = yout(:,i-1) + h/12*(23*f-16*fmin1+5*fmin2);
    fmin2 = fmin1;
    fmin1 = f;
    
end

tout = tout';
yout = yout';
end


function [tout,yout] = AM3(ode,tspan,y0,h)
%AM3 Adam Moulton 3rd order
%
% PROTOTYPE:
%     [tout,yout] = AM3(ode,tspan,y0,h)
%
% INPUT:
%    ode[dim]           rhs right hand side [unit]
%    tspan[dim]         limits of integration [unit]
%    y0[dim]            initial state [unit]
%    h[dim]             time step [unit]
%
% OUTPUT:
%    tout[dim]           vector of times [unit]
%    yout[dim]           vector of state at multiple times [unit]
%
% CONTRIBUTORS:
%   Davide Zamblera - Politecnico di Milano - Space eng.
% 

%check input
if isrow(y0)
    y0 = y0';
end

f0 = ode(tspan(1),y0);
if isrow(f0) && ~isscalar(f0)
    error('Function should return column array.')
end
opts = optimset('Diagnostics','off', 'Display','off');

tout = tspan(1):h:tspan(2);
n = length(tout);
yout = zeros(length(y0),n);


yout(:,1) = y0;

% initial build up
fmin1 = ode(tout(:,1),yout(:,1));
impl = @(y) y - y0 - h*ode(tout(2),y);
yout(:,2) = fsolve(impl,y0,opts);
f = ode(tout(2),yout(:,2));

for i = 3:n
    impl  = @(y) y - yout(:,i-1) - h/12*(5*ode(tout(i),y) + 8*f - fmin1);
    yout(:,i) = fsolve(impl,yout(:,i-1),opts);
    fmin1 = f;
    f = ode(tout(i),yout(:,i));
    
end

tout = tout';
yout = yout';
end


function [tout,yout] = ABM3(ode,tspan,y0,h)
%ABM3 Adam Bashforth Moulton 3rd order
%
% PROTOTYPE:
%     [tout,yout] = AB3(ode,tspan,y0,h)
%
% INPUT:
%    ode[dim]           rhs right hand side [unit]
%    tspan[dim]         limits of integration [unit]
%    y0[dim]            initial state [unit]
%    h[dim]             time step [unit]
%
% OUTPUT:
%    tout[dim]           vector of times [unit]
%    yout[dim]           vector of state at multiple times [unit]
%
% CONTRIBUTORS:
%   Davide Zamblera - Politecnico di Milano - Space eng.
%

%check input
if isrow(y0)
    y0 = y0';
end

f0 = ode(tspan(1),y0);
if isrow(f0) && ~isscalar(f0)
    error('Function should return column array.')
end

tout = tspan(1):h:tspan(2);
n = length(tout);
yout = zeros(length(y0),n);


yout(:,1) = y0;

% initial build up
fmin2 = ode(tout(:,1),yout(:,1));
ypred = y0 + h*fmin2;
yout(:,2) = y0 + h*ode(tout(2),ypred);

fmin1 = ode(tout(2),yout(:,2));
ypred = yout(:,2) + h/2*(3*fmin1 - fmin2);
yout(:,3) = yout(:,2) + h/2*(ode(tout(3),ypred) + fmin1);

for i = 4:n
    f = ode(tout(i-1),yout(:,i-1));
    ypred = yout(:,i-1) + h/12*(23*f - 16*fmin1 + 5*fmin2);
    yout(:,i) =  yout(:,i-1) + h/12*(5*ode(tout(i),ypred) + 8*f - fmin1);
    fmin2 = fmin1;
    fmin1 = f;

    
end

tout = tout';
yout = yout';
end


function [tout,yout] = BDF3(ode,tspan,y0,h)
%BDF3 Backward Difference Formula 3rd order
%
% PROTOTYPE:
%     [tout,yout] = BDF3(ode,tspan,y0,h)
%
% INPUT:
%    ode[dim]           rhs right hand side [unit]
%    tspan[dim]         limits of integration [unit]
%    y0[dim]            initial state [unit]
%    h[dim]             time step [unit]
%
% OUTPUT:
%    tout[dim]           vector of times [unit]
%    yout[dim]           vector of state at multiple times [unit]
%
% CONTRIBUTORS:
%   Davide Zamblera - Politecnico di Milano - Space eng.
% 

%check input
if isrow(y0)
    y0 = y0';
end

f0 = ode(tspan(1),y0);
if isrow(f0) && ~isscalar(f0)
    error('Function should return column array.')
end
opts = optimset('Diagnostics','off', 'Display','off');

tout = tspan(1):h:tspan(2);
n = length(tout);
yout = zeros(length(y0),n);


yout(:,1) = y0;

% initial build up
ymin2 = y0;
impl = @(y) y - y0 - h*ode(tout(2),y);
yout(:,2) = fsolve(impl,y0,opts);
ymin1 = yout(:,2);
impl = @(y) y - 4/3*ymin1 + 1/3*ymin2 - 2/3*h*ode(tout(3),y);
yout(:,3) = fsolve(impl,ymin1,opts);
y = yout(:,3);

for i = 4:n
    impl  = @(ypl1) ypl1 - 18/11*y + 9/11*ymin1 - 2/11*ymin2 - 6/11*h*ode(tout(i),ypl1);
    yout(:,i) = fsolve(impl,yout(:,i-1),opts);
    ymin2 = ymin1;
    ymin1 = y;
    y = yout(:,i);
    
end

tout = tout';
yout = yout';
end

