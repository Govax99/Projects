clc
clear
close all
%% Symmetric Data
Ix = 0.0504;
Iy = 0.0504;
Iz = 0.0109;
J = [Ix 0 0; 0 Iy 0; 0 0 Iz];
invJ = inv(J);
omegax0 = 0.45;
omegay0 = 0.52;
omegaz0 = 0.55;
omega0 = [omegax0, omegay0, omegaz0]';
omegar = 0;
Ir = 0;
lam = (Iz - Ix)/Ix * omegaz0;
omegaxAn = @(t) omegax0*cos(lam*t) - omegay0*sin(lam*t);
omegayAn = @(t) omegax0*sin(lam*t) + omegay0*cos(lam*t);
omegazAn = @(t) omegaz0;
%% after simulation
time = out.tout;
plot(time,omegaxAn(time));
hold on;
plot(time,omegayAn(time));
plot([time(1),time(end)],[omegaz0,omegaz0]);
title('Analytical Solution');
legend({'\omega_x','\omega_y','\omega_z'});
xlabel('time [s]');
ylabel('\omega [rad/s]')
hold off;
figure;
omega = out.simout(:,1:3);
hn = out.simout(:,4);
plot(time,omega(:,1),time,omega(:,2),time,omega(:,3));
legend({'\omega_x','\omega_y','\omega_z'});
title('Simulink Simulation');
xlabel('time [s]');
ylabel('\omega [rad/s]')

%% angular momentum and energy
T = zeros(length(time),1);
for i = 1:length(time)
    T(i)= 0.5* (Ix*omega(i,1)^2 + Iy*omega(i,2)^2 + Iz*omega(i,3)^2);
end
subplot(1,2,1);
plot(time,T);
title('Kinetic Energy');
subplot(1,2,2)
plot(time,hn);
title('Angular Momentum');