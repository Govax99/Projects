clear;
close all;
clc;

%% TEST P
k = 13;
[P,dP] = coeffP(k,deg2rad(60));
P = P(2:end,:);
dP = dP(2:end,:);
alfaG = 0;
a = 6371.2;
Re = 6378;
H = 500;
r = Re + H;
[g,h] = ghnorm(2005);
phi = deg2rad(-180:1:180);
theta = deg2rad(-90:1:90);
b = zeros(length(phi),length(theta));
for i = 1:length(phi)
    for j = 1:length(theta)
        b_I = MAGNETIC_FIELD(r,theta(j),phi(i),k,alfaG,g,h,a);
        b(i,j) = norm(b_I);
    end
end
slice = 20:2:52;
slice = slice * 10^-6;
b = b';
[X,Y] = meshgrid(rad2deg(phi),rad2deg(theta));
figure;
img = imread("earth.jpg");
imagesc([-180 180],[90 -90],img);
hold on;
contour(X,Y,b,slice);
colorbar;
set(gca,'ydir','normal');