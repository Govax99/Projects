function b_I = MAGNETIC_FIELD(r,theta,phi,k,alfaG,g,h,a)
%calculate in ECI coordinates the Earth magnetic field b via a spherical
%harmonic model (quasi-Schmidt normalized formulas)
%
% PROTOTYPE:
%   b_I = MAGNETIC_FIELD(r,theta,phi,k,alfaG,g,h,R)
%
% INPUT:
%    r[-]            radius at which b is computed [Km]
%    theta[-]        latitude at which b is computed [rad]
%    phi[-]          longitude at which Pnm is calculated [rad]
%    k[-]            order of the model
%    alfaG[-]        initial longitude of observer meridian [rad]
%	 g[13x14] 	     first gaussian coefficent (n->rows ,m->columns) [T]
%	 h[13x14] 	     second gaussian coefficent (n->rows ,m->columns) [T]
%    a[-]            reference radius of the Earth [from paper 6371.2 Km]
%
% OUTPUT:
%	 b_I[3x1] 	     Earth magnetic field in ECI coordinates
%
% CONTRIBUTORS:
%   Davide Zamblera - Politecnico di Milano - Space eng.
%
% VERSIONS
% 2021-11-06: First version
%
b_r = 0;
b_th = 0;
b_phi = 0;
theta = pi/2 - theta;
[P,dP] = coeffP(k,theta);
P = P(2:end,:);
dP = dP(2:end,:);
for n = 1:k
    s1 = 0;
    s2 = 0;
    s3 = 0;
    for m = 0:n
        s1 = s1 + (g(n,m+1) * cos(m*phi) + h(n,m+1) * sin(m*phi)) * P(n,m+1);
        s2 = s2 + (g(n,m+1) * cos(m*phi) + h(n,m+1) * sin(m*phi)) * dP(n,m+1);
        s3 = s3 + m*(-g(n,m+1) * sin(m*phi) + h(n,m+1) * cos(m*phi)) * P(n,m+1);
    end
    b_r = b_r + (a/r)^(n+2) * (n+1) * s1;
    b_th = b_th + (a/r)^(n+2) * s2;
    b_phi = b_phi + (a/r)^(n+2) * s3;
end
b_th = - b_th;
b_phi = -1/sin(theta) * b_phi;

delta = pi/2 - theta;
alfa = phi + alfaG;

b1 = (b_r * cos(delta) + b_th * sin(delta))*cos(alfa) - b_phi * sin(alfa);
b2 = (b_r*cos(delta)+b_th*sin(delta))*sin(alfa)+b_phi*cos(alfa);
b3 = (b_r*sin(delta) - b_th*cos(delta));
b_I = [b1,b2,b3]';
end