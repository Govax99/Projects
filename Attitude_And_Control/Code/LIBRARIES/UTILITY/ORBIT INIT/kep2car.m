function [r,v] = kep2car(a,en,i,OM,om,M,mu)
%kep2car transformation from keplerian parameterss to cartesian coordinates
%
% PROTOTYPE:
%   [a,en,i,OM,om,M] = kep2car(r,v,mu)
%
% INPUT:
%	 a[1] 	    semi-major axis [Km]
% 	 en[1]      eccentricity (norm) [-]
%	 i[1] 	    inclination [rad]
%	 OM[1] 	    longitude of the ascending node [rad]
%	 om[1] 	    argument of periapsis [rad]
%	 M[1] 	    mean anomaly [rad]
%    mu[1]      planetary constants of the planet (mu = mass * G) [km^3/s^2]
%
% OUTPUT:
%    r[3x1]     position vector [Km]
%    v[3x1]     velocity vector [Km/s]
%
% CONTRIBUTORS:
%   Davide Zamblera - Politecnico di Milano - Space eng.
%
% VERSIONS
% 2021-10-11: First version
%
if nargin < 7
    mu = astroConstants(13);
end
if M < pi
    E0 = real(M + en/2);
else
    E0 = real(M - en/2);
end
E = 0;
if M ~= 0
    E = fzero(@(E) E-en*sin(E)-M,E0);
end
f = 2*atan(sqrt((1+en)/(1-en))*tan(E/2));
hn = sqrt(mu*a*(1-en^2));
rn = (a*(1-en^2))/(1+en*cos(f));
r1 = [rn*cos(f),rn*sin(f),0]';
v1 = [-mu/hn*sin(f),mu/hn*(en+cos(f)),0]';


A1 = [cos(om) sin(om) 0;-sin(om) cos(om) 0;0 0 1];


A2 = [1 0 0;0 cos(i) sin(i);0 -sin(i) cos(i)];


A3 = [cos(OM) sin(OM) 0;-sin(OM) cos(OM) 0; 0 0 1];

A = (A1*A2*A3)';
r = A*r1;
v = A*v1;
r(isnan(r)) = 0;
end

