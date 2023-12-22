function [a,en,i,OM,om,M,f] = car2kep(r,v,mu)
%car2kep transformation from cartesian coordinates to keplerian parameters
%
% PROTOTYPE:
%   [a,en,i,OM,om,M] = car2kep(r,v,mu)
%
% INPUT:
%    r[3x1]     position vector [Km]
%    v[3x1]     velocity vector [Km/s]
%    mu[1]      planetary constants of the planet (mu = mass * G) [km^3/s^2]
%
% OUTPUT:
%	 a[1] 	    semi-major axis [Km]
% 	 en[1]      Eccentricity (norm) [-]
%	 i[1] 	    inclination [rad]
%	 OM[1] 	    longitude of the ascending node [rad]
%	 om[1] 	    argument of periapsis [rad]
%	 M[1] 	    mean anomaly [rad]
%
% CONTRIBUTORS:
%   Davide Zamblera - Politecnico di Milano - Space eng.
%
% VERSIONS
% 2021-10-11: First version
%
if nargin < 3
    mu = astroConstants(13);
end
rn = norm(r);
vn = norm(v);
h = cross(r,v);
hn = norm(h);
a = -mu/(2*(0.5*vn^2-mu/rn));
vr = dot(v,r)/rn;
e = 1/mu*((vn^2-mu/rn)*r-rn*vr*v);
en = norm(e);
i = acos(h(3)/hn); %in i defined in (0,pi) no check
k = [0 0 1]';
if i ~= 0 && i~=pi
    N = cross(k,h);
    Nn = norm(N);
else
    N = [1 0 0]';
    Nn = norm(N);
end
OM = acos(N(1)/Nn);
if N(2) < 0
    OM = 2*pi - OM;
end
if abs(en) > 1e-10
    om = acos(dot(e,N)/(en*Nn));
    f = acos(dot(r,e)/(rn*en));
else
    e = [1 0 0]';
    om = 0;
    f = acos(dot(r,e)/(rn));
end
if e(3) < 0
    om = 2*pi - om;
elseif e(3) == 0
    if e(2) < 0
        om = 2*pi - om;
    end
end
if dot(v,r)<0
    f = 2*pi - f;
end
M = th2M(f,en);
end

