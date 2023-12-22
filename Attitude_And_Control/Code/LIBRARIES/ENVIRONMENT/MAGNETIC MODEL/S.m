function Snm = S(n,m)
%calculate recursively the function S of passage from Gaussian-Legendre to Schmidt
%
% PROTOTYPE:
%   Snm = S(n,m)
%
% INPUT:
%    n[-]            degree of the function Pnm
%    m[-]            order of the function Pnm
%
% OUTPUT:
%	 Snm[-] 	     transform normalization factor 
%
% CONTRIBUTORS:
%   Davide Zamblera - Politecnico di Milano - Space eng.
%
% VERSIONS
% 2021-11-06: First version
%
% ATT-> conditions to verify, they were not on the theory
if n == 1 && m == 0
    Snm = 1;
    return;
end
if m == 0
    Snm = S(n-1,0)*(2*n - 1)/n;
    return;
end
Snm = S(n,m-1)*sqrt(((m == 1) + 1)*(n-m+1)/(n+m));
return;
