function [P,dP] = coeffP(k,theta)
%calculate recursively the Gaussian quasi-normalized associated Legendre function of degree n and order m
%
% PROTOTYPE:
%   Pnm = P(n,m,theta)
%
% INPUT:
%    n[-]            degree of the function Pnm
%    m[-]            order of the function Pnm
%    theta[-]        co-latitude at which Pnm is calculated [rad]
%
% OUTPUT:
%	 Pnm[-] 	     gaussian-legendre function 
%
% CONTRIBUTORS:
%   Davide Zamblera - Politecnico di Milano - Space eng.
%
% VERSIONS
% 2021-11-06: First version
%
% ATT-> conditions to verify, they were not on the theory
P = zeros(k+1);
dP = zeros(k+1);
P11=1; P10=P11;
dP11=0; dP10=dP11;
P20=P10; dP20=dP10;
for m=0:k
    for n=1:k
        if m<=n
            % Calculate Legendre polynomials and derivatives recursively
            if n==m
                P2 = sin(theta)*P11;
                dP2 = sin(theta)*dP11 + cos(theta)*P11;
                P11=P2; P10=P11; P20=0;
                dP11=dP2; dP10=dP11; dP20=0;
            elseif n==1
                P2 = cos(theta)*P10;
                dP2 = cos(theta)*dP10 - sin(theta)*P10;
                P20=P10; P10=P2;
                dP20=dP10; dP10=dP2;
            else
                K = ((n-1)^2-m^2)/((2*n-1)*(2*n-3));
                P2 = cos(theta)*P10 - K*P20;
                dP2 = cos(theta)*dP10 - sin(theta)*P10 - K*dP20;
                P20=P10; P10=P2;
                dP20=dP10; dP10=dP2;
            end
            P(n+1,m+1)=P2;
            dP(n+1,m+1)=dP2;
        end
    end
end
P(1,1) = 1;