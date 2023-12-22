function [g,h] = ghnorm(year)
%ghnorm transforms the data (g,h) from igrf.txt in Gaussian-Legendre norm to
%quasi-Schmidt normalization
%
% PROTOTYPE:
%   [g,h] = ghnorm(year)
%
% INPUT:
%    year[-]            year of simulation in range [1900-2025]
%
% OUTPUT:
%	 g[13x14] 	        first gaussian coefficent (n->rows ,m->columns) [T]
%	 h[13x14] 	        second gaussian coefficent (n->rows ,m->columns) [T]
%
% CONTRIBUTORS:
%   Davide Zamblera - Politecnico di Milano - Space eng.
%
% VERSIONS
% 2021-11-06: First version
%
years = {'1900' '1905' '1910' '1915' '1920' '1925' '1930' '1935' '1940' '1945' '1950' '1955' '1960' '1965' '1970' '1975' '1980' '1985' '1990' '1995'   '2000'    '2005'   '2010'    '2015'   '2020'};
I = find(strcmp(years, num2str(year)));
%extract data
file1 = fopen("igrf.txt");
C = textscan(file1,strcat('%s%f%f',repmat('%*f',1,I-1),'%f',repmat('%*f',1,25-I),'%f'),'MultipleDelimsAsOne',true, 'Delimiter',' [;', 'HeaderLines',4);
len = length(cell2mat(C(:,2)));
gh = C(:,1);
n = cell2mat(C(:,2));
m = cell2mat(C(:,3));
val = cell2mat(C(:,4));
sv = cell2mat(C(:,5));
N = 13;
%construct matrixes g,h
g=zeros(N,N+1);
h=zeros(N,N+1);
hsv=zeros(N,N+1);
gsv=zeros(N,N+1);
for x=1:len
    if strcmp(gh{1,1}{x,1},'g')
        g(n(x),m(x)+1) = val(x);
        gsv(n(x),m(x)+1) = sv(x);
    else
    h(n(x),m(x)+1) = val(x);
    hsv(n(x),m(x)+1) = sv(x);
    end
end
%using hsv gsv prediction for years (2021-2025)
delta = year - 2020;
if delta > 0 && delta < 6
    h = h + delta*hsv;
    g = g + delta*gsv;
end
%normalization transform
for x = 1:max(n)
    for y = 1:max(m)
        h(x,y) = S(x,y-1)*h(x,y);
        g(x,y) = S(x,y-1)*g(x,y);
    end
end
g = g*10^-9;
h = h*10^-9;
fclose(file1);
