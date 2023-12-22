function M = th2M(th,e)
E = qck(2*atan(sqrt((1-e)/(1+e))*tan(th/2)));
M = qck(E - e*sin(E));
end

