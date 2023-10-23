function res = Energy(x,params)
    res = params.m*(params.l*x(2))^2/2 - params.m*params.l*params.g*cos(x(1));
end