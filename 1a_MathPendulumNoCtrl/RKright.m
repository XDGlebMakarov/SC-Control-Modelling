function dxdt = RKright(x,params)
    dxdt = [x(2);-params.g/params.l*sin(x(1))];
end