function dxdt = RKright(x,u,params)
    dxdt = [x(2);-params.g/params.l*sin(x(1))] + [0;u];
end