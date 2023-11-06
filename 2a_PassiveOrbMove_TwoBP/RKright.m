function dxdt = RKright(x,params)
    dxdt = [x(4:6);
            -params.mu/norm(x(1:3))^3*x(1:3)];
end