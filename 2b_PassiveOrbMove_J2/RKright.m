function dxdt = RKright(x,params)
    r = x(1:3);
    nr = norm(r);
    dxdt = [x(4:6);
            -params.mu/nr^3*r + params.delta/nr^5*(r - 5*r(3)^2/nr^2*r + 2*r(3)*[0;0;1])];
end