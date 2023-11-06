function x1 = RKstep(x,dt,params)
    k1 = RKright(x        ,params);
    k2 = RKright(x+k1*dt/2,params);
    k3 = RKright(x+k2*dt/2,params);
    k4 = RKright(x+k3*dt  ,params);
    x1 = x + (k1+2*k2+2*k3+k4)*dt/6;
end