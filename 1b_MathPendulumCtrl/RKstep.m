function x1 = RKstep(x,u,dt,params)
    k1 = RKright(x        ,u,params);
    k2 = RKright(x+k1*dt/2,u,params);
    k3 = RKright(x+k2*dt/2,u,params);
    k4 = RKright(x+k3*dt  ,u,params);
    x1 = x + (k1+2*k2+2*k3+k4)*dt/6;
end