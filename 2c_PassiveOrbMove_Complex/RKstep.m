function x1 = RKstep(x,t,dt,params)
    jd = params.jd0 + t;
    k1 = RKright(x        , jd     , params);
    k2 = RKright(x+k1*dt/2, jd+dt/2, params);
    k3 = RKright(x+k2*dt/2, jd+dt/2, params);
    k4 = RKright(x+k3*dt  , jd+dt  , params);
    x1 = x + (k1+2*k2+2*k3+k4)*dt/6;
end