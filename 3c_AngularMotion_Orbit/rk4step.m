function [x1,Q1,w1] = rk4step(x,Q,w,dt,params)
    [kx1,kQ1,kw1] = rk4right(x         ,Q         ,w         ,params);
    [kx2,kQ2,kw2] = rk4right(x+kx1*dt/2,Q+kQ1*dt/2,w+kw1*dt/2,params);
    [kx3,kQ3,kw3] = rk4right(x+kx2*dt/2,Q+kQ2*dt/2,w+kw2*dt/2,params);
    [kx4,kQ4,kw4] = rk4right(x+kx3*dt  ,Q+kQ3*dt  ,w+kw3*dt  ,params);
    x1 = x + (kx1+2*kx2+2*kx3+kx4)*dt/6;
    Q1 = Q + (kQ1+2*kQ2+2*kQ3+kQ4)*dt/6;
    Q1 = Q1/norm(Q1); % нормировка кватерниона
    w1 = w + (kw1+2*kw2+2*kw3+kw4)*dt/6;
end