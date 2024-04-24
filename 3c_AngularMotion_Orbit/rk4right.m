function [dx,dQ,dw] = rk4right(x,Q,w,params)
    r = x(1:3);
    nr = norm(r);
    dx = [x(4:6);
          -params.mu/(nr^3)*r];
    dQ = [-Q(2),-Q(3),-Q(4);
           Q(1),-Q(4), Q(3);
           Q(4), Q(1),-Q(2);
          -Q(3), Q(2), Q(1)]*w/2;
    r_BF = quattrans(Q,r);
    dw = params.J\(-cross(w,params.J*w) + 3*params.mu/(nr^5)*cross(r_BF,params.J*r_BF));
end