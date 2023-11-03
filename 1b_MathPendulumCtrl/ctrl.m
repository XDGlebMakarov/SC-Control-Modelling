function u = ctrl(x,req,k1,k2,un,um,params)
    u = k1*(x(1)-req) - k2*x(2) + params.g*params.l*sin(req) + ...
            normrnd(0,un);
    if abs(u) > um
        u = u*um/abs(u);
    end
end