function [i,e,p,Omega,omega,u] = xyz2elem(r,v,mu)
    c = cross(r,v);
    i = acos(c'*[0;0;1]/norm(c)); % наклонение
    l_ = cross([0;0;1],c); l = l_/norm(l_);
    Omega = atan2(l(2),l(1)); % аргумент восходящего узла
    f = cross(v,c) - mu*r/norm(r);
    n_ = cross(c,l); n = n_/norm(n_);
    omega = atan2(f'*n,f'*l); % аргумент перицентра
    u = atan2(r'*n,r'*l); % аргумент широты
    p = c'*c/mu; % фокальный параметр
    e = norm(f)/mu; % эксцентриситет
end