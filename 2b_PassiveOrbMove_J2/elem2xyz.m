function [r,v] = elem2xyz(i,e,p,Omega,omega,u,mu)
    nu = u-omega;
    R = p/(1+e*cos(nu));
    r_orb = R*[cos(nu);sin(nu);0];
    v_orb = sqrt(mu/p)*e*sin(nu)*[cos(nu);sin(nu);0] + sqrt(mu*p)/R*[-sin(nu);cos(nu);0];
    A = rotmat(3,omega)*rotmat(1,i)*rotmat(3,Omega);
    r = A'*r_orb;
    v = A'*v_orb;
end