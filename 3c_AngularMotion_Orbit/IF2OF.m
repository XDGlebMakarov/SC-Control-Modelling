function res = IF2OF(Omega,i,u)
    res = rotmat(3,Omega)*rotmat(1,i)*rotmat(3,u);
end