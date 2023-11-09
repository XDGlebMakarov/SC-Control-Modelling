function dxdt = RKright(x,jd,params)
    r = x(1:3);
    GCRF2ITRF = frame_transformation.simpleGCRF2ITRF(jd);
    r_GCRF = GCRF2ITRF'*r;
    addacc = models.getGravNxN(r_GCRF, params.kc, params.ks, params.ncg, ...
                               params.R_Earth, params.mu, params.N_harm);
    nr = norm(r);
    dxdt = [x(4:6);
            -params.mu/nr^3*r + GCRF2ITRF*addacc];
end