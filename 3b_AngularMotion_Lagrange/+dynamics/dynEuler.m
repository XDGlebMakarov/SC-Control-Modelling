function dangs = dynEuler(angs,w)
    dangs = zeros(3,1);
    theta = angs(2);
    phi = angs(3);
    if sin(theta) == 0
        dangs(1) = 0;
        dangs(3) = w(3);
    else
        dangs(1) = (w(1)*sin(phi) + w(2)*cos(phi))/sin(theta);
        dangs(3) = w(3) - dangs(1)*cos(theta);
    end
    dangs(2) = w(1)*cos(phi) - w(2)*sin(phi);
end