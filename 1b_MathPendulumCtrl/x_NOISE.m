function res = x_NOISE(x,xn)
    res = x + [normrnd(0,xn(1));
               normrnd(0,xn(2))];
end