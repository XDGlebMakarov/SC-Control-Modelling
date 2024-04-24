function res = quattrans(q,v)
    res = qmult(qmult([q(1);-q(2:4)],[0;v]),q);
    res = res(2:4);
end