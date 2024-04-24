function res = quattrans(q,v)
    res = qmult(qmult(qconj(q),[0;v]),q);
    res = res(2:4);
end