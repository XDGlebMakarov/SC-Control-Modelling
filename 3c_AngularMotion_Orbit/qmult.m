function q = qmult(q1,q2)
    q1s = q1(1);
    q1v = q1(2:4);
    q2s = q2(1);
    q2v = q2(2:4);
    q = [q1s*q2s - q1v'*q2v;
         cross(q1v,q2v) + q1s*q2v + q2s*q1v];
end