function dA = dynCosMatrix(A,w)
    dA = -vec2skewmat(w)*A;
end