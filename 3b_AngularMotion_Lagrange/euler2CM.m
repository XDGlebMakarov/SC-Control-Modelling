function A = euler2CM(arr)
    A = rotmatrix(3,arr(3))*rotmatrix(1,arr(2))*rotmatrix(3,arr(1));
end