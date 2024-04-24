function B = ortGS(A)
    sz = size(A);
    B = zeros(sz);
    for i=1:sz(2)
        b = A(:,i);
        for j = 1:i-1
            b = b - proj(A(:,i),B(:,j));
        end
        B(:,i) = b/norm(b);
    end
end

function p = proj(a,b)
    p = (a'*b)/(b'*b)*b;
end