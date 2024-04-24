function R = rotmatrix(ax,alpha)
    switch ax
        case 1
            R = [1, 0         ,0         ;
                 0, cos(alpha),sin(alpha);
                 0,-sin(alpha),cos(alpha)];
        case 2
            R = [cos(alpha),0,-sin(alpha);
                 0         ,1, 0         ; 
                 sin(alpha),0, cos(alpha)];
        case 3
            R = [ cos(alpha),sin(alpha),0;
                 -sin(alpha),cos(alpha),0;
                  0         ,0         ,1];
    end
end