function [dobj,dw] = rk4right(obj,w,params)
    sz = size(obj);
    if prod(sz == [3 1])
        dobj = dynamics.dynEuler(obj,w);
    elseif prod(sz == [3 3])
        dobj = dynamics.dynCosMatrix(obj,w);
    elseif prod(sz == [4 1])
        dobj = dynamics.dynQuat(obj,w);
    else
        throw(MException('myComponent:inputError', ...
              'Ошибка размера терминов ориентации: [%d %d]',sz(1),sz(2)))
    end
    dw = params.J\(-cross(w,params.J*w));
end