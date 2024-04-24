function [dobj,dw] = rk4right(obj,w,params)
    sz = size(obj);
    e_g = [0;0;-1];
    if prod(sz == [3 1])
        dobj = dynamics.dynEuler(obj,w);
        e_g_BF = angle2dcm(obj(1),obj(2),obj(3),'ZXZ')*e_g;
    elseif prod(sz == [3 3])
        dobj = dynamics.dynCosMatrix(obj,w);
        e_g_BF = obj*e_g;
    elseif prod(sz == [4 1])
        dobj = dynamics.dynQuat(obj,w);
        e_g_BF = quattrans(obj,e_g);
    else
        throw(MException('myComponent:inputError', ...
              'Ошибка размера терминов ориентации: [%d %d]',sz(1),sz(2)))
    end
    dw = params.J\(-cross(w,params.J*w) + params.m*params.g*params.L*...
                                                    cross([0;0;1],e_g_BF));
end