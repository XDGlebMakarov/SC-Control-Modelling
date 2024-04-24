function [obj1,w1] = rk4step(obj,w,dt,params)
    [ko1,kw1] = rk4right(obj         ,w         ,params);
    [ko2,kw2] = rk4right(obj+ko1*dt/2,w+kw1*dt/2,params);
    [ko3,kw3] = rk4right(obj+ko2*dt/2,w+kw2*dt/2,params);
    [ko4,kw4] = rk4right(obj+ko3*dt  ,w+kw3*dt  ,params);
    obj1 = obj + (ko1+2*ko2+2*ko3+ko4)*dt/6;
    sz = size(obj1);
    if prod(sz == [3 3])
        obj1 = ortGS(obj1); % ортогонализация Грамма-Шмидта
    elseif prod(sz == [4 1])
        obj1 = obj1/norm(obj1); % нормировка кватерниона
    elseif prod(sz == [3 1])
        for k=1:3
            if abs(obj1(k)) > pi
                obj1(k) = obj1(k) - sign(obj1(k))*2*pi;
            end
        end
    else
        throw(MException('myComponent:inputError', ...
              'Ошибка размера терминов ориентации: [%d %d]',sz(1),sz(2)))
    end
    w1 = w + (kw1+2*kw2+2*kw3+kw4)*dt/6;
end