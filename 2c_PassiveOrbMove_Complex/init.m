params = struct('mu',398600.4416e9, ... % м^3/c^2, грав параметр
                'R_Earth',6371.302e3, ... % м, радиус Земли
                'data0',[2023 7 9 0 55 31], ... % дата начала интегрирования
                'jd0',0, ... % начальная юлианская дата
                'model','EGM2008', ... % модель грав поля, EGM2008/EGM96
                'N_harm',10, ... % число гармоник грав поля
                'kc',[], ... % коэффициенты при cos
                'ks',[], ... % коэффициенты при sin
                'ncg',[], ... % матрица модели грав поля
                'i',pi/6, ... % рад, наклонение
                'e',0.5, ... % эксцентриситет
                'p',11000e3, ... % м, фокальный параметр
                'Omega',pi/4, ... % рад, долгота восходящего узла
                'omega',pi/3, ... % рад, аргумент перицентра
                'nu0',pi/6); % начальная истинная аномалия
[params.kc, params.ks, params.ncg] = ...
            models.loadCoefsGravNxN(pwd, params.N_harm, params.model);
params.jd = time_transformation.date2JD(params.data0);
[r0, v0] = elem2xyz(params.i, params.e, params.p, params.Omega, ...
                    params.omega, params.nu0, params.mu);
x_START = [r0;  % радиус-вектор, м
           v0]; % скорость, м/с
T = 1e5; % с
dt = 0.5;