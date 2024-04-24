params = struct('J',diag([0.15,0.15,0.10]), ... % тензор инерции, кг*м^2
                'mu',398600.4416, ... % км^3/c^2, грав параметр
                'i',63.4*pi/180, ... % рад, наклонение
                'e',0, ... % эксцентриситет
                'p',6600, ... % км, фокальный параметр
                'Omega',pi/4, ... % рад, долгота восходящего узла
                'omega',pi/3, ... % рад, аргумент перицентра
                'u0',pi/2); % начальный аргумент широты
% Углы относительно ОрбСК
psi0 = 0; % рад
theta0 = pi/6; % рад
phi0 = 0; % рад
[r0,v0] = elem2xyz(params.i, params.e, params.p, params.Omega, ...
                    params.omega, params.u0, params.mu);
START_X = [r0;  % радиус-вектор, км
           v0]; % скорость, км/с
START_Q_orb = angle2quat(psi0,theta0,phi0,'ZXZ')';
START_Q = qmult(dcm2quat(IF2OF(params.Omega,params.i,params.u0))',START_Q_orb);
START_W = 1e-3*[1;1;1]; % угловая скорость в ССК, рад/с

T = 1e2; % с
dt = 0.005;