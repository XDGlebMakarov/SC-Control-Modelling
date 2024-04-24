params = struct('J',diag([0.15,0.15,0.2]), ... % тензор инерции, кг*м^2
                'm',1, ... % масса, кг
                'g',9.81, ... % м/с^2, ускорение свободного падения
                'L',0.1); % м, расстояние до центра масс
psi0 = 0; % рад
theta0 = pi/6; % рад
phi0 = 0; % рад
START_CONDS = containers.Map();
START_CONDS('Euler') = [psi0;theta0;phi0];
START_CONDS('CosMatrix') = angle2dcm(psi0,theta0,phi0,'ZXZ');
START_CONDS('Quat') = angle2quat(psi0,theta0,phi0,'ZXZ')';
START_W = [0;0;1]; % рад/с в ССК

T = 10; % с
dt = 0.005;