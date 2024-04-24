params = struct('J',diag([0.15,0.15,0.10])); % тензор инерции, кг*м^2
m = 3; % кг
psi0 = 0; % рад
theta0 = pi/6; % рад
phi0 = 0; % рад
START_CONDS = containers.Map();
START_CONDS('Euler') = [psi0;theta0;phi0];
START_CONDS('CosMatrix') = angle2dcm(psi0,theta0,phi0,'ZXZ');
START_CONDS('Quat') = angle2quat(psi0,theta0,phi0,'ZXZ')';
START_W = [0;0;1]; % рад/с в ССК

T = 1e2; % с
dt = 0.005;