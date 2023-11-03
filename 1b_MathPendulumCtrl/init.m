params = struct('m',1 + normrnd(0,1e-4), ... % кг
                'l',1 + normrnd(0,1e-4), ... % м
                'g',10 + normrnd(0,1e-4));   % м/с^2
x_START = [pi/4; % phi
           0];   % omega, 1/с
phi_REQ = pi/6;
x_NOISE_SIGMA = [pi*1e-4;
                 pi*1e-5];
T = 50; % с
dt = 1e-2;
dt_ctrl = 5e-2;
lyap_k1 = 5e-1*params.g/(params.l*sqrt(2));
lyap_k2 = 2;
u_MAX = 2*pi;
u_NOISE_SIGMA = pi*1e-6;