close all
clear
clc

init

t = 0:dt:T;
N = length(t);
angs = zeros(3,N);
angs(:,1) = START_CONDS('Euler');
CM = zeros(3,3,N);
CM(:,:,1) = START_CONDS('CosMatrix');
quat = zeros(4,N);
quat(:,1) = START_CONDS('Quat');
w = zeros(3,N);
w(:,1) = START_W;

FI = zeros(7,N);
FI(:,1) = firstints(angs(:,1),CM(:,:,1),quat(:,1),w(:,1),params);

tic()
for i=1:N-1
    [angs(:,i+1),w(:,i+1)] = rk4step(angs(:,i),w(:,i),dt,params);
    [CM(:,:,i+1),~       ] = rk4step(CM(:,:,i),w(:,i),dt,params);
    [quat(:,i+1),~       ] = rk4step(quat(:,i),w(:,i),dt,params);
    FI(:,i+1) = firstints(angs(:,i+1),CM(:,:,i+1),quat(:,i+1),w(:,i+1),params);
end
toc()

%%
psi_angs = 180/pi*angs(1,:);
theta_angs = 180/pi*angs(2,:);
phi_angs = 180/pi*angs(3,:);
[psi_CM  ,theta_CM  ,phi_CM  ] = dcm2angle(CM,'ZXZ');
[psi_quat,theta_quat,phi_quat] = quat2angle(quat','ZXZ');
psi_CM = 180/pi*psi_CM';
psi_quat = 180/pi*psi_quat';
theta_CM = 180/pi*theta_CM';
theta_quat = 180/pi*theta_quat';
phi_CM = 180/pi*phi_CM';
phi_quat = 180/pi*phi_quat';

figure('units','normalized','outerposition',[0 0 1 1]);
subplot(2,3,1)
plot(t,psi_angs,'r',LineWidth=1.5)
hold on
plot(t,psi_CM,'g',LineWidth=1.5)
plot(t,psi_quat,'b',LineWidth=1.5)
grid on
legend({'Euler','DCM','Quat'})
title('\psi, ^\circ')

subplot(2,3,2)
plot(t,theta_angs,'r',LineWidth=1.5)
hold on
plot(t,theta_CM,'g',LineWidth=1.5)
plot(t,theta_quat,'b',LineWidth=1.5)
grid on
legend({'Euler','DCM','Quat'})
title('\theta, ^\circ')

subplot(2,3,3)
plot(t,phi_angs,'r',LineWidth=1.5)
hold on
plot(t,phi_CM,'g',LineWidth=1.5)
plot(t,phi_quat,'b',LineWidth=1.5)
grid on
legend({'Euler','DCM','Quat'})
title('\phi, ^\circ')

subplot(2,3,4)
plot(t,psi_angs-psi_CM,'r',LineWidth=1.5)
hold on
plot(t,psi_angs-psi_quat,'g',LineWidth=1.5)
plot(t,psi_CM-psi_quat,'b',LineWidth=1.5)
grid on
legend({'Euler-DCM','Euler-Quat','DCM-Quat'})
title('\Delta\psi, ^\circ')

subplot(2,3,5)
plot(t,theta_angs-theta_CM,'r',LineWidth=1.5)
hold on
plot(t,theta_angs-theta_quat,'g',LineWidth=1.5)
plot(t,theta_CM-theta_quat,'b',LineWidth=1.5)
grid on
legend({'Euler-DCM','Euler-Quat','DCM-Quat'})
title('\Delta\theta, ^\circ')

subplot(2,3,6)
plot(t,phi_angs-phi_CM,'r',LineWidth=1.5)
hold on
plot(t,phi_angs-phi_quat,'g',LineWidth=1.5)
plot(t,phi_CM-phi_quat,'b',LineWidth=1.5)
grid on
legend({'Euler-DCM','Euler-Quat','DCM-Quat'})
title('\Delta\phi, ^\circ')

figure
plot(t,w(1,:),'r',LineWidth=1.5)
hold on
plot(t,w(2,:),'g',LineWidth=1.5)
plot(t,w(3,:),'b',LineWidth=1.5)
grid on
legend({'\omega_x','\omega_y','\omega_z'})
title('\omega')

figure('units','normalized','outerposition',[0 0 1 1]);
subplot(1,3,1)
plot(t,FI(1,:),'r',LineWidth=1.5)
hold on
plot(t,FI(2,:),'g',LineWidth=1.5)
plot(t,FI(3,:),'b',LineWidth=1.5)
grid on
xlabel('t, s')
legend({'Euler','DCM','Quat'})
title('h',FontSize=20)

subplot(1,3,2)
plot(t,FI(4,:),LineWidth=1.5)
grid on
xlabel('t, s')
title('K^{BF}_z',FontSize=20)

subplot(1,3,3)
plot(t,FI(5,:),'r',LineWidth=1.5)
hold on
plot(t,FI(6,:),'g',LineWidth=1.5)
plot(t,FI(7,:),'b',LineWidth=1.5)
grid on
xlabel('t, s')
legend({'Euler','DCM','Quat'})
title('K^{IF}_z',FontSize=20)

sgtitle('First Integrals')

function Fi = firstints(angs,CM,q,w,params)
    CM_quat = quat2dcm(q');
    h_angs = w'*params.J*w/2 + params.m*params.g*params.L*cos(angs(2));
    h_CM   = w'*params.J*w/2 + params.m*params.g*params.L*CM(3,3);
    h_quat = w'*params.J*w/2 + params.m*params.g*params.L*CM_quat(3,3);
    K_bf = params.J*w;
    K_if_angs = angle2dcm(angs(1),angs(2),angs(3),'ZXZ')'*K_bf;
    K_if_CM = CM'*K_bf;
    K_if_quat = CM_quat'*K_bf;
    Fi = [h_angs;h_CM;h_quat;K_bf(3);K_if_angs(3);K_if_CM(3);K_if_quat(3)];
end