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

FI = zeros(10,N);
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
psi_angs = angs(1,:);
theta_angs = angs(2,:);
phi_angs = angs(3,:);
[psi_CM  ,theta_CM  ,phi_CM  ] = dcm2angle(CM,'ZXZ');
[psi_quat,theta_quat,phi_quat] = quat2angle(quat','ZXZ');
psi_CM = psi_CM';
psi_quat = psi_quat';
theta_CM = theta_CM';
theta_quat = theta_quat';
phi_CM = phi_CM';
phi_quat = phi_quat';

figure('units','normalized','outerposition',[0 0 1 1]);
subplot(2,3,1)
plot(t,psi_angs,'r',LineWidth=1.5)
hold on
plot(t,psi_CM,'g',LineWidth=1.5)
plot(t,psi_quat,'b',LineWidth=1.5)
grid on
legend({'Euler','DCM','Quat'})
title('\psi')

subplot(2,3,2)
plot(t,theta_angs,'r',LineWidth=1.5)
hold on
plot(t,theta_CM,'g',LineWidth=1.5)
plot(t,theta_quat,'b',LineWidth=1.5)
grid on
legend({'Euler','DCM','Quat'})
title('\theta')

subplot(2,3,3)
plot(t,phi_angs,'r',LineWidth=1.5)
hold on
plot(t,phi_CM,'g',LineWidth=1.5)
plot(t,phi_quat,'b',LineWidth=1.5)
grid on
legend({'Euler','DCM','Quat'})
title('\phi')

subplot(2,3,4)
plot(t,psi_angs-psi_CM,'r',LineWidth=1.5)
hold on
plot(t,psi_angs-psi_quat,'g',LineWidth=1.5)
plot(t,psi_CM-psi_quat,'b',LineWidth=1.5)
grid on
legend({'Euler-DCM','Euler-Quat','DCM-Quat'})
title('\Delta\psi')

subplot(2,3,5)
plot(t,theta_angs-theta_CM,'r',LineWidth=1.5)
hold on
plot(t,theta_angs-theta_quat,'g',LineWidth=1.5)
plot(t,theta_CM-theta_quat,'b',LineWidth=1.5)
grid on
legend({'Euler-DCM','Euler-Quat','DCM-Quat'})
title('\Delta\theta')

subplot(2,3,6)
plot(t,phi_angs-phi_CM,'r',LineWidth=1.5)
hold on
plot(t,phi_angs-phi_quat,'g',LineWidth=1.5)
plot(t,phi_CM-phi_quat,'b',LineWidth=1.5)
grid on
legend({'Euler-DCM','Euler-Quat','DCM-Quat'})
title('\Delta\phi')

figure
plot(t,w(1,:),'r',LineWidth=1.5)
hold on
plot(t,w(2,:),'g',LineWidth=1.5)
plot(t,w(3,:),'b',LineWidth=1.5)
grid on
legend({'\omega_x','\omega_y','\omega_z'})
title('\omega')

dh = FI(1,2:end)-FI(1,1:end-1);
dK_angs = FI(2: 4,2:end)-FI(2: 4,1:end-1);
dK_CM   = FI(5: 7,2:end)-FI(5: 7,1:end-1);
dK_quat = FI(8:10,2:end)-FI(8:10,1:end-1);
figure('units','normalized','outerposition',[0 0 1 1]);
subplot(2,2,1)
plot(t(2:end),dh,LineWidth=1.5)
grid on
xlabel('t, s')
ylabel('dE')
title('\Delta Kinetic energy')

subplot(2,2,2)
plot(t(2:end),dK_angs(1,:),'r',LineWidth=1.5)
hold on
plot(t(2:end),dK_CM(1,:),'g',LineWidth=1.5)
plot(t(2:end),dK_quat(1,:),'b',LineWidth=1.5)
grid on
xlabel('t, s')
ylabel('dK_x')
legend({'Euler','DCM','Quat'})
title('\Delta Kinetic moment K_x')

subplot(2,2,3)
plot(t(2:end),dK_angs(2,:),'r',LineWidth=1.5)
hold on
plot(t(2:end),dK_CM(2,:),'g',LineWidth=1.5)
plot(t(2:end),dK_quat(2,:),'b',LineWidth=1.5)
grid on
xlabel('t, s')
ylabel('dK_y')
legend({'Euler','DCM','Quat'})
title('\Delta Kinetic moment K_y')

subplot(2,2,4)
plot(t(2:end),dK_angs(3,:),'r',LineWidth=1.5)
hold on
plot(t(2:end),dK_CM(3,:),'g',LineWidth=1.5)
plot(t(2:end),dK_quat(3,:),'b',LineWidth=1.5)
grid on
xlabel('t, s')
ylabel('dK_z')
legend({'Euler','DCM','Quat'})
title('\Delta Kinetic moment K_z')

function Fi = firstints(angs,CM,q,w,params)
    h = w'*params.J*w/2;
    K_bf = params.J*w;
    K_angs = euler2CM(angs)'*K_bf;
    K_CM = CM'*K_bf;
    K_quat = quattrans(qconj(q),K_bf);
    Fi = [h;K_angs;K_CM;K_quat];
end