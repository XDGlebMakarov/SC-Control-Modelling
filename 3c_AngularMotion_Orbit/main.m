close all
clear
clc

init

t = 0:dt:T;
N = length(t);
x = zeros(6,N);
x(:,1) = START_X;
Q_orb = zeros(4,N);
Q_orb(:,1) = START_Q_orb;
Q = zeros(4,N);
Q(:,1) = START_Q;
w = zeros(3,N);
w(:,1) = START_W;

FI = zeros(10,N);
FI(:,1) = firstints(x(:,1),Q(:,1),w(:,1),params);

tic()
for i=1:N-1
    [x(:,i+1),Q(:,i+1),w(:,i+1)] = rk4step(x(:,i),Q(:,i),w(:,i),dt,params);
    [inc,~,~,Omega,~,u] = xyz2elem(x(1:3,i+1),x(4:6,i+1),params.mu);
    Q_orb(:,i+1) = qmult(qconj(dcm2quat(IF2OF(Omega,inc,u))'),Q(:,i+1));
    FI(:,i+1) = firstints(x(:,i+1),Q(:,i+1),w(:,i+1),params);
end
toc()

%%
[psi,theta,phi] = quat2angle(Q_orb','ZXZ');
psi = 180/pi*psi;
theta = 180/pi*theta;
phi = 180/pi*phi;

figure('units','normalized','outerposition',[0 0 1 1]);
subplot(1,3,1)
plot(t,psi,LineWidth=1.5)
grid on
xlabel('t, с')
title('\psi, ^\circ',FontSize=20)

subplot(1,3,2)
plot(t,theta,LineWidth=1.5)
grid on
xlabel('t, с')
title('\theta, ^\circ',FontSize=20)

subplot(1,3,3)
plot(t,phi,LineWidth=1.5)
grid on
xlabel('t, с')
title('\phi, ^\circ',FontSize=20)
sgtitle('Euler angles OF2BF')

figure
plot(t,w(1,:),'r',LineWidth=1.5)
hold on
plot(t,w(2,:),'g',LineWidth=1.5)
plot(t,w(3,:),'b',LineWidth=1.5)
grid on
legend({'\omega_x','\omega_y','\omega_z'})
title('\omega')

h = FI(1,:);
c = FI(2:4,:);
f = FI(5:7,:);
K = FI(8:10,:);
figure('units','normalized','outerposition',[0 0 1 1]);
subplot(2,5,1)
plot(t,h,LineWidth=1.5)
grid on
xlabel('t, s')
title('h',FontSize=20)

subplot(2,5,2)
plot(t,c(1,:),'r',LineWidth=1.5)
grid on
xlabel('t, s')
title('c_x',FontSize=20)

subplot(2,5,3)
plot(t,c(2,:),'g',LineWidth=1.5)
grid on
xlabel('t, s')
title('c_y',FontSize=20)

subplot(2,5,4)
plot(t,c(3,:),'b',LineWidth=1.5)
grid on
xlabel('t, s')
title('c_z',FontSize=20)

subplot(2,5,5)
plot(t,f(1,:),'r',LineWidth=1.5)
grid on
xlabel('t, s')
title('f_x',FontSize=20)

subplot(2,5,6)
plot(t,f(2,:),'g',LineWidth=1.5)
grid on
xlabel('t, s')
title('f_y',FontSize=20)

subplot(2,5,7)
plot(t,f(3,:),'b',LineWidth=1.5)
grid on
xlabel('t, s')
title('f_z',FontSize=20)

subplot(2,5,8)
plot(t,K(1,:),'r',LineWidth=1.5)
grid on
xlabel('t, s')
title('K_x',FontSize=20)

subplot(2,5,9)
plot(t,K(2,:),'g',LineWidth=1.5)
grid on
xlabel('t, s')
title('K_y',FontSize=20)

subplot(2,5,10)
plot(t,K(3,:),'b',LineWidth=1.5)
grid on
xlabel('t, s')
title('K_z',FontSize=20)

sgtitle('First Integrals')

function Fi = firstints(x,Q,w,params)
    J = params.J;
    mu = params.mu;
    DCM = quat2dcm(Q');
    r = x(1:3);
    nr = norm(r);
    v = x(4:6);
    c = cross(r,v);
    h = w'*J*w/2 + 3/2*mu/nr^5*(DCM*r)'*J*(DCM*r) - w'*J*DCM*c/nr^2;
    f = cross(v,c) - mu*r/nr;
    K = J*w;
    Fi = [h;c;f;K];
end