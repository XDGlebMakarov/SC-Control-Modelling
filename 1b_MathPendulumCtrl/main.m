% close all
clear
clc

init

t = 0:dt:T;
N = length(t);
x = zeros(2,N);
x(:,1) = x_START;
E = zeros(1,N);
E(1) = Energy(x(:,1),params);
dE = zeros(1,N);
t_ctrl = 0:dt_ctrl:T;
u = zeros(1,length(t_ctrl));
j = 1;
u(1) = ctrl(x_START,phi_REQ,lyap_k1,lyap_k2,u_NOISE_SIGMA,u_MAX,params);

for i=2:N
    if t(i) >= t_ctrl(j+1)
        j = j+1;
        u(j) = ctrl(x(:,i-1),phi_REQ,lyap_k1,lyap_k2,u_NOISE_SIGMA,u_MAX,params);
    end
    x(:,i-1) = x_NOISE(x(:,i-1),x_NOISE_SIGMA);
    x(:,i) = RKstep(x(:,i-1),u(:,j),dt,params);
    E(i) = Energy(x(:,i),params);
    dE(i) = E(i)/E(i-1) - 1;
end

%%
figure('units','normalized','outerposition',[0 0 1 1]);
subplot(2,3,1)
plot(t,180/pi*x(1,:))
grid on
xlabel('t, s')
ylabel('\phi, ^\circ')
title('\phi(t)')

subplot(2,3,2)
plot(t,180/pi*x(2,:))
grid on
xlabel('t, s')
ylabel('\omega, ^\circ/s')
title('\omega(t)')

subplot(2,3,3)
grid on
hold on
plot([0 T],180/pi*u_MAX*ones(1,2),'r--',LineWidth=1)
plot([0 T],-180/pi*u_MAX*ones(1,2),'r--',LineWidth=1)
plot(t_ctrl,180/pi*u,'k-',LineWidth=1.5)
ylim(1.1*180/pi*u_MAX*[-1 1])
xlabel('t, s')
ylabel('u, ^\circ/s^2')
title('u(t)')

subplot(2,3,4)
plot(t,E)
grid on
xlabel('t, s')
ylabel('E, E_0-units')
title('E(t)')

subplot(2,3,5)
plot(t,dE)
grid on
xlabel('t, s')
ylabel('dE/E')
title('dE/E(t)')

subplot(2,3,6)
plot(180/pi*x(1,:),180/pi*x(2,:))
grid on
xlabel('\phi, ^\circ')
ylabel('\omega, ^\circ/s')
title('Phase plane \phi-\omega')

sgtitle('Math Pendulum (with control)')