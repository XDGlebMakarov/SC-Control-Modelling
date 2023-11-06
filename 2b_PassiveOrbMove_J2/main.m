close all
clear
clc

init

t = 0:dt:T;
N = length(t);
x = zeros(6,N);
x(:,1) = x_START;
FI = zeros(12,N);
FI(:,1) = firstints(x_START(1:3),x_START(4:6),params.mu,params.delta);
%{
1 - h
2,3,4 - c(x,y,z)
5,6,7 - f(x,y,z)
8 - i
9 - e
10 - p
11 - Omega
12 - omega
%}

tic()
for i=2:N
    x(:,i) = RKstep(x(:,i-1),dt,params);
    FI(:,i) = firstints(x(1:3,i),x(4:6,i),params.mu,params.delta);
end
toc()

%%

figure('units','normalized','outerposition',[0 0 1 1]);
tts = ["h","c_x","c_y","c_z","f_x","f_y","f_z","i","e","p","\Omega","\omega"];
for i=1:8
    subplot(2,4,i)
    plot(t,FI(i,:),'k-',LineWidth=1)
    grid on
    xlabel('t, s')
    title(tts(i))
end

sgtitle('Passive orbital move - 2 body problem, FI')

figure('units','normalized','outerposition',[0 0 1 1]);
for i=1:2
    subplot(2,4,i)
    plot(t,FI(8+i,:),'k-',LineWidth=1)
    grid on
    xlabel('t, s')
    title(tts(8+i))
end
for i=5:6
    subplot(2,4,i)
    plot(t,FI(6+i,:),'k-',LineWidth=1)
    grid on
    xlabel('t, s')
    title(tts(6+i))
end

subplot(2,4,[3 4 7 8])
plot3(x(1,:),x(2,:),x(3,:))
grid on
xlabel('x, km')
ylabel('y, km')
zlabel('z, km')
title('Trajectory')

sgtitle('Passive orbital move - 2 body problem, FI & Trajectory')

function Fi = firstints(r,v,mu,d)
    nr = norm(r);
    c = cross(r,v);
    sinlat = pi/2 - norm(cross([0;0;1],r))/nr;
    h = v'*v/2 - mu/nr + d/nr^3*(1/3-sinlat^2);
    f = cross(v,c) - mu*r/nr;
    [p1,p2,p3,p4,p5,~] = xyz2elem(r,v,mu);
    Fi = [h;c;f;p1;p2;p3;p4;p5];
end