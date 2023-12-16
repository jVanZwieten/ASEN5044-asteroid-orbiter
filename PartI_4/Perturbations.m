clc
clear
close all
format short

global mu rsa phi0 rho Am omega umin umax
mu = 4.892E-9;
rsa = [1.5E8 0 0]';
phi0 = 1E14;
rho = 0.4;
Am = (1/62)*10^(-6);
omega = 2*pi/(4.296057*3600);
umin = 0;
umax = 1024;

delTobs = 600;
tEnd = 72*60*60;

delTint = 60; % s
tEnd = 72*60*60; % 72h -> s

r0 = [0 -1 0]';
rdot0 = [0 0 sqrt(mu/norm(r0))]';
state0 = [r0; rdot0];

NL_state = zeros(6,length(1:delTint:tEnd)+1);
NL_state(:,1) = state0;

% compute NL states at 60s time step
for i=1:tEnd/delTint
    NL_state(:,i+1) = numerical.rk4_state(NL_state(:,i),delTint);
end

dx0 = [1e-5 1e-5 1e-5 1e-7 1e-7 1e-7]';


asrp = f.solarRadPress();
u = asrp(4:6);

dx = zeros(6,length(1:delTobs:tEnd)+1);
linear_state = zeros(6,length(1:delTobs:tEnd)+1);
linear_state(:,1) = state0;
dx(:,1) = dx0;

% linearize about nominal orbit at each 600s time step
for i = 1:tEnd/600
    x = NL_state(1,10*(i-1)+1);
    y = NL_state(2,10*(i-1)+1);
    z = NL_state(3,10*(i-1)+1);
    r = norm(NL_state(1:3,10*(i-1)+1));

    A = [zeros(3)  eye(3);
        zeros(3) zeros(3)];
    A(4:6,1:3) = [3*mu*x^2/r^5 - mu/r^3    3*mu*x*y/r^5          3*mu*x*z;
                 3*mu*y*x/r^5       3*mu*y^2/r^5 - mu/r^3   3*mu*y*z/r^5;
                 3*mu*z*x/r^5             3*mu*y*z/r^5    3*mu*z^2/r^5 - mu/r^3];

    B = [zeros(3,3); eye(3,3)];

    Ahat = [A B;
        zeros(length(B(1,:)),length(A(1,:))),zeros(length(B(1,:)),length(B(1,:)))];
    expmA = expm(Ahat*delTobs);

    F = expmA(1:length(A(:,1)),1:length(A(1,:)));
    G = expmA(1:length(A(:,1)),length(A(1,:))+1:length(A(1,:))+length(B(1,:)));

    linear_state(:,i+1) = F*linear_state(:,i);
    dx(:,i+1) = F*dx(:,i);

end



tVec = (0:delTobs:tEnd)/60/60;
%%
% plot nonlinear and perturbed states vs time
figure()
subplot(6,1,1)
plot(tVec, NL_state(1,1:10:end))
hold on
plot(tVec, NL_state(1,1:10:end)+dx(1,:))
title('Linear and Nonlinear Propagation of States','FontSize',16)
ylabel("x (km)",'FontSize',13)
legend('Nonlinear State','Linear State',location = 'southeast')

subplot(6,1,2)
plot(tVec, NL_state(2,1:10:end))
hold on
plot(tVec, NL_state(2,1:10:end)+dx(2,:))
ylabel("y (km)",'FontSize',13)

subplot(6,1,3)
plot(tVec, NL_state(3,1:10:end))
hold on
plot(tVec, NL_state(3,1:10:end)+dx(3,:))
ylabel("z (km)",'FontSize',13)

subplot(6,1,4)
plot(tVec, NL_state(4,1:10:end))
hold on
plot(tVec, NL_state(4,1:10:end)+dx(4,:))
ylabel("$\dot{x}$ ($\frac{km}{s}$)", 'Interpreter','latex','FontSize',14)

subplot(6,1,5)
plot(tVec, NL_state(5,1:10:end))
hold on
plot(tVec, NL_state(5,1:10:end)+dx(5,:))
ylabel("$\dot{y}$ ($\frac{km}{s}$)", 'Interpreter','latex','FontSize',14)

subplot(6,1,6)
plot(tVec, NL_state(6,1:10:end))
hold on
plot(tVec, NL_state(6,1:10:end)+dx(6,:))
xlabel("Time (h)",'FontSize',13)
ylabel("$\dot{z}$ ($\frac{km}{s}$)", 'Interpreter','latex','FontSize',14)



% plot perturbed states vs time
figure()
sgtitle('Linearized perturbations','FontSize',16)
subplot(6,1,1)
plot(tVec,dx(1,:))
ylabel('\Deltax (km)','FontSize',13)

subplot(6,1,2)
plot(tVec,dx(2,:))
ylabel('\Deltay (km)','FontSize',13)

subplot(6,1,3)
plot(tVec,dx(3,:))
ylabel('\Deltaz (km/s)','FontSize',13)

subplot(6,1,4)
plot(tVec,dx(4,:))
ylabel("$\dot{x}$ ($\frac{km}{s}$)", 'Interpreter','latex','FontSize',14)


subplot(6,1,5)
plot(tVec,dx(5,:))
ylabel("$\dot{y}$ ($\frac{km}{s}$)", 'Interpreter','latex','FontSize',14)

subplot(6,1,6)
plot(tVec,dx(6,:))
ylabel("$\dot{z}$ ($\frac{km}{s}$)", 'Interpreter','latex','FontSize',14)
xlabel('Time (hours)','FontSize',13)
