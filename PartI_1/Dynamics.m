clc
clear
close all
format longg

setGlobalVariables()
global X_0 delT_integration t_end

X = numerical.propagate(X_0, delT_integration, t_end);

tVec = (0:delT_integration:t_end)/60/60;

figure()
subplot(6,1,1)
plot(tVec, X(1,1:end-1))
title('Nonlinear States vs Time','FontSize',18)
ylabel("x (km)",'FontSize',15)

subplot(6,1,2)
plot(tVec, X(2,1:end-1))
ylabel("y (km)",'FontSize',15)

subplot(6,1,3)
plot(tVec, X(3,1:end-1))
ylabel("z (km)",'FontSize',15)

subplot(6,1,4)
plot(tVec, X(4,1:end-1))
ylabel("$\dot{x}$ ($\frac{km}{s}$)", 'Interpreter','latex','fontweight','bold','FontSize',16)

subplot(6,1,5)
plot(tVec, X(5,1:end-1))
ylabel("$\dot{y}$ ($\frac{km}{s}$)", 'Interpreter','latex','fontweight','bold','FontSize',16)

subplot(6,1,6)
plot(tVec, X(6,1:end-1))
xlabel("Time (h)",'FontSize',15)
ylabel("$\dot{z}$ ($\frac{km}{s}$)", 'Interpreter','latex','fontweight','bold','FontSize',16)