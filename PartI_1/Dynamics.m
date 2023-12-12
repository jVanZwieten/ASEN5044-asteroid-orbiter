clc
clear
close all
format longg

setGlobalVariables()
global X_0 delT_integration t_end

X = numerical.propagate(X_0, delT_integration, t_end);

tVec = (1:delT_integration:t_end)/60/60;

figure()
subplot(6,1,1)
plot(tVec, X(1,1:end-1))
title('Nonlinear States vs Time')
xlabel("Time (h)")
ylabel("x (km)")

subplot(6,1,2)
plot(tVec, X(2,1:end-1))
xlabel("Time (h)")
ylabel("y (km)")

subplot(6,1,3)
plot(tVec, X(3,1:end-1))
xlabel("Time (h)")
ylabel("z (km)")

subplot(6,1,4)
plot(tVec, X(4,1:end-1))
xlabel("Time (h)")
ylabel("$\dot{x}$ (km/s)", 'Interpreter', 'latex')

subplot(6,1,5)
plot(tVec, X(5,1:end-1))
xlabel("Time (h)")
ylabel("$\dot{y}$ (km/s)", 'Interpreter', 'latex')

subplot(6,1,6)
plot(tVec, X(6,1:end-1))
xlabel("Time (h)")
ylabel("$\dot{z}$ (km/s)", 'Interpreter', 'latex')