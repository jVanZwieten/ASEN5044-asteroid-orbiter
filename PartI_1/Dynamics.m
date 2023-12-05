clc
clear
close all
format longg

global muA rsa phi0 rho Am
muA = 4.892E-9;
rsa = [1.5E8 0 0]';
phi0 = 1E14;
rho = 0.4;
Am = (1/62)*10^(-6);

delTint = 60;
tEnd = 72*60*60;

r0 = [0 -1 0]';
rdot0 = [0 0 sqrt(muA/norm(r0))]';
state0 = [r0; rdot0];

state = zeros(6,length(1:delTint:tEnd)+1);
state(:,1) = state0;

for i=1:tEnd/delTint
    state(:,i+1) = rk4_state(state(:,i),delTint);
end

tVec = (1:delTint:tEnd)/60/60;

figure()
subplot(6,1,1)
plot(tVec,state(1,1:end-1))

subplot(6,1,2)
plot(tVec,state(2,1:end-1))

subplot(6,1,3)
plot(tVec,state(3,1:end-1))

subplot(6,1,4)
plot(tVec,state(4,1:end-1))

subplot(6,1,5)
plot(tVec,state(5,1:end-1))

subplot(6,1,6)
plot(tVec,state(6,1:end-1))










