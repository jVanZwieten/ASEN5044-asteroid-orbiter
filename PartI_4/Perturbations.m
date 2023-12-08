clc
clear
close all
format short

addpath(genpath(fileparts(pwd)))

global mu rsa phi0 rho Am
mu = 4.892E-9;
rsa = [1.5E8 0 0]';
phi0 = 1E14;
rho = 0.4;
Am = (1/62)*10^(-6);

delTobs = 600;
delTint = 60; % s
tEnd = 72*60*60; % 72h -> s
sec2hr = 3600;

r0 = [0 -1 0]';
rdot0 = [0 0 sqrt(mu/norm(r0))]';
state0 = [r0; rdot0];

NL_state = zeros(6,length(1:delTint:tEnd)+1);
NL_state(:,1) = state0;

% compute NL states at 60s time step
for i=1:tEnd/delTint
    NL_state(:,i+1) = rk4_state(NL_state(:,i),delTint);
end

dx0 = [1e-5 1e-5 1e-5 1e-7 1e-7 1e-7]';


asrp = f.solarRadPress();
u = asrp(4:6);

dx = zeros(6,length(1:delTobs:tEnd)+1);
linear_state = zeros(6,length(1:delTobs:tEnd)+1);
linear_state(:,1) = state0;
dx(:,1) = dx0;

% linearize about nominal orbit at each 600s time step
for i = 1:tEnd/delTobs
    x = NL_state(1,10*(i-1)+1);
    y = NL_state(2,10*(i-1)+1);
    z = NL_state(3,10*(i-1)+1);
    r = norm(NL_state(1:3,10*(i-1)+1));

    [A,B] = CTsys.dynMat(x,y,z,r);
    [F,G] = DTsys.dynMat(A,B,delTobs);

    linear_state(:,i+1) = F*linear_state(:,i);
    dx(:,i+1) = F*dx(:,i);

end

% plot perturbed states vs time
tVec = (0:delTobs:tEnd)/sec2hr;
figure()
sgtitle('Linearized perturbations')
subplot(6,1,1)
plot(tVec,dx(1,:))
ylabel('\Deltax (km)')

subplot(6,1,2)
plot(tVec,dx(2,:))
ylabel('\Deltay (km)')

subplot(6,1,3)
plot(tVec,dx(3,:))
ylabel('\Deltaz (km/s)')

subplot(6,1,4)
plot(tVec,dx(4,:))
ylabel('\Deltaxdot (km/s)')


subplot(6,1,5)
plot(tVec,dx(5,:))
ylabel('\Deltaydot (km/s)')

subplot(6,1,6)
plot(tVec,dx(6,:))
ylabel('\Deltazdot (km)')
xlabel('Time (hours)')




