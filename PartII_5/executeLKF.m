close all
format longg
clc
clear

%addpath(genpath(fileparts(pwd)))

data = load("orbitdetermination-finalproj_data_2023_11_14.mat");
setGlobalVariables()



dx0 = [1e-7 1e-7 1e-7 1e-9 1e-9 1e-9]';    %[1e-7 1e-7 1e-7 1e-9 1e-9 1e-9]'; 
asrp = f.solarRadPress();
asrp = asrp(4:6);

gammaW = [zeros(3); eye(3)];
% gammaW = CTsys.Gamma();   Gamma is not defined correctly right now
W = sigma_w^2;
Q = DTsys.noiseMat(W,delT_observation);



% :::::Define Nonlinear matrices:::::

[NL_state,noisy_NL_state] = LinearizedKalmanFilter.genNLState(dx0,gammaW);
[y_noiseless,y_table,y_actual_noisy] = LinearizedKalmanFilter.genNLMeas(NL_state,noisy_NL_state,data);
%[y_actual_noiseless,y_actual_noisy]=LinearizedKalmanFilter.genNLMeas(noisy_NL_state,data);



% :::::Calculate filtered measurements and NEES test:::::

Qkf=Q;  %cov(noisy_NL_state'-NL_state');
Qkf(1,1) = 5e3*Qkf(1,1);
Qkf(2:3,2:3) = 1e7*Qkf(2:3,2:3);
Qkf(4,4) = 1e3*Qkf(4,4); % only changes velocity covariance
Qkf(5:6,5:6) = 5e4*Qkf(5:6,5:6);
Nsimruns = 10;
nMeas = length(1:delT_observation:t_end)+1;
NEES_samps = zeros(Nsimruns,nMeas);
NIS_samps = zeros(Nsimruns,nMeas);

for k = 1:Nsimruns 
    [xP,P_P,filt_total_state,NEES_hist,NIS_hist] = LinearizedKalmanFilter.LKF(NL_state,noisy_NL_state,dx0,y_table,y_actual_noisy,data,Qkf);
    NEES_samps(k,:)=NEES_hist;
    NIS_samps(k,:)=NIS_hist;
end



% :::::Plots:::::

% Psig = P_P(1,1,:);
time = (0:delT_observation:t_end)/3600;
NLtime = (0:delT_integration:(t_end+60))/3600;
% 
% figure()
% plot(time,filt_total_state(1,:),'red')
% hold on
% plot(time,NL_state(1,1:10:end-1),'bx')
% 
% plot(time,2*sqrt(Psig(1,1:end-1))+filt_total_state(1,:),'black --')
% plot(time,-2*sqrt(Psig(1,1:end-1))+filt_total_state(1,:),'black --')


% plot typical noisy traj vs noiseless traj
figure()
subplot(611)
sgtitle('Typical Noisy Truth Trajectory')
plot(NLtime,NL_state(1,:),NLtime,noisy_NL_state(1,:),'--')
ylabel('x (km)')
legend('Noiseless','Noisy')

subplot(612)
plot(NLtime,NL_state(2,:),NLtime,noisy_NL_state(2,:),'--')
ylabel('y (km)')

subplot(613)
plot(NLtime,NL_state(3,:),NLtime,noisy_NL_state(3,:),'--')
ylabel('z (km)')

subplot(614)
plot(NLtime,NL_state(4,:),NLtime,noisy_NL_state(4,:),'--')
ylabel('xdot (km/s)')

subplot(615)
plot(NLtime,NL_state(5,:),NLtime,noisy_NL_state(5,:),'--')
ylabel('ydot (km/s)')

subplot(616)
plot(NLtime,NL_state(6,:),NLtime,noisy_NL_state(6,:),'--')
ylabel('zdot (km/s)')
xlabel('Time (hours)')

% plot noisy simulated data
figure()
for ii = 1:3
    meastime = y_actual_noisy(find(y_actual_noisy(:,2)==ii),1);
    noisy_u = y_actual_noisy(find(y_actual_noisy(:,2)==ii),3);
    noisy_v = y_actual_noisy(find(y_actual_noisy(:,2)==ii),4);
    noiseless_u = y_noiseless(find(y_noiseless(:,2)==ii),3);
    noiseless_v = y_noiseless(find(y_noiseless(:,2)==ii),4);
    subplot(211)
    plot(meastime/3600,noisy_u,'x')
    hold on
    plot(meastime/3600,noiseless_u,'o')
    ylabel('u (pixels)')

    subplot(212)    
    plot(meastime/3600,noisy_v,'x')
    hold on
    plot(meastime/3600,noiseless_v,'o')
    ylabel('v (pixels)')
end
sgtitle('Simulated NL Measurements (lmks 1-3)')
legend('Noisy','Noiseless')

xlabel('Time(hours)')


% plot state estimation errors
state_est_error = filt_total_state - NL_state(:,1:10:end);
figure()
sgtitle('Position state estimation error')
subplot(311)
plot(time,state_est_error(1,:))
hold on
plot(time,2*sqrt(reshape(P_P(1,1,1:end-1),[],1)),'black --')
plot(time,-2*sqrt(reshape(P_P(1,1,1:end-1),[],1)),'black --')
ylabel('x (km)')
legend('Error','\pm2\sigma bounds')
%xlim([0 10])
%ylim([-1e-4 1e-4])

subplot(312)
plot(time,state_est_error(2,:))
hold on
plot(time,2*sqrt(reshape(P_P(2,2,1:end-1),[],1)),'black --')
plot(time,-2*sqrt(reshape(P_P(2,2,1:end-1),[],1)),'black --')
ylabel('y (km)')
%xlim([0 10])
%ylim([-1e-3 1e-3])

subplot(313)
plot(time,state_est_error(3,:))
hold on
plot(time,2*sqrt(reshape(P_P(3,3,1:end-1),[],1)),'black --')
plot(time,-2*sqrt(reshape(P_P(3,3,1:end-1),[],1)),'black --')
ylabel('z (km)')
xlabel('Time (hours)')
%xlim([0 10])
%ylim([-1e-3 1e-3])

figure()
sgtitle('Velocity state estimation error')
subplot(311)
plot(time,state_est_error(4,:))
hold on
plot(time,2*sqrt(reshape(P_P(4,4,1:end-1),[],1)),'black --')
plot(time,-2*sqrt(reshape(P_P(4,4,1:end-1),[],1)),'black --')
ylabel('xdot (km)')
legend('Error','\pm2\sigma bounds')
%xlim([0 10])
%ylim([-1e-4 1e-4])

subplot(312)
plot(time,state_est_error(5,:))
hold on
plot(time,2*sqrt(reshape(P_P(5,5,1:end-1),[],1)),'black --')
plot(time,-2*sqrt(reshape(P_P(5,5,1:end-1),[],1)),'black --')
ylabel('ydot (km)')
%xlim([0 10])
%ylim([-1e-3 1e-3])

subplot(313)
plot(time,state_est_error(6,:))
hold on
plot(time,2*sqrt(reshape(P_P(6,6,1:end-1),[],1)),'black --')
plot(time,-2*sqrt(reshape(P_P(6,6,1:end-1),[],1)),'black --')
ylabel('zdot (km)')
xlabel('Time (hours)')
%xlim([0 10])
%ylim([-1e-3 1e-3])

% plot(NLtime,NL_state(1,:),time,filt_total_state(1,:))
% rows = 1:length(y_table);
% 
% figure()
% plot(y_table(rows,3),y_table(rows,4),'blue x')
% hold on
% plot(linear_meas(rows,3),linear_meas(rows,4),'red o')
% xlim([umin(1) umax(1)])
% ylim([umin(1) umax(1)])

% plotter = [];
% for i=1:length(NL_state)
%     if(mod())
% end

% plot(1:length(linear_state),linear_state(1,:))
% hold on
% plot(1:length(NL_state),NL_state(1,:))

%%DO NEES TEST:
epsNEESbar = mean(NEES_samps,1);
alphaNEES = 0.05; %%significance level
Nnx = Nsimruns*6;
%%compute intervals:
r1x = chi2inv(alphaNEES/2, Nnx )./ Nsimruns
r2x = chi2inv(1-alphaNEES/2, Nnx )./ Nsimruns

figure()
plot(epsNEESbar,'ro','MarkerSize',6,'LineWidth',2),hold on
plot(r1x*ones(size(epsNEESbar)),'r--','LineWidth',2)
plot(r2x*ones(size(epsNEESbar)),'r--','LineWidth',2)
ylabel('NEES statistic, $\bar{\epsilon}_x$','Interpreter','latex', 'FontSize',14)
xlabel('time step, k','FontSize',14)
title('NEES Estimation Results','FontSize',14)
legend('NEES @ time k', 'r_1 bound', 'r_2 bound'),grid on
ylim([r1x-2 r2x+2])


%%DO NIS TEST:
epsNISbar = mean(NIS_samps,1);
alphaNIS = 0.05; 
Nny = Nsimruns*2;
%%compute intervals:
r1y = chi2inv(alphaNIS/2, Nny )./ Nsimruns
r2y = chi2inv(1-alphaNIS/2, Nny )./ Nsimruns

figure()
plot(epsNISbar,'bo','MarkerSize',6,'LineWidth',2),hold on
plot(r1y*ones(size(epsNISbar)),'b--','LineWidth',2)
plot(r2y*ones(size(epsNISbar)),'b--','LineWidth',2)
ylabel('NIS statistic, $\bar{\epsilon}_y$','Interpreter','latex','FontSize',14)
xlabel('time step, k','FontSize',14)
title('NIS Estimation Results','FontSize',14)
legend('NIS @ time k', 'r_1 bound', 'r_2 bound'),grid on
%ylim([r1y-2 r2y+2])