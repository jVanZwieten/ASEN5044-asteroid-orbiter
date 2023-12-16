close all
format longg
clc
clear

data = load("orbitdetermination-finalproj_data_2023_11_14.mat");
setGlobalVariables()



dx0 = zeros(6,1); %[1e-5 1e-5 1e-5 1e-7 1e-7 1e-7]'; 
gammaW = [zeros(3); zeros(3)];
% gammaW = CTsys.Gamma();   Gamma is not defined correctly right now
W = sigma_w^2;
Q = DTsys.noiseMat(W,delT_observation);



% :::::Define Nonlinear matrices:::::

[NL_state,noisy_NL_state] = LinearizedKalmanFilter.genNLState(dx0,gammaW);
[y_sim,y_table] = LinearizedKalmanFilter.genNLMeasFromRealMeas(NL_state,data);



% :::::Calculate filtered measurements and NEES test:::::

Qkf=Q;
%Qkf(4:6,4:6) = 1e-30*eye(3); % only changes velocity covariance
[xP,P_P,filt_total_state,~] = LinearizedKalmanFilter.LKF(NL_state,noisy_NL_state,dx0,y_sim,y_table,data,Qkf);



% :::::Plots:::::

time = (0:delT_observation:t_end)/3600;
NLtime = (0:delT_integration:(t_end+60))/3600;

% plot typical noisy traj vs noiseless traj
% figure()
% subplot(611)
% sgtitle('Typical Noisy Truth Trajectory')
% plot(NLtime,NL_state(1,:),NLtime,noisy_NL_state(1,:))
% ylabel('x (km)')
% legend('Noiseless','Noisy')
% 
% subplot(612)
% plot(NLtime,NL_state(2,:),NLtime,noisy_NL_state(2,:))
% ylabel('y (km)')
% 
% subplot(613)
% plot(NLtime,NL_state(3,:),NLtime,noisy_NL_state(3,:))
% ylabel('z (km)')
% 
% subplot(614)
% plot(NLtime,NL_state(4,:),NLtime,noisy_NL_state(4,:))
% ylabel('xdot (km/s)')
% 
% subplot(615)
% plot(NLtime,NL_state(5,:),NLtime,noisy_NL_state(5,:))
% ylabel('ydot (km/s)')
% 
% subplot(616)
% plot(NLtime,NL_state(6,:),NLtime,noisy_NL_state(6,:))
% ylabel('zdot (km/s)')
% xlabel('Time (hours)')

% plot noisy simulated data
% figure()
% for ii = 1:5
%     meastime = y_table(find(y_table(:,2)==ii),1);
%     noisy_u = y_table(find(y_table(:,2)==ii),3);
%     noisy_v = y_table(find(y_table(:,2)==ii),4);
%     plot(meastime/3600,noisy_u,'x')
% end

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



% plot nonlinear and estimated state

figure()
sgtitle('Position')
subplot(311)
plot(time,NL_state(1,1:10:end))
hold on
plot(time,filt_total_state(1,:))
plot(time,2*sqrt(reshape(P_P(1,1,1:end-1),[],1))+filt_total_state(1,:)','black --')
plot(time,-2*sqrt(reshape(P_P(1,1,1:end-1),[],1))+filt_total_state(1,:)','black --')
ylabel('x (km)')
legend('Nominal Orbit','Total Filter State','\pm2\sigma bounds')
%ylim([-1e-4 1e-4])

subplot(312)
plot(time,NL_state(2,1:10:end))
hold on
plot(time,filt_total_state(2,:))
plot(time,2*sqrt(reshape(P_P(2,2,1:end-1),[],1))+filt_total_state(2,:)','black --')
plot(time,-2*sqrt(reshape(P_P(2,2,1:end-1),[],1))+filt_total_state(2,:)','black --')
ylabel('y (km)')
%ylim([-1e-3 1e-3])

subplot(313)
plot(time,NL_state(3,1:10:end))
hold on
plot(time,filt_total_state(3,:))
plot(time,2*sqrt(reshape(P_P(3,3,1:end-1),[],1))+filt_total_state(3,:)','black --')
plot(time,-2*sqrt(reshape(P_P(3,3,1:end-1),[],1))+filt_total_state(3,:)','black --')
ylabel('z (km)')
xlabel('Time (hours)')
%ylim([-1e-3 1e-3])

figure()
sgtitle('Velocity')
subplot(311)
plot(time,NL_state(4,1:10:end))
hold on
plot(time,filt_total_state(4,:))
plot(time,2*sqrt(reshape(P_P(4,4,1:end-1),[],1))+filt_total_state(4,:)','black --')
plot(time,-2*sqrt(reshape(P_P(4,4,1:end-1),[],1))+filt_total_state(4,:)','black --')
ylabel('xdot (km)')
legend('Nominal Orbit','Total Filter State','\pm2\sigma bounds')
%ylim([-1e-4 1e-4])

subplot(312)
plot(time,NL_state(5,1:10:end))
hold on
plot(time,filt_total_state(5,:))
plot(time,2*sqrt(reshape(P_P(5,5,1:end-1),[],1))+filt_total_state(5,:)','black --')
plot(time,-2*sqrt(reshape(P_P(5,5,1:end-1),[],1))+filt_total_state(5,:)','black --')
ylabel('ydot (km)')
%ylim([-1e-3 1e-3])

subplot(313)
plot(time,NL_state(6,1:10:end))
hold on
plot(time,filt_total_state(6,:))
plot(time,2*sqrt(reshape(P_P(6,6,1:end-1),[],1))+filt_total_state(6,:)','black --')
plot(time,-2*sqrt(reshape(P_P(6,6,1:end-1),[],1))+filt_total_state(6,:)','black --')
ylabel('zdot (km)')
xlabel('Time (hours)')
%ylim([-1e-3 1e-3])

RSS_error = sqrt(state_est_error(1,:).^2 + state_est_error(2,:).^2 + state_est_error(3,:).^2);

figure()
plot(time,RSS_error)
xlabel('Time (hours)')
ylabel('\epsilon (km)')
title('RSS Position Error')

