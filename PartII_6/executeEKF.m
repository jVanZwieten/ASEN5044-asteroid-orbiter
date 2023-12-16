clc
clear
close all
format longg
setGlobalVariables()

rng(117)

dx0 = [1e-5 1e-5 1e-5 1e-7 1e-7 1e-7]';
gammaW = [zeros(3); eye(3)];

% no noise test
data = load("orbitdetermination-finalproj_data_2023_11_14.mat");
% X_truth = numerical.propagate(X_0, delT_integration, t_end + delT_observation); % +delT_observation accounts for  t = 0
% X_perturbed = numerical.propagate(X_0+dx0, delT_integration, t_end + delT_observation);
% X_noisy = X_perturbed + noiseMaker(zeros(6,1),kron([0 0; 0 1],(sigma_w^2)*eye(3)),length(X_perturbed))';
[X_truth,X_noisy] = LinearizedKalmanFilter.genNLState(dx0,gammaW);
Q = DTsys.noiseMat(sigma_w^2, delT_observation);
% utilities.plot3D(X_truth(1:3, :), [-1 1; -1 1; -1 1], "no process noise orbit")

integrationObservationRatio = delT_observation/delT_integration;
landmarkPositions = propagateLandmarksInInertialFrame(data.pos_lmks_A, delT_observation, t_end, w_A);
X_truthAtObservationEpochs = X_truth(:, (1 + integrationObservationRatio):integrationObservationRatio:end);
X_noisyAtObservationEpochs = X_noisy(:, (1 + integrationObservationRatio):integrationObservationRatio:end);
R = sigma_u^2*eye(p);
[Y_simulated,~] = measurement.SimulateYData(X_truthAtObservationEpochs, landmarkPositions, data.R_CtoN, R, delT_observation);
[Y_perturbed,Y_noisy] = measurement.SimulateYData(X_noisyAtObservationEpochs, landmarkPositions, data.R_CtoN, R, delT_observation);

W = 1e-18;
% Qkf = DTsys.noiseMat(W, delT_observation);
Qkf = Q;
Qkf(1:3,1:3) = 10*Q(1:3,1:3);
Qkf(4,4) = 8e6*Q(4,4); % only changes velocity covariance
Qkf(5:6,5:6) = 1e9*Q(5:6,5:6);
Nsimruns = 10;
nMeas = length(1:delT_observation:t_end)+1;
NEES_samps = zeros(Nsimruns,nMeas);
NIS_samps = zeros(Nsimruns,nMeas);

for k = 1:Nsimruns 
    [X_estimate, P_estimate,NEES_hist,NIS_hist] = ExtendedKalmanFilter.Filter(X_0, P_0, Qkf, Y_noisy, R, delT_observation, landmarkPositions, data.R_CtoN);
    NEES_samps(k,:)=NEES_hist;
    NIS_samps(k,:)=NIS_hist;
end



T = (-600:delT_observation:t_end)/60/60;

% INPUTS FOR THIS FUNCTION DON'T MATCH REQUIRED INPUTS -->
% utilities.confidencePlots(T, X_estimate, P_estimate, 2, "State Estimate, no noise", XLabels, XUnits)

utilities.errorPlots(T(:, 2:end-1), X_truthAtObservationEpochs, X_estimate(:, 2:end), P_estimate(:, :, 2:end), 2, "Error, no noise", XLabels, XUnits)

epsNEESbar = mean(NEES_samps,1);
alphaNEES = 0.05; %%significance level
Nnx = k*n;
r1x = chi2inv(alphaNEES/2, Nnx )./ Nsimruns;
r2x = chi2inv(1-alphaNEES/2, Nnx )./ Nsimruns;

figure()
plot(epsNEESbar,'ro','MarkerSize',6,'LineWidth',2),hold on
plot(r1x*ones(size(epsNEESbar)),'r--','LineWidth',2)
plot(r2x*ones(size(epsNEESbar)),'r--','LineWidth',2)
ylabel('NEES statistic, $\bar{\epsilon}_x$','Interpreter','latex', 'FontSize',14)
xlabel('time step, k','FontSize',14)
title('NEES Estimation Results','FontSize',14)
legend('NEES @ time k', 'r_1 bound', 'r_2 bound'),grid on
ylim([r1x-2 r2x+2])

% DO NIS TEST:
epsNISbar = mean(NIS_samps,1);
alphaNIS = 0.05; 
Nny = Nsimruns*3;
%%compute intervals:
r1y = chi2inv(alphaNIS/2, Nny )./ Nsimruns;
r2y = chi2inv(1-alphaNIS/2, Nny )./ Nsimruns;

figure()
plot(epsNISbar,'bo','MarkerSize',6,'LineWidth',2),hold on
plot(r1y*ones(size(epsNISbar)),'b--','LineWidth',2)
plot(r2y*ones(size(epsNISbar)),'b--','LineWidth',2)
ylabel('NIS statistic, $\bar{\epsilon}_y$','Interpreter','latex','FontSize',14)
xlabel('time step, k','FontSize',14)
title('NIS Estimation Results','FontSize',14)
legend('NIS @ time k', 'r_1 bound', 'r_2 bound'),grid on
% ylim([r1y-2 r2y+2])

% process noise test
% X_truthWithNoise = mvnrnd(X_truth', Q)';
% X_truthAtObservationEpochs = X_truthWithNoise(:, (1 + integrationObservationRatio):integrationObservationRatio:end);
% Y_simulatedProcessNoSignalNoise = measurement.SimulateYData(X_truthAtObservationEpochs, landmarkPositions, data.R_CtoN, zeros(2), delT_observation);
% 
% [X_estimateProcessNoise, P_processNoise] = ExtendedKalmanFilter.Filter(X_0, P_0, Q, Y_simulatedProcessNoSignalNoise, R, delT_observation, landmarkPositions, data.R_CtoN);
% utilities.confidencePlots(T, X_estimateProcessNoise, P_processNoise, 2, "State Estimate, process noise", XLabels, XUnits)
% utilities.errorPlots(T(:, 2:end), X_truthAtObservationEpochs, X_estimateProcessNoise(:, 2:end), P_processNoise(:, :, 2:end), 2, "Error, process noise", XLabels, XUnits)
% 
% % signal noise test
% Y_simulatedProcessSignalNoise = [Y_simulatedProcessNoSignalNoise(1:2, :); mvnrnd(Y_simulatedProcessNoSignalNoise(3:4, :)', R)'];
% [X_estimateProcessSignalNoise, P_processSignalNoise] = ExtendedKalmanFilter.Filter(X_0, P_0, Q, Y_simulatedProcessNoSignalNoise, R, delT_observation, landmarkPositions, data.R_CtoN);
% utilities.confidencePlots(T, X_estimateProcessSignalNoise, P_processSignalNoise, 2, "State Estimate, process & signal noise", XLabels, XUnits)
% utilities.errorPlots(T(:, 2:end), X_truthAtObservationEpochs, X_estimateProcessSignalNoise(:, 2:end), P_processSignalNoise(:, :, 2:end), 2, "Error, process & signal noise", XLabels, XUnits)



% % :::::Plots:::::
% 
% time = (0:delT_observation:t_end+600)/3600;
% NLtime = (0:delT_integration:(t_end+660))/3600;
% 
% NL_state = X_truth;
% noisy_NL_state = X_noisy;
% filt_total_state = X_estimate;
% P_P = P_estimate;
% 
% % plot typical noisy traj vs noiseless traj
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
% 
% % % plot noisy simulated data
% % figure()
% % for ii = 1:5
% %     meastime = y_table(find(y_table(:,2)==ii),1);
% %     noisy_u = y_table(find(y_table(:,2)==ii),3);
% %     noisy_v = y_table(find(y_table(:,2)==ii),4);
% %     noiseless_u = y_noiseless(find(y_noiseless(:,2)==ii),3);
% %     noiseless_v = y_noiseless(find(y_noiseless(:,2)==ii),4);
% %     plot(meastime/3600,noisy_u,'x')
% %     hold on
% %     plot(meastime/3600,noiseless_u,'o')
% % end
% 
% % plot state estimation errors
% state_est_error = filt_total_state - noisy_NL_state(:,1:10:end);
% 
% 
% figure()
% sgtitle('Position state estimation error')
% subplot(311)
% plot(time,state_est_error(1,:))
% hold on
% plot(time,2*sqrt(reshape(P_P(1,1,:),[],1)),'black --')
% plot(time,-2*sqrt(reshape(P_P(1,1,:),[],1)),'black --')
% %ylim([-1e-4 1e-4])
% 
% subplot(312)
% plot(time,state_est_error(2,:))
% hold on
% plot(time,2*sqrt(reshape(P_P(2,2,:),[],1)),'black --')
% plot(time,-2*sqrt(reshape(P_P(2,2,:),[],1)),'black --')
% %ylim([-1e-3 1e-3])
% 
% subplot(313)
% plot(time,state_est_error(3,:))
% hold on
% plot(time,2*sqrt(reshape(P_P(3,3,:),[],1)),'black --')
% plot(time,-2*sqrt(reshape(P_P(3,3,:),[],1)),'black --')
% %ylim([-1e-3 1e-3])
% 
% ybounds = [-.02 .02];
% figure()
% sgtitle('Velocity state estimation error')
% subplot(311)
% plot(time,state_est_error(4,:))
% hold on
% plot(time,2*sqrt(reshape(P_P(4,4,:),[],1)),'black --')
% plot(time,-2*sqrt(reshape(P_P(4,4,:),[],1)),'black --')
% ylim(ybounds)
% 
% subplot(312)
% plot(time,state_est_error(5,:))
% hold on
% plot(time,2*sqrt(reshape(P_P(5,5,:),[],1)),'black --')
% plot(time,-2*sqrt(reshape(P_P(5,5,:),[],1)),'black --')
% ylim(ybounds)
% 
% subplot(313)
% plot(time,state_est_error(6,:))
% hold on
% plot(time,2*sqrt(reshape(P_P(6,6,:),[],1)),'black --')
% plot(time,-2*sqrt(reshape(P_P(6,6,:),[],1)),'black --')
% ylim(ybounds)
% 
% 
% 
% %%DO NEES TEST:
% epsNEESbar = mean(NEES_samps,1);
% alphaNEES = 0.05; %%significance level
% Nnx = k*n;
% %%compute intervals:
% r1x = chi2inv(alphaNEES/2, Nnx )./ Nsimruns
% r2x = chi2inv(1-alphaNEES/2, Nnx )./ Nsimruns
% 
% figure()
% plot(epsNEESbar,'ro','MarkerSize',6,'LineWidth',2),hold on
% plot(r1x*ones(size(epsNEESbar)),'r--','LineWidth',2)
% plot(r2x*ones(size(epsNEESbar)),'r--','LineWidth',2)
% ylabel('NEES statistic, $\bar{\epsilon}_x$','Interpreter','latex', 'FontSize',14)
% xlabel('time step, k','FontSize',14)
% title('NEES Estimation Results','FontSize',14)
% legend('NEES @ time k', 'r_1 bound', 'r_2 bound'),grid on
% % ylim([r1x-2 r2x+2])
% 
% 
% %%DO NIS TEST:
% epsNISbar = mean(NIS_samps,1);
% alphaNIS = 0.05; 
% Nny = Nsimruns*3;
% %%compute intervals:
% r1y = chi2inv(alphaNIS/2, Nny )./ Nsimruns
% r2y = chi2inv(1-alphaNIS/2, Nny )./ Nsimruns
% 
% figure()
% plot(epsNISbar,'bo','MarkerSize',6,'LineWidth',2),hold on
% plot(r1y*ones(size(epsNISbar)),'b--','LineWidth',2)
% plot(r2y*ones(size(epsNISbar)),'b--','LineWidth',2)
% ylabel('NIS statistic, $\bar{\epsilon}_y$','Interpreter','latex','FontSize',14)
% xlabel('time step, k','FontSize',14)
% title('NIS Estimation Results','FontSize',14)
% legend('NIS @ time k', 'r_1 bound', 'r_2 bound'),grid on
% % ylim([r1y-2 r2y+2])
