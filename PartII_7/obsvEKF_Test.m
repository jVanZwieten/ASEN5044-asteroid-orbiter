clc
clear
close all
format longg
setGlobalVariables()



data = load("orbitdetermination-finalproj_data_2023_11_14.mat");
NL_state = numerical.propagate(X_0, delT_integration, t_end);
Q = DTsys.noiseMat(sigma_w^2, delT_observation);

integrationObservationRatio = delT_observation/delT_integration;
landmarkPositions = propagateLandmarksInInertialFrame(data.pos_lmks_A, delT_observation, t_end, w_A);
X_truthAtObservationEpochs = NL_state(:, 1:integrationObservationRatio:end);
Y_given = data.y_table';
R = sigma_u^2*eye(p);

[filt_total_state, P_P] = ExtendedKalmanFilter.Filter(X_0, P_0, Q, Y_given, R, delT_observation, landmarkPositions, data.R_CtoN);
filt_total_state = filt_total_state(:,1:end-1);

state_est_error = X_truthAtObservationEpochs - filt_total_state;


% :::::Plots:::::

time = (0:delT_observation:t_end)/3600;
NLtime = (0:delT_integration:(t_end+60))/3600;


% plot state estimation errors
figure()
sgtitle('Position state estimation error')
subplot(311)
plot(time,state_est_error(1,:))
hold on
plot(time,2*sqrt(reshape(P_P(1,1,1:end-1),[],1)),'black --')
plot(time,-2*sqrt(reshape(P_P(1,1,1:end-1),[],1)),'black --')
%ylim([-1e-4 1e-4])

subplot(312)
plot(time,state_est_error(2,:))
hold on
plot(time,2*sqrt(reshape(P_P(2,2,1:end-1),[],1)),'black --')
plot(time,-2*sqrt(reshape(P_P(2,2,1:end-1),[],1)),'black --')
%ylim([-1e-3 1e-3])

subplot(313)
plot(time,state_est_error(3,:))
hold on
plot(time,2*sqrt(reshape(P_P(3,3,1:end-1),[],1)),'black --')
plot(time,-2*sqrt(reshape(P_P(3,3,1:end-1),[],1)),'black --')
%ylim([-1e-3 1e-3])

figure()
sgtitle('Velocity state estimation error')
subplot(311)
plot(time,state_est_error(4,:))
hold on
plot(time,2*sqrt(reshape(P_P(4,4,1:end-1),[],1)),'black --')
plot(time,-2*sqrt(reshape(P_P(4,4,1:end-1),[],1)),'black --')
%ylim([-1e-4 1e-4])

subplot(312)
plot(time,state_est_error(5,:))
hold on
plot(time,2*sqrt(reshape(P_P(5,5,1:end-1),[],1)),'black --')
plot(time,-2*sqrt(reshape(P_P(5,5,1:end-1),[],1)),'black --')
%ylim([-1e-3 1e-3])

subplot(313)
plot(time,state_est_error(6,:))
hold on
plot(time,2*sqrt(reshape(P_P(6,6,1:end-1),[],1)),'black --')
plot(time,-2*sqrt(reshape(P_P(6,6,1:end-1),[],1)),'black --')
%ylim([-1e-3 1e-3])



% plot nonlinear and estimated state

figure()
sgtitle('Position')
subplot(311)
plot(time,NL_state(1,1:10:end))
hold on
plot(time,filt_total_state(1,:))
plot(time,2*sqrt(reshape(P_P(1,1,1:end-1),[],1))+filt_total_state(1,:)','black --')
plot(time,-2*sqrt(reshape(P_P(1,1,1:end-1),[],1))+filt_total_state(1,:)','black --')
%ylim([-1e-4 1e-4])

subplot(312)
plot(time,NL_state(2,1:10:end))
hold on
plot(time,filt_total_state(2,:))
plot(time,2*sqrt(reshape(P_P(2,2,1:end-1),[],1))+filt_total_state(2,:)','black --')
plot(time,-2*sqrt(reshape(P_P(2,2,1:end-1),[],1))+filt_total_state(2,:)','black --')
%ylim([-1e-3 1e-3])

subplot(313)
plot(time,NL_state(3,1:10:end))
hold on
plot(time,filt_total_state(3,:))
plot(time,2*sqrt(reshape(P_P(3,3,1:end-1),[],1))+filt_total_state(3,:)','black --')
plot(time,-2*sqrt(reshape(P_P(3,3,1:end-1),[],1))+filt_total_state(3,:)','black --')
%ylim([-1e-3 1e-3])

figure()
sgtitle('Velocity')
subplot(311)
plot(time,NL_state(4,1:10:end))
hold on
plot(time,filt_total_state(4,:))
plot(time,2*sqrt(reshape(P_P(4,4,1:end-1),[],1))+filt_total_state(4,:)','black --')
plot(time,-2*sqrt(reshape(P_P(4,4,1:end-1),[],1))+filt_total_state(4,:)','black --')
%ylim([-1e-4 1e-4])

subplot(312)
plot(time,NL_state(5,1:10:end))
hold on
plot(time,filt_total_state(5,:))
plot(time,2*sqrt(reshape(P_P(5,5,1:end-1),[],1))+filt_total_state(5,:)','black --')
plot(time,-2*sqrt(reshape(P_P(5,5,1:end-1),[],1))+filt_total_state(5,:)','black --')
%ylim([-1e-3 1e-3])

subplot(313)
plot(time,NL_state(6,1:10:end))
hold on
plot(time,filt_total_state(6,:))
plot(time,2*sqrt(reshape(P_P(6,6,1:end-1),[],1))+filt_total_state(6,:)','black --')
plot(time,-2*sqrt(reshape(P_P(6,6,1:end-1),[],1))+filt_total_state(6,:)','black --')
%ylim([-1e-3 1e-3])

