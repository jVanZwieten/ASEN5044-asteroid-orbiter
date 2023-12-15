clc
clear
close all
format longg
setGlobalVariables()

global X_0 XLabels XUnits P_0 delT_integration delT_observation t_end w_A p sigma_w sigma_u
rng(117)

% no noise test
data = load("orbitdetermination-finalproj_data_2023_11_14.mat");
X_truthNoProcessNoise = numerical.propagate(X_0, delT_integration, t_end + delT_observation); % +delT_observation accounts for  t = 0
Q = DTsys.noiseMat(sigma_w^2, delT_observation);
utilities.plot3D(X_truthNoProcessNoise(1:3, :), [-1 1; -1 1; -1 1], "no process noise orbit")

integrationObservationRatio = delT_observation/delT_integration;
landmarkPositions = propagateLandmarksInInertialFrame(data.pos_lmks_A, delT_observation, t_end, w_A);
X_truthAtObservationEpochs = X_truthNoProcessNoise(:, (1 + integrationObservationRatio):integrationObservationRatio:end);
Y_simulatedNoSignalNoise = measurement.SimulateYData(X_truthAtObservationEpochs, landmarkPositions, data.R_CtoN, zeros(2), delT_observation);
R = sigma_u^2*eye(p);

[X_estimateNoNoise, P_noNoise] = ExtendedKalmanFilter.Filter(X_0, P_0, Q, Y_simulatedNoSignalNoise, R, delT_observation, landmarkPositions, data.R_CtoN);

T = (-600:delT_observation:t_end)/60/60;
utilities.confidencePlots(T, X_estimateNoNoise, P_noNoise, 2, "State Estimate, no noise", XLabels, XUnits)

utilities.errorPlots(T(:, 2:end), X_truthAtObservationEpochs, X_estimateNoNoise(:, 2:end), P_noNoise(:, :, 2:end), 2, "Error, no noise", XLabels, XUnits)


% process noise test
X_truthWithNoise = mvnrnd(X_truthNoProcessNoise', Q)';
X_truthAtObservationEpochs = X_truthWithNoise(:, (1 + integrationObservationRatio):integrationObservationRatio:end);
Y_simulatedProcessNoSignalNoise = measurement.SimulateYData(X_truthAtObservationEpochs, landmarkPositions, data.R_CtoN, zeros(2), delT_observation);

[X_estimateProcessNoise, P_processNoise] = ExtendedKalmanFilter.Filter(X_0, P_0, Q, Y_simulatedProcessNoSignalNoise, R, delT_observation, landmarkPositions, data.R_CtoN);
utilities.confidencePlots(T, X_estimateProcessNoise, P_processNoise, 2, "State Estimate, process noise", XLabels, XUnits)
utilities.errorPlots(T(:, 2:end), X_truthAtObservationEpochs, X_estimateProcessNoise(:, 2:end), P_processNoise(:, :, 2:end), 2, "Error, process noise", XLabels, XUnits)

% signal noise test
Y_simulatedProcessSignalNoise = [Y_simulatedProcessNoSignalNoise(1:2, :); mvnrnd(Y_simulatedProcessNoSignalNoise(3:4, :)', R)'];
[X_estimateProcessSignalNoise, P_processSignalNoise] = ExtendedKalmanFilter.Filter(X_0, P_0, Q, Y_simulatedProcessNoSignalNoise, R, delT_observation, landmarkPositions, data.R_CtoN);
utilities.confidencePlots(T, X_estimateProcessSignalNoise, P_processSignalNoise, 2, "State Estimate, process & signal noise", XLabels, XUnits)
utilities.errorPlots(T(:, 2:end), X_truthAtObservationEpochs, X_estimateProcessSignalNoise(:, 2:end), P_processSignalNoise(:, :, 2:end), 2, "Error, process & signal noise", XLabels, XUnits)
