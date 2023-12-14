clc
clear
close all
format longg
setGlobalVariables()

global X_0 P_0 delT_integration delT_observation t_end w_A p sigma_w sigma_u


data = load("orbitdetermination-finalproj_data_2023_11_14.mat");
X_truthNoProcessNoise = numerical.propagate(X_0, delT_integration, t_end);
Q = DTsys.noiseMat(sigma_w^2, delT_observation);

integrationObservationRatio = delT_observation/delT_integration;
landmarkPositions = propagateLandmarksInInertialFrame(data.pos_lmks_A, delT_observation, t_end, w_A);
X_truthAtObservationEpochs = X_truthNoProcessNoise(:, (1 + integrationObservationRatio):integrationObservationRatio:end);
Y_simulatedNoSignalNoise = measurement.SimulateYData(X_truthAtObservationEpochs, landmarkPositions, data.R_CtoN, zeros(2), delT_observation);
R = sigma_u^2*eye(p);

[X_estimateNoNoise, P_noNoise] = ExtendedKalmanFilter.Filter(X_0, P_0, Q, Y_simulatedNoSignalNoise, R, delT_observation, landmarkPositions, data.R_CtoN);

e = X_truthAtObservationEpochs(:, 2:end) - X_estimateNoNoise;