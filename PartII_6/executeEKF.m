clc
clear
close all
format longg
setGlobalVariables()

global X_0 delT_integration delT_observation t_end w_A
P0 = [10/1000*eye(3) zeros(3)
      zeros(3) 0.5/1000/1000*eye(3)];

data = load("orbitdetermination-finalproj_data_2023_11_14.mat");
landmarkPositions = propagateLandmarksInInertialFrame(data.pos_lmks_A, delT_observation, t_end, w_A);
X_truthNoProcessNoise = numerical.propagate(X_0, delT_integration, t_end);

integrationObservationRatio = delT_observation/delT_integration;
Y_simulatedNoSignalNoise = measurement.SimulateYData(X_truthNoProcessNoise(:, 2:integrationObservationRatio:end), landmarkPositions, data.R_CtoN, zeros(2), delT_observation);