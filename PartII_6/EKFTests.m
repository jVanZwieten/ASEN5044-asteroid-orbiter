classdef EKFTests
    properties(Constant)
        data = load("orbitdetermination-finalproj_data_2023_11_14.mat");
        W = 1.4e-11; % tune on this property
        dX_0 = [1e-5 1e-5 1e-5 1e-7 1e-7 1e-7]';
        MonteCarloIterations = 10;
        
    end
    
    methods(Static)
        function setup()
            clc
            clear
            close all
            format longg
            rng(117)
            setGlobalVariables()
            
            global delT_observation delT_integration t_end w_A T integrationObservationDtRatio landmarkPositions Qkf
            T = (-600:delT_observation:t_end)/60/60;
            integrationObservationDtRatio = delT_observation/delT_integration;
            
            landmarkPositions = propagateLandmarksInInertialFrame(EKFTests.data.pos_lmks_A, delT_observation, t_end, w_A);
            Qkf = DTsys.noiseMat(EKFTests.W^2, delT_observation);
        end
        
        function X_observe = generateXTrueObserved(X_0, noise)
            global delT_integration delT_observation t_end
            if nargin > 1 && noise
                gamma = [zeros(3); eye(3)];
            else
                gamma = zeros(6,3);
            end
            
            X_integration = numerical.propagate(X_0, delT_integration, t_end + delT_observation, gamma);
            X_observe = EKFTests.XToObservedTime(X_integration);
        end
        
        function X_observe = XToObservedTime(X_integration)
            global integrationObservationDtRatio
            X_observe = X_integration(:, (1 + integrationObservationDtRatio):integrationObservationDtRatio:end);
        end
        
        function Y = YFromX(X, noise)
            global R delT_observation landmarkPositions
            [Y] = measurement.SimulateYData(X, landmarkPositions, EKFTests.data.R_CtoN, R, delT_observation);
            if nargin > 1 && noise
                Y = mvnrnd(Y', R)';
            end
        end
        
        function runEKFAndReport(X_truth, X_0, P_0, Y, title)
            global delT_observation t_end w_A XLabels XUnits R
            
            T = (-600:delT_observation:t_end)/60/60;
            data = load("orbitdetermination-finalproj_data_2023_11_14.mat");
            landmarkPositions = propagateLandmarksInInertialFrame(data.pos_lmks_A, delT_observation, t_end, w_A);
            W = 3e-12; % tune on this property
            Qkf = DTsys.noiseMat(W^2, delT_observation);
            
            [X_estimate, P_estimate] = ExtendedKalmanFilter.Filter(X_0, P_0, Qkf, Y, R, delT_observation, landmarkPositions, data.R_CtoN);
            utilities.confidencePlots(T, X_estimate, P_estimate, 2, "State Estimate, " + title, XLabels, XUnits)
            utilities.errorPlots(T(:, 2:end), X_truth, X_estimate(:, 2:end), P_estimate(:, :, 2:end), 2, "Error, " + title, XLabels, XUnits)
        end
        
        function truthModelTest(X_truths, X_0, P_0, Ys, title, boundNeesNis)
            global Qkf R delT_observation landmarkPositions T XLabels XUnits
            monteCarloIterations = size(X_truths, 3);
            totalSteps = size(X_truths, 2);
            NEES = zeros(monteCarloIterations, totalSteps);
            NIS = zeros(monteCarloIterations, totalSteps);
            
            for i = 1:monteCarloIterations
                [X, P, NIS(i,:)] = ExtendedKalmanFilter.Filter(X_0, P_0, Qkf, Ys{i}, R, delT_observation, landmarkPositions, EKFTests.data.R_CtoN);
                NEES(i,:) = utilities.calculateNEES(X_truths(:, :, i), X(:, 2:end), P(:, :, 2:end));
            end
            
            utilities.confidencePlots(T, X, P, 2, "State Estimate, " + title, XLabels, XUnits)
            utilities.errorPlots(T(:, 2:end), X_truths(:, :, end), X(:, 2:end), P(:, :, 2:end), 2, "Error, " + title, XLabels, XUnits)
            
            if nargin < 6
                boundNeesNis = false;
            end
            utilities.NEESPlot(NEES, title, boundNeesNis)
            utilities.NISPlot(NIS, title, boundNeesNis)
        end
        
        function noiseless()
            EKFTests.setup();
            global X_0 P_0
            
            X_noiseless = EKFTests.generateXTrueObserved(X_0, false);
            Y_noiseless = EKFTests.YFromX(X_noiseless);
            EKFTests.runEKFAndReport(X_noiseless, X_0, P_0, Y_noiseless, "Noiseless")
        end
        
        function perturbation()
            EKFTests.setup();
            global X_0  P_0
            
            X_perturbedObserved = EKFTests.generateXTrueObserved(X_0 + EKFTests.dX_0, false);
            [Y_purturbed] = EKFTests.YFromX(X_perturbedObserved);
            EKFTests.runEKFAndReport(X_perturbedObserved, X_0, P_0, Y_purturbed, "Perturbed")
        end
        
        function processNoise()
            EKFTests.setup();
            global X_0 P_0
            
            [X_processNoiseObserved, Y_processNoiseObserved] = EKFTests.MonteCarloTruthSimulations(X_0, EKFTests.MonteCarloIterations);
            EKFTests.truthModelTest(X_processNoiseObserved, X_0, P_0, Y_processNoiseObserved, "Process Noise")
        end

        function [X, Y] = generateProcessNoiseSim(X_0, signalNoise)
            if nargin < 2
                signalNoise = false;
            end

            X = EKFTests.generateXTrueObserved(X_0, true);
            Y = EKFTests.YFromX(X, signalNoise);
        end

        function [X, Y] = MonteCarloTruthSimulations(X_0, iterations, signalNoise)
            if nargin < 3
                signalNoise = false;
            end

            [X_processNoiseObserved0, Y_processNoiseObserved0] = EKFTests.generateProcessNoiseSim(X_0, signalNoise);
            X = zeros(size(X_processNoiseObserved0, 1), size(X_processNoiseObserved0, 2), iterations);
            X(:, :, 1) = X_processNoiseObserved0;
            Y = cell(iterations);
            Y{1} = Y_processNoiseObserved0;
            
            for i = 2:iterations
                [X(:, :, i), Y{i}] = EKFTests.generateProcessNoiseSim(X_0, signalNoise);
            end
        end

        function processNoisePerturbation()
            EKFTests.setup()
            global X_0 P_0

            [X_processNoiseObserved, Y_processNoiseObserved] = EKFTests.MonteCarloTruthSimulations(X_0 + EKFTests.dX_0, EKFTests.MonteCarloIterations);
            EKFTests.truthModelTest(X_processNoiseObserved, X_0, P_0, Y_processNoiseObserved, "Process Noise")
        end

        function fullNoise()
            EKFTests.setup()
            global X_0 P_0 Qkf R T landmarkPositions delT_observation

            [X_truth, Y] = EKFTests.MonteCarloTruthSimulations(X_0, EKFTests.MonteCarloIterations);

            X_truthTypical = X_truth(:, :, 1);
            Y_typical = Y{1};

            Y_noisy = cell(length(Y));
            for i = 1:length(Y)
                Y_noisy{i} = Y{i};
                Y_noisy{i}(3:4, :) = mvnrnd(Y{i}(3:4, :)', R)';
            end

            X_Noiseless = EKFTests.generateXTrueObserved(X_0);
            utilities.plotSimulatedTrajectories(T(2:end), X_Noiseless, X_truthTypical)
            
            Y_noisyTypical = Y_noisy{1};
            utilities.plotSimulatedYNoise(Y_typical,  Y_noisyTypical)

            monteCarloIterations = size(X_truth, 3);
            totalSteps = size(X_truth, 2);
            NEES = zeros(monteCarloIterations, totalSteps);
            NIS = zeros(monteCarloIterations, totalSteps);
            
            for i = 1:monteCarloIterations
                [X, P, NIS(i,:)] = ExtendedKalmanFilter.Filter(X_0, P_0, Qkf, Y_noisy{i}, R, delT_observation, landmarkPositions, EKFTests.data.R_CtoN);
                NEES(i,:) = utilities.calculateNEES(X_truth(:, :, i), X(:, 2:end), P(:, :, 2:end));
            end
            utilities.plotStateEstimationErrors(T(2:end), X_truth(:, :, end), X(:, 2:end), P(:, :, 2:end))

            utilities.NEESPlot(NEES, "")
            utilities.NISPlot(NIS, "")
        end
    end
end