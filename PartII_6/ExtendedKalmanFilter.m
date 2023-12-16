classdef ExtendedKalmanFilter
    methods(Static)
        function [X, P, NIS_hist] = Filter(X_initial, P_initial, Q, Y, R, dT, landmarkPositionsThroughTime_inertial, rotations_cameraToInertial)
            global n
            assert(size(X_initial, 1) == n && size(X_initial, 2) == 1)
            assert(size(P_initial, 1) == n && size(P_initial, 2) == n && size(P_initial, 3) == 1)

            totalSteps = Y(1, end)/dT + 1; % +1 accounts for t=0

            X = zeros(n, totalSteps + 1);
            X(:, 1) = X_initial; % X_0
            P = zeros(n, n, totalSteps + 1);
            P(:, :, 1) = P_initial;

            NIS_hist = zeros(1,totalSteps);

            for k = 2:(totalSteps + 1)
                t = (k - 2)*dT; % t = 0 corresponds to X(:, 2), X_1, dT seconds after X_0
                Y_epoch = Y(2:end, Y(1, :) == t);
                correspondingLandmarks = ExtendedKalmanFilter.correspondingLandmarks(landmarkPositionsThroughTime_inertial, Y_epoch(1, :), t);
                Y_epoch = reshape(Y_epoch(2:3, :), 1, [])';
                rotation_cameraToInertial = rotations_cameraToInertial(:, :, k - 1);

                [X(:, k), P(:, :, k), NIS] = ExtendedKalmanFilter.propagateExtendedKalmanFilter(X(:, k - 1), P(:, :, k - 1), Q, dT, Y_epoch, R, correspondingLandmarks, rotation_cameraToInertial);
                NIS_hist(k - 1) = NIS/length(correspondingLandmarks);
            end
        end
        
        function landmarkPositions = correspondingLandmarks(allLandmarkPositions, landmarkIndicies, t)
            landmarkPositions = zeros(3, length(landmarkIndicies));
            
            landmarks_epoch = allLandmarkPositions(2:end, allLandmarkPositions(1, :) == t);
            for i = 1:length(landmarkIndicies)
                landmarkPositions(:, i) = landmarks_epoch(2:end, landmarks_epoch(1, :) == landmarkIndicies(i));
            end
        end

        function [Xestimate_k, P_k, NEES_k, NIS_k] = propagateExtendedKalmanFilter(X_kPrevious, P_kPrevious, omegaQomega, dT, Y_k, R, landmarkPositions, rotation_cameraToInertial)
            global n

            gammaW = [zeros(3); eye(3)];

            XinitialEstimate_k = numerical.rk4_state(X_kPrevious, dT, gammaW);

            Ftilde_k = eye(n) + dT*CTsys.AEvaluated(X_kPrevious);
            P_kInitial = Ftilde_k*P_kPrevious*Ftilde_k' + omegaQomega;

            Yestimate = measurement.Y_epoch(XinitialEstimate_k, landmarkPositions, rotation_cameraToInertial);
            Htilde_k = ExtendedKalmanFilter.Htilde(XinitialEstimate_k, landmarkPositions, rotation_cameraToInertial);
            error_k = Y_k - Yestimate;

            R_k = kron(eye(size(landmarkPositions, 2)), R);
            K_k = P_kInitial*Htilde_k'*(Htilde_k*P_kInitial*Htilde_k' + R_k)^-1;

            Xestimate_k = XinitialEstimate_k + K_k*error_k;
            P_k = (eye(n) - K_k*Htilde_k)*P_kInitial;
            chol(P_k);  % Should crash if we do something dumb

            ex = -Xestimate_k;
            ey = error_k;
            NEES_k = ex'*pinv(P_k)*ex;
            NIS_k = ey'*pinv(Htilde_k*P_kInitial*Htilde_k'+R_k)*ey;
        end
        
        function Htilde_k = Htilde(X_estimate, landmarkPositions, rotation_cameraToInertial)
            global n p
            Htilde_k = zeros(2*size(landmarkPositions, 2), n);
            for i = 1:size(landmarkPositions, 2)
                landmarkPosition = landmarkPositions(:, i);
                Htilde_k((i*p - 1):(i*p), 1:6) = CTsys.CEvaluated(X_estimate, landmarkPosition, rotation_cameraToInertial);
            end
        end
    end
end