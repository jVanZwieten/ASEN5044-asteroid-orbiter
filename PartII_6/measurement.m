classdef measurement
    methods(Static)
        function Y = h(X, landmarkPosition, rotation_cameraToInertial)
            global f_camera u_0 v_0
            orbiterPosition = X(1:3);
            i_c = rotation_cameraToInertial(:, 1);
            j_c = rotation_cameraToInertial(:, 2);
            k_c = rotation_cameraToInertial(:, 3);

            Y = zeros(2, 1);
            Y(1) = f_camera*((landmarkPosition - orbiterPosition)'*i_c)/((landmarkPosition - orbiterPosition)'*k_c) + u_0;
            Y(2) = f_camera*((landmarkPosition - orbiterPosition)'*j_c)/((landmarkPosition - orbiterPosition)'*k_c) + v_0;
        end

        function YStack = Y_epoch(X, landmarkPositions, rotation_cameraToInertial)
            global p

            YStack = zeros(2*size(landmarkPositions, 2), 1);
            for landmarkIndex = 1:size(landmarkPositions, 2)
                YStack((landmarkIndex*p - 1):landmarkIndex*p, 1) = measurement.h(X, landmarkPositions(:, landmarkIndex), rotation_cameraToInertial);
            end
        end

        function [Y, Y_noisy] = SimulateYData(X, landmarkPositionsThroughTime_inertial, rotations_cameraToInertial, R, dT)
            Y = zeros(4, size(landmarkPositionsThroughTime_inertial, 2));
            t_end = landmarkPositionsThroughTime_inertial(1, end);

            timeSteps = t_end/dT + 1; % +1 to account for t = 0

            i = 0;
            for k = 1:timeSteps
                t = (k - 1)*dT;
                landmarks_epoch = landmarkPositionsThroughTime_inertial(2:end, landmarkPositionsThroughTime_inertial(1, :) == t);
                rotation_camera = rotations_cameraToInertial(:, :, k);
                for landmarkIndex = 1:size(landmarks_epoch, 2)
                    X_k = X(:, k);
                    landmarkPosition = landmarks_epoch(2:end, landmarks_epoch(1, :) == landmarkIndex);
                    u = measurement.h(X_k, landmarkPosition, rotation_camera);

                    if(isVisible(u, landmarkPosition, X_k(1:3), rotation_camera(:, 3)))
                        y = [t
                            landmarkIndex
                            u];
    

                        i = i + 1;
                        Y(:, i) = y;
                    end
                end
            end
            Y = Y(:, 1:i);

            signalNoise = mvnrnd([0 0], R, size(Y, 2))';
            Y_noisy = Y + [zeros(2, size(Y, 2)); signalNoise];
        end
    end
end