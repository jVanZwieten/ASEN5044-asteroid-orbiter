function landmarkPositions = propagateLandmarksInInertialFrame(landmarkPositions_initial, delT, t_end, w_A)
    totalSteps = t_end/delT + 1; % +1 to account for t = 0

    numberLandmarks = size(landmarkPositions_initial, 2);

    landmarkPositions = zeros(5, numberLandmarks*totalSteps); % this does not include initial state; assuming t = 0s has propagated 1 step

    i = 1;
    for step = 1:totalSteps
        rotation = w_A*((step-1)*delT);
        R_asteroidToInertial = [
            cos(rotation) -sin(rotation) 0
            sin(rotation) cos(rotation) 0
            0 0 1];

        for landmarkIndex = 1:numberLandmarks
            landmarkPositions(:, i) = [
                (step - 1)*delT
                landmarkIndex
                R_asteroidToInertial*landmarkPositions_initial(:, landmarkIndex)];

            i = i + 1; % Matlab doesn't have an increment operator... *sigh*
        end
    end
end