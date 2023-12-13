setGlobalVariables()
global delT_integration delT_observation w_A X_0
data = load("orbitdetermination-finalproj_data_2023_11_14.mat");

t_end = 3*60*60;
tSteps = t_end/delT_observation;
landmarkPositions = propagateLandmarksInInertialFrame(data.pos_lmks_A, delT_observation, t_end, w_A);

landmarks = 3;
colors = jet(tSteps);

figure
view(3)
hold on
for landmarkIndex = 1:size(data.pos_lmks_A, 2)
    landmarkOverTime = landmarkPositions(3:end, landmarkPositions(2, :) == landmarkIndex);
    for timeStep = 1:(size(landmarkOverTime, 2) - 1)
        plot3(landmarkOverTime(1, timeStep:(timeStep + 1)), landmarkOverTime(2, timeStep:(timeStep + 1)), landmarkOverTime(3, timeStep:(timeStep + 1)), ...
            'color', colors(timeStep, :))
    end
end

X = numerical.propagate(X_0, delT_integration, t_end);
colors = jet(size(X, 2));
for timeStep = 1:(size(X, 2) - 1)
    plot3(X(1, timeStep:(timeStep + 1)), X(2, timeStep:(timeStep + 1)), X(3, timeStep:(timeStep + 1)), ...
        'color', colors(timeStep, :))
end

xlim([-1 1])
ylim([-1 1])
zlim([-1 1])