global mu_A T_rotationA w_A r_sa phi_0 rho Am delT_integration delT_observation t_end r_0 rDot_0 X_0 XLabels XUnits P_0 n Q sigma_w f_camera u_0 v_0 u_min v_min u_max v_max sigma_u sigma_v p R
mu_A = 4.892E-9;                        % km^3/s^2
T_rotationA = 4.296057;                 % h
w_A = 2*pi/(T_rotationA*60*60);         % rad/s

r_sa = [1.5E8 0 0]';                    % km
phi_0 = 1E14;                           % kg*km/s^2
rho = 0.4;                              % coefficient of reflectivity
Am = (1/62)*10^(-6);                    % km^2/kg, area to mass ratio

delT_integration = 60;                  % s
delT_observation = 600;                 % s
t_end = 72*60*60;                       % 72h -> s

r_0 = [0 -1 0]';                        % km
rDot_0 = [0 0 sqrt(mu_A/norm(r_0))]';   % km/s
X_0 = [r_0; rDot_0];                    % initial state
XLabels = ["x", "y", "z", "v_x", "v_y", "v_z"];
XUnits = ["km", "km", "km", "km/s", "km/s", "km/s"];
P_0 = [0.01/(1000)^2*eye(3) zeros(3)      % initial state covarience
      zeros(3) 0.001/(1e6)^2*eye(3)];
sigma_w = 1e-12;                         % km/s^2, process noise standard deviation
n = size(X_0, 1);
Q = DTsys.noiseMat(sigma_w^2, delT_observation);

f_camera = 2089.7959;                   % pixels, camera focal length
u_0 = 512;                              % pixels
u_min = 0;                              % pixels
u_max = 1024;                           % pixels
sigma_u = 0.25;                         % pixels^2
v_0 = 512;                              % pixels
v_min = 0;                              % pixels
v_max = 1024;                           % pixels
sigma_v = 0.25;                         % pixels^2
p = 2;
R = sigma_u^2*eye(p);