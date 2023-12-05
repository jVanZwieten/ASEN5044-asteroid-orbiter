%%% Final Project - Orbit Determination

clc, clear
% load in data from canvas
load("C:\Users\jared\MATLAB Drive\ASEN5044\Project" + ...
    "\orbitdetermination-finalproj_data_2023_11_14.mat")

% % plot landmarks in 3D (uncomment to see plot)
% figure(1)
% scatter3(pos_lmks_A(1,:),pos_lmks_A(2,:),pos_lmks_A(3,:))

% Constants
mu_A = 4.892e-9; % gravitational parameter of asteroid (km^3/s^2)
f_C = 2089.7959; % camera focal length (pixels)
sigma_u = 0.25;  % u-measurement std dev (pixels)
sigma_v = 0.25;  % v-measurement std dev (pixels)
u_max = 1024;    % max pixel coordinates
v_max = 1024;    % max pixel coordinates
u0 = 512;        % principal point coordinates
v0 = 512;        % principal point coordinates
omega_A = 2*pi/(4.296057*3600); % asteroid rotation rate



%%% Part I - nonlinear dynamics propagation

% Set up nominal initial state
pos_0 = [0 -1 0]';                    % position vector (km)
vel_0 = [0 0 sqrt(mu_A/norm(pos_0))]';   % velocity vector (km/s)
x0 = [pos_0; vel_0];                     % state vector

% Define ode45 integration time span
dt = 600;                  % this can be 60s or 600s (measurement dt)
t_end = y_table(end,1);    % last time in measurement data
tspan = 0:dt:t_end;

% solve NL dynamics using ode45 and plot
[t_NL, x_NL] = ode45(@Osiris_EOM,tspan,x0);
plot3(x_NL(:,1),x_NL(:,2),x_NL(:,3),'LineWidth',1.5);
hold on
plot3(0,0,0,'.','MarkerSize',20)
grid on
xlabel('X (km)')
ylabel('Y (km)')
zlabel('Z (km)')


%%% using NL measurement model, generate perfect measurements to compare
r = x_NL(:,1:3);

% Convert x_NL to the rotating asteroid frame to plot and compare with 
% plot given in project instructions

for t = 1:length(tspan)
    theta = omega_A*(t-1)*dt;
    R_AtoN(:,:,t)= [cos(theta) -sin(theta) 0;
                   sin(theta)  cos(theta) 0;
                     0            0      1];
end


% should R_AtoN be transposed to convert from inertial to rotating??? 


y_sim=[];

for ind = 1:1 % what time point you're checking at
    i_C = R_CtoN(:,1,ind);
    j_C = R_CtoN(:,2,ind);
    k_C = R_CtoN(:,3,ind);
    r = R_AtoN(:,:,ind)*x_NL(ind,1:3)';

    % determine which landmarks are in view and in front
    for jj = 1:11  % which landmarks to check if in FOV
        l = pos_lmks_A(:,jj);
        u = f_C*(dot(l-r,i_C)/dot(l-r,k_C)) + u0;
        v = f_C*(dot(l-r,j_C)/dot(l-r,k_C)) + v0;
        landmark_in_FOV = (u>0 && u<u_max) && (v>0 && v<v_max) && (dot(l-r,k_C)>0);
       
        if landmark_in_FOV == true
           landmark_in_front = dot(l,k_C) < 0;
           if landmark_in_front == true
               % if landmark is in FOV and in front, store the time, 
               % landmark ID, pixel coordinates
               time = (ind-1)*dt;
               ID = jj;
               y_sim = [y_sim; time ID u v];
           end
        end
    end
end
