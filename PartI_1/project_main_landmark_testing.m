%%% Final Project - Orbit Determination

close all
clc, clear

% load in data from canvas
%load("C:\Users\jared\MATLAB Drive\ASEN5044\Project" + ...
%    "\orbitdetermination-finalproj_data_2023_11_14.mat")
data = load("orbitdetermination-finalproj_data_2023_11_14.mat");

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
pos_0 = [0 -1 0]';                       % position vector (km)
vel_0 = [0 0 sqrt(mu_A/norm(pos_0))]';   % velocity vector (km/s)
x0 = [pos_0; vel_0];                     % state vector

% Define ode45 integration time span
dt = 60;          % this can be 60s or 600s (measurement dt)
t_end = 5*86400;  % 5 days
DCO = 3*86400;    % data cutoff (DCO) = 3 days
tspan = 0:dt:t_end;
meas_tspan = 0:600:DCO;

% % solve NL dynamics using ode45 and plot
opts = odeset('RelTol',1e-12,'AbsTol',1e-12);
[t_NL, x_NL] = ode45(@Osiris_EOM,tspan,x0,opts);
% figure()
% plot3(x_NL(:,1),x_NL(:,2),x_NL(:,3)) %,'LineWidth',1.5);
% hold on
% plot3(0,0,0,'.','MarkerSize',20)
% grid on
% xlabel('X (km)')
% ylabel('Y (km)')
% zlabel('Z (km)')

% plot states vs time
for s = 1:6
    states = ["x (km)","y (km)","z (km)","xdot (km/s)","ydot (km/s)","zdot (km/s)"];
    subplot(6,1,s)
    plot(tspan/3600,x_NL(:,s),'LineWidth',2)
    ylabel(states(s))
    xlabel('Time (hours)')
    sgtitle('Inertial States vs Time')
    hold on
end



% subplot(612)
% plot(tspan/3600,x_NL(:,2),'LineWidth',2)
% 
% subplot(612)
% plot(tspan/3600,x_NL(:,3),'LineWidth',2)
% 
% subplot(612)
% plot(tspan/3600,x_NL(:,4),'LineWidth',2)
% 
% subplot(612)
% plot(tspan/3600,x_NL(:,5),'LineWidth',2)
% 
% subplot(612)
% plot(tspan/3600,x_NL(:,6),'LineWidth',2)


%%% using NL measurement model, generate perfect measurements to compare
R_AtoN = zeros(3,3,length(tspan));
R_AtoN(:,:,1) = [1 0 0; 0 1 0; 0 0 1];
r_rotating(:,1) = x_NL(1,1:3)';

% Convert x_NL to the rotating asteroid frame to plot and compare with 
% plot given in project instructions

for t = 1:length(tspan)-1
    theta = omega_A*(t)*dt;
    R_AtoN(:,:,t+1)= [cos(theta) -sin(theta) 0;
                      sin(theta)  cos(theta) 0;
                          0           0      1];

    r_rotating(:,t+1) = R_AtoN(:,:,t+1)*x_NL(t+1,1:3)';
end

% figure()   
% plot3(r_rotating(1,:),r_rotating(2,:),r_rotating(3,:))

y_sim=[];

% why isnt it working :/
% use inertial position vector


for ind = 1:length(meas_tspan) % what time point you're checking at
    i_C = data.R_CtoN(:,1,ind);
    j_C = data.R_CtoN(:,2,ind);
    k_C = data.R_CtoN(:,3,ind);
    r = x_NL(10*(ind-1) + 1,1:3)'; % every tenth position

    % determine which landmarks are in view and in front
    for jj = 1:length(data.pos_lmks_A) % which landmarks to check if in FOV
        l = R_AtoN(:,:,10*(ind-1)+1)*data.pos_lmks_A(:,jj); % convert lmk to inertial
        u = f_C*(dot(l-r,i_C)/dot(l-r,k_C)) + u0;
        v = f_C*(dot(l-r,j_C)/dot(l-r,k_C)) + v0;
        landmark_in_FOV = (u>0 && u<u_max) && (v>0 && v<v_max) && (dot(l-r,k_C)>0);
       
        if landmark_in_FOV == true
           landmark_in_front = dot(l,k_C) < 0;
           if landmark_in_front == true
               % if landmark is in FOV and in front, store the time, 
               % landmark ID, pixel coordinates
               time = (10*(ind-1))*dt;
               ID = jj;
               y_sim = [y_sim; time ID u v];
           end
        end
    end
end

% plot measurements 

y_sim = sortrows(y_sim,2); % sort measurements by ID
lmk=[];
markers = {'o','+','*','.','x','s','d','v','^','>'};
for k = 1:5
    i=1; % reset i
    figure() % open new figure

    while i <= length(data.pos_lmks_A)

         % need to look for k sets of 10 landmarks and plot u,v values vs
         % time with legend
         ID_index = find(y_sim(:,2) == i+10*(k-1));
         first = ID_index(1);
         last = ID_index(end);
         u = y_sim(first:last,3);
         v = y_sim(first:last,4);
         lmk = [lmk; string(i+10*(k-1))];
         time = y_sim(first:last,1);
         subplot(211)
         plot(time/3600,u,markers{i}) 
         hold on
         legend(lmk,'FontSize',11)
         subplot(212)
         plot(time/3600,v,markers{i})
         hold on
         legend(lmk,'FontSize',11) 
         
         % if i gets to 11 (limit each plot to 10 landmarks)
         if mod(i,10) == 0
            % label plots and break
            subplot(211)
            ylabel('u (pixels)','FontSize',15)
            sgtitle('Simulated Nonlinear Measurements','FontSize',18)
            subplot(212)
            xlabel('Time (hours)','FontSize',15)
            ylabel('v (pixels)','FontSize',15)
            lmk=[];
            break
         end
         
         i = i+1; % increment to next landmark
    end
end
