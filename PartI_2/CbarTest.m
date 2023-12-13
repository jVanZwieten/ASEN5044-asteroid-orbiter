clc
clear
close all
format longg

addpath(genpath(fileparts(pwd)))

data = load("orbitdetermination-finalproj_data_2023_11_14.mat");
setGlobalVariables()
global mu_A f_camera w_A delT_integration delT_observation t_end
u0 = [512 512]';            % center of optical plane
umin = [0 0]';              % Pixels
umax = [1024 1024]';        % Pixels

tVec = 0:delT_observation:t_end;

l = data.pos_lmks_A(:,1);
r0 = [0 -1 0]';
rdot0 = [0 0 sqrt(mu_A/norm(r0))]';
state0 = [r0; rdot0];

NL_state = zeros(6,length(1:delT_integration:t_end)+2);
NL_state(:,1) = state0;

dx0 = [1e-5 1e-5 1e-5 1e-7 1e-7 1e-7]';
dx = zeros(6,length(1:delT_observation:t_end)+2);
dx(:,1) = dx0;
total_state = zeros(6,length(1:delT_observation:t_end)+2);
total_state(:,1) = state0 + dx0;

y = nan(2,length(1:delT_observation:t_end)+1);
dy = nan(2,length(1:delT_observation:t_end)+1);
yAct = nan(2,length(1:delT_observation:t_end)+1);
NL_y = nan(2,length(1:delT_observation:t_end)+1);

j = 1;
for i=1:(t_end/delT_integration)+1
    NL_state(:,i+1) = numerical.rk4_state(NL_state(:,i),delT_integration);
    

    time = (i-1)*delT_integration;

    if(~mod(time,delT_observation))
        nlx = NL_state(1,i);
        nly = NL_state(2,i);
        nlz = NL_state(3,i);
        r = NL_state(1:3,i);

        [A,B] = CTsys.dynMat(nlx,nly,nlz,norm(r));
        [F,G] = DTsys.dynMat(A,B,CTsys.Gamma,delT_observation);
       
        dx(:,j+1) = F*dx(:,j);
        total_state(:,j) = NL_state(:,j) + dx(:,j);
        
        Rcn = data.R_CtoN(:,:,j);
        ic = Rcn(:,1);
        jc = Rcn(:,2);
        kc = Rcn(:,3);
        
        theta = w_A*time;
        Rna = [cos(theta) -sin(theta) 0;
               sin(theta) cos(theta) 0;
               0 0 1];
        lrot = Rna*l;

        [H,M] = CTsys.measMat(r,lrot,ic,jc,kc);

        NL_y(:,j) = [((f_camera*(lrot-r)'*ic)/((lrot-r)'*kc)) + u0(1);
                     ((f_camera*(lrot-r)'*jc)/((lrot-r)'*kc)) + u0(2)];

        dy(:,j) = H*dx(:,j);
        y(:,j) = NL_y(:,j) + dy(:,j);

        if(~isVisible(y(:,j),lrot,r,kc))
            y(:,j) = nan(2,1);
            dy(:,j) = nan(2,1);
        end

        if(~isVisible(NL_y(:,j),lrot,r,kc))
            NL_y(:,j) = nan(2,1);
        end

        lmks = data.y_table(find(data.y_table(:,1)==time),2:4);    
        lmk1 = lmks(find(lmks(:,1)==1),2:3);

        if(~isempty(lmk1))
            yAct(:,j) = lmk1;
        end

        j = j+1;
    end
end

figure()
plot(tVec/3600, dy(1,:),'x')
title('Perturbation measurements')
xlabel('Time (hours)')
ylabel('\Deltau (pixels)')
hold on


figure()
plot(tVec/3600, dy(2,:),'x')
title('Perturbation measurements')
xlabel('Time (hours)')
ylabel('\Deltav (pixels)')
hold on


figure()
plot(tVec/3600,NL_y(1,:),'o')
xlabel('Time (hours)')
ylabel('u (pixels)')
hold on
plot(tVec/3600,y(1,:),'x')
legend('NL Meas','Linear Meas')


