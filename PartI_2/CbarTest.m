clc
clear
close all
format longg

addpath(genpath(fileparts(pwd)))

data = load("orbitdetermination-finalproj_data_2023_11_14.mat");

global mu fo wA u0 umin umax
mu = 4.892E-9;              % gravity parameter of asteroid
fo = 2089.7959;             % Pixels
wA = 2*pi/(4.296057*3600);  % rad/s; asteroid rotation rate
u0 = [512 512]';            % center of optical plane
umin = [0 0]';              % Pixels
umax = [1024 1024]';        % Pixels

delTint = 60;               % s
delTobs = 600;              % s
tEnd = 72*60*60;            % 72h -> s
tVec = 0:delTobs:tEnd;

l = data.pos_lmks_A(:,1);
r0 = [0 -1 0]';
rdot0 = [0 0 sqrt(mu/norm(r0))]';
state0 = [r0; rdot0];

NL_state = zeros(6,length(1:delTint:tEnd)+2);
NL_state(:,1) = state0;

dx0 = [1e-5 1e-5 1e-5 1e-7 1e-7 1e-7]';
dx = zeros(6,length(1:delTobs:tEnd)+2);
dx(:,1) = dx0;
total_state = zeros(6,length(1:delTobs:tEnd)+2);
total_state(:,1) = state0 + dx0;

y = nan(2,length(1:delTobs:tEnd)+1);
dy = nan(2,length(1:delTobs:tEnd)+1);
yAct = nan(2,length(1:delTobs:tEnd)+1);
NL_y = nan(2,length(1:delTobs:tEnd)+1);

j = 1;
for i=1:(tEnd/delTint)+1
    NL_state(:,i+1) = rk4_state(NL_state(:,i),delTint);
    

    time = (i-1)*delTint;

    if(~mod(time,delTobs))
        nlx = NL_state(1,i);
        nly = NL_state(2,i);
        nlz = NL_state(3,i);
        r = NL_state(1:3,i);

        [A,B] = CTsys.dynMat(nlx,nly,nlz,norm(r));
        [F,G] = DTsys.dynMat(A,B,delTobs);
       
        dx(:,j+1) = F*dx(:,j);
        total_state(:,j) = NL_state(:,j) + dx(:,j);
        
        Rcn = data.R_CtoN(:,:,j);
        ic = Rcn(:,1);
        jc = Rcn(:,2);
        kc = Rcn(:,3);
        
        theta = wA*time;
        Rna = [cos(theta) -sin(theta) 0;
               sin(theta) cos(theta) 0;
               0 0 1];
        lrot = Rna*l;

        [H,M] = CTsys.measMat(r,lrot,ic,jc,kc);

        NL_y(:,j) = [((fo*(lrot-r)'*ic)/((lrot-r)'*kc)) + u0(1);
                     ((fo*(lrot-r)'*jc)/((lrot-r)'*kc)) + u0(2)];

        dy(:,j) = H*dx(:,j);
        y(:,j) = NL_y(:,j) + dy(:,j);

        if(~isVisible(y(:,j),lrot,r,kc))
            y(:,j) = nan(2,1);
            dy(:,j) = nan(2,1);
        end

        if(~isVisible(NL_y(:,j),lrot,r,kc))
            NL_y(:,j) = nan(2,1);
            % y(:,j) = nan(2,1);
            % dy(:,j) = nan(2,1);
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

% plot(y(2,:)-yAct(2,:))
% 
% figure()
% plot(y(1,:)-NL_y(1,:))
% hold on
% plot(y(2,:)-NL_y(2,:))
% 
% figure()
% plot(NL_y(1,:)-yAct(1,:))
% hold on
% plot(NL_y(2,:)-yAct(2,:))


figure()
plot(tVec/3600,NL_y(1,:),'o')
xlabel('Time (hours)')
ylabel('u (pixels)')
hold on
plot(tVec/3600,y(1,:),'x')
legend('NL Meas','Linear Meas')

% figure()
% plot(NL_y(1,:),NL_y(2,:),'red o')
% 
% figure()
% plot(yAct(1,:),yAct(2,:),'green o')

