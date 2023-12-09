clc
clear
close all
format longg

addpath(genpath(fileparts(pwd)))

data = load('orbitdetermination-finalproj_data_2023_11_14.mat');

global mu fo wA u0 umin umax
mu = 4.892E-9;
fo = 2089.7959;             % Pixels
wA = 2*pi/(4.296057*3600);
u0 = [512 512]';
umin = [0 0]';
umax = [1024 1024]';

delTint = 60;               % s
delTobs = 600;              % s
tEnd = 72*60*60;            % 72h -> s

l = data.pos_lmks_A(:,10);
% l = [.25 0 0]';             % For Debugging
r0 = [0 -1 0]';
rdot0 = [0 0 sqrt(mu/norm(r0))]';
state0 = [r0; rdot0];
asrp = f.solarRadPress();
asrp = asrp(4:6);

NL_state = zeros(6,length(1:delTint:tEnd)+2);
NL_state(:,1) = state0;

linear_state = zeros(6,length(1:delTobs:tEnd)+1);
linear_state(:,1) = state0;

meas = nan(2,length(1:delTobs:tEnd)+1);
yAct = nan(2,length(1:delTobs:tEnd)+1);
NL_y = nan(2,length(1:delTobs:tEnd)+1);

j = 1;

for i=1:(tEnd/delTint)+1
    NL_state(:,i+1) = rk4_state(NL_state(:,i),delTint);

    time = (i-1)*delTint;

    if(~mod(time,delTobs))
        x = NL_state(1,i);
        y = NL_state(2,i);
        z = NL_state(3,i);
        NL_r = norm(NL_state(1:3,i));

        [A,B] = CTsys.dynMat(x,y,z,NL_r);
        [F,G] = DTsys.dynMat(A,B,delTobs);

        linear_state(:,j+1) = F*linear_state(:,j)+G*asrp;

        r = linear_state(1:3,j);

        Rcn = data.R_CtoN(:,:,j);
        ic = Rcn(:,1);
        jc = Rcn(:,2);
        kc = Rcn(:,3);
        
        theta = wA*time;
        Rna = [cos(theta) -sin(theta) 0;
               sin(theta) cos(theta) 0;
               0 0 1];
        lrot = Rna*l;
%         lrot = l;             % For Debugging

        [H,M] = CTsys.measMat(NL_r,lrot,ic,jc,kc);
        meas(:,j) = H*linear_state(:,j)+M*u0;
        NL_y(:,j) = [((fo*(lrot-r)'*ic)/((lrot-r)'*kc)) + u0(1);
                     ((fo*(lrot-r)'*jc)/((lrot-r)'*kc)) + u0(2)];

        
%   Keep these commented until H matrix works right for debugging

%         if(~isVisible(meas(:,j),lrot,r,kc))
%             meas(:,j) = nan(2,1);
%         end
%         if(~isVisible(NL_y(:,j),lrot,r,kc))
%             NL_y(:,j) = nan(2,1);
%         end

        lmks = data.y_table(find(data.y_table(:,1)==time),2:4);    
        lmk10 = lmks(find(lmks(:,1)==10),2:3);

        if(~isempty(lmk10))
            yAct(:,j) = lmk10;
        end

        j = j+1;
    end
end

% y
% yAct

% figure()
% plot(meas(1,:)-yAct(1,:))
% hold on
% plot(meas(2,:)-yAct(2,:))
% 
% figure()
% plot(meas(1,:)-NL_y(1,:))
% hold on
% plot(meas(2,:)-NL_y(2,:))
% 
% figure()
% plot(NL_y(1,:)-yAct(1,:))
% hold on
% plot(NL_y(2,:)-yAct(2,:))


% figure()
% plot(meas(1,:),meas(2,:),'blue o')
% 
% figure()
% plot(NL_y(1,:),NL_y(2,:),'red o')
% 
% figure()
% plot(yAct(1,:),yAct(2,:),'green o')
% 

figure()
plot(1:(tEnd/delTobs)+1,meas(1,:))

figure()
plot(1:(tEnd/delTobs)+1,meas(2,:))

% figure()
% plot(1:(tEnd/delTobs)+1,NL_y(1,:))

% figure()
% plot(1:(tEnd/delTobs)+1,NL_y(2,:))
% 
% 
% figure()
% plot(1:(tEnd/delTobs)+1,yAct(2,:))

