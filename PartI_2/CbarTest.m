clc
clear
close all
format longg

addpath(genpath(fileparts(pwd)))

data = load('orbitdetermination-finalproj_data_2023_11_14.mat');

global mu fo wA
mu = 4.892E-9;
fo = 2089.7959;             % Pixels
wA = 2*pi/(4.296057*3600);

delTint = 60;               % s
delTobs = 600;              % s
tEnd = 72*60*60;            % 72h -> s

l = data.pos_lmks_A(:,1);
% l = [.25 0 0]';             % For Debugging
r0 = [0 -1 0]';
rdot0 = [0 0 sqrt(mu/norm(r0))]';
state0 = [r0; rdot0];

NL_state = zeros(6,length(1:delTint:tEnd)+2);
NL_state(:,1) = state0;

u0 = [512 512]';

y = zeros(2,length(1:delTobs:tEnd)+1);
yAct = nan(2,length(1:delTobs:tEnd)+1);
j = 1;

for i=1:(tEnd/delTint)+1
    NL_state(:,i+1) = rk4_state(NL_state(:,i),delTint);

    time = (i-1)*delTint;

    if(~mod(time,delTobs))
        Rcn = data.R_CtoN(:,:,j);
        ic = Rcn(:,1);
        jc = Rcn(:,2);
        kc = Rcn(:,3);
        
        theta = wA*time;
        Rna = [cos(theta) -sin(theta) 0;
               sin(theta) cos(theta) 0;
               0 0 1];
        lrot = Rna'*l;

        [H,M] = CTsys.measMat(NL_state(1:3,i),lrot,ic,jc,kc);
        y(:,j) = H*NL_state(:,i)+M*u0;

        lmks = data.y_table(find(data.y_table(:,1)==time),2:4);    
        lmk1 = lmks(find(lmks(:,1)==1),2:3);

        if(~isempty(lmk1))
            yAct(:,j) = lmk1;
        end

        j = j+1;
    end
end

% y
% yAct

figure()
plot(y(1,:)-yAct(1,:))
hold on
plot(y(2,:)-yAct(2,:))

figure()
plot(y(1,:),y(2,:),'blue o')

figure()
plot(yAct(1,:),yAct(2,:),'red o')
