clc
clear
close all
format short

global muA rsa phi0 rho Am
muA = 4.892E-9;
rsa = [1.5E8 0 0]';
phi0 = 1E14;
rho = 0.4;
Am = (1/62)*10^(-6);

delTobs = 600;

A = [0 0 0 1 0 0;
     0 0 0 0 1 0;
     0 0 0 0 0 1;
     -muA 0 0 0 0 0;
     0 2*muA 0 0 0 0;
     0 0 -muA 0 0 0];

B = [zeros(3,3); eye(3,3)];

Ahat = [A B;
        zeros(length(B(1,:)),length(A(1,:))),zeros(length(B(1,:)),length(B(1,:)))];
expmA = expm(Ahat*delTobs);

F = expmA(1:length(A(:,1)),1:length(A(1,:)))
G = expmA(1:length(A(:,1)),length(A(1,:))+1:length(A(1,:))+length(B(1,:)))


tEnd = 72*60*60;

r0 = [0 -1 0]';
rdot0 = [0 0 sqrt(muA/norm(r0))]';
x0 = [r0; rdot0];
delx0 = [1e-5 1e-5 1e-5 1e-7 1e-7 1e-7]';

asrp = f.solarRadPress();
u = asrp(4:6);

x = zeros(6,length(1:delTobs:tEnd)+1);
x(:,1) = x0;
perx = zeros(6,length(1:delTobs:tEnd)+1);
perx(:,1) = x0+delx0;

for i=1:tEnd/delTobs
    x(:,i+1) = F*x(:,i);%+G*u;
    perx(:,i+1) = F*perx(:,i);%+G*u;
end

tVec = (1:delTobs:tEnd)/60/60;
delx = perx-x;

val = length(x)-1;

figure()
subplot(6,1,1)
plot(1:val,delx(1,1:val))

subplot(6,1,2)
plot(1:val,delx(2,1:val))

subplot(6,1,3)
plot(1:val,delx(3,1:val))

subplot(6,1,4)
plot(1:val,delx(4,1:val))

subplot(6,1,5)
plot(1:val,delx(5,1:val))

subplot(6,1,6)
plot(1:val,delx(6,1:val))





