clc
clear
close all
format short

mu = 4.892E-9;
delT = 600;

A = [0 0 0 1 0 0;
     0 0 0 0 1 0;
     0 0 0 0 0 1;
     -mu 0 0 0 0 0;
     0 2*mu 0 0 0 0;
     0 0 -mu 0 0 0];

B = [zeros(3,3); eye(3,3)];

Ahat = [A B;
        zeros(length(B(1,:)),length(A(1,:))),zeros(length(B(1,:)),length(B(1,:)))];
expmA = expm(Ahat*delT);

F = expmA(1:length(A(:,1)),1:length(A(1,:)))
G = expmA(1:length(A(:,1)),length(A(1,:))+1:length(A(1,:))+length(B(1,:)))