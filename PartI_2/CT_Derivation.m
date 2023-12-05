clc
clear
close all

syms x y z xd yd zd mu
assume(mu,"real")

r = sqrt(x^2 + y^2 + z^2);

F = [xd yd zd -mu*x/r^3 -mu*y/r^3 -mu*z/r^3]';

pFpX = jacobian(F,[x y z xd yd zd]);

x = 0;
y = -1;
z = 0;
xd = 0;
yd = 0;
zd = sqrt(mu/norm([x y z]));

Abar = subs(pFpX)
Bbar = [zeros(3,3); eye(3,3)]


