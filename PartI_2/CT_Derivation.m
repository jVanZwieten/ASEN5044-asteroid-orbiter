clc
clear
close all

syms x y z xd yd zd r mu real

r = norm([x y z]);

F = [xd yd zd -mu*x/r^3 -mu*y/r^3 -mu*z/r^3]';

pFpX = jacobian(F,[x y z xd yd zd]);

x = 0;
y = -1;
z = 0;
xd = 0;
yd = 0;
zd = sqrt(mu/norm([x y z]));

Abar = subs(pFpX)
Bbar = [zeros(3,3); eye(3,3)];

%%
clear
syms x y z xd yd zd r f l1 l2 l3 ic1 ic2 ic3 jc1 jc2 jc3 kc1 kc2 kc3 real

l = [l1 l2 l3]';
r = [x y z]';
ic = [ic1 ic2 ic3]';
jc = [jc1 jc2 jc3]';
kc = [kc1 kc2 kc3]';


H = [f*((l-r)'*ic)/((l-r)'*kc);
     f*((l-r)'*jc)/((l-r)'*kc)];

pHpX = jacobian(H,[x y z xd yd zd]);

Cbar = subs(pHpX)



