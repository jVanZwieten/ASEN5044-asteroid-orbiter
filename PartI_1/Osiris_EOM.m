function [x_dot] = Osiris_EOM(t,x0)
%   Equations of motion that govern Osiris Rex in its close proximity to
%   the asteroid Bennu. Solar radiation (SRP) is the predominant peturbing 
%   force and is the only one included here.

% Constants:
r_S = [1.5e8 0 0]'; % Position of the sun relative to Bennu's center
phi0 = 1e14; % SRP constant (kg*km/s^2)
rho = 0.4; % Coefficient of reflectivity (-)
Am = (1/62)*1e-6; % Area to Mass ratio of Osiris Rex (km^2/kg)
mu_A = 4.892e-9;

% reshape incoming initial state vector to 6x1 of the form x0=[r;v]
x0 = reshape(x0,[],1);

% extract position and velocity vectors from x0
r = x0(1:3);
v = x0(4:6);

% initialize x_dot
x_dot = zeros(6,1);

% compute x_dot vector (derivative of x0)
x_dot(1:3) = v;
x_dot(4:6) = -(mu_A/(norm(r))^3) * r + ...
    -(phi0/(norm(r_S))^2)*(1 + (4/9)*rho)*(Am)*r_S/norm(r_S);


end