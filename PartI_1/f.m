classdef f
    methods(Static)

        function a2B = twoBody(state,gammaW)
            global mu_A sigma_w

            wTilde = noiseMaker(zeros(3,1),(sigma_w^2)*eye(3),1)';
            r = state(1:3);
            rdot = state(4:6);
            a2B = zeros(6,1);

            a2B(1:3) = rdot;
            a2B(4:6) = -(mu_A/(norm(r))^3)*r;
            a2B(1:6) = a2B + gammaW*wTilde;
            
         end

        function aSRP = solarRadPress()
            global r_sa phi_0 rho Am
            rsaHat = r_sa/norm(r_sa);

            aSRP = zeros(6,1);

            aSRP(4:6) = -(phi_0/norm(r_sa)^2)*(1+(4*rho/9))*Am*rsaHat;
        end

    end
end