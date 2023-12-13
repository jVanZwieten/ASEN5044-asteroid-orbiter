classdef f
    methods(Static)

        function a2B = twoBody(state)
            global mu_A
            
            r = state(1:3);
            rdot = state(4:6);
            a2B = zeros(6,1);

            a2B(1:3) = rdot;
            a2B(4:6) = -(mu_A/(norm(r))^3)*r;

         end

        function aSRP = solarRadPress()
            global r_sa phi_0 rho Am
            rsaHat = r_sa/norm(r_sa);

            aSRP = zeros(6,1);

            aSRP(4:6) = -(phi_0/norm(r_sa)^2)*(1+(4*rho/9))*Am*rsaHat;
        end

    end
end