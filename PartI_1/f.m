classdef f
    methods(Static)

        function a2B = twoBody(state)
            global mu
            
            r = state(1:3);
            rdot = state(4:6);
            a2B = zeros(6,1);

            a2B(1:3) = rdot;
            a2B(4:6) = -(mu/(norm(r))^3)*r;

         end

        function aSRP = solarRadPress()
            global rsa phi0 rho Am
            rsaHat = rsa/norm(rsa);

            aSRP = zeros(6,1);

            aSRP(4:6) = -(phi0/norm(rsa)^2)*(1+(4*rho/9))*Am*rsaHat;
        end

    end
end