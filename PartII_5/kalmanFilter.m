classdef kalmanFilter
    methods(Static)

        function [xM1,PM1] = timeUpd(xP0,u0,PP0,F,G,Q)
            xM1 = F*xP0 + G*u0;
            PM1 = F*PP0*F' + Q;
        end

        function [xP1,PP1] = measUpd(n,xM1,y1,PM1,H1,R)
            K1 = PM1*H1'*pinv(H1*PM1*H1'+R);
            xP1 = xM1 + K1*(y1 - H1*xM1);
            PP1 = (eye(n,n)-K1*H1)*PM1;

            chol(PP1);  % Should crash if we do something dumb
        end

    end
end