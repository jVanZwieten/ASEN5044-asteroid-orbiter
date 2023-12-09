classdef CTsys
    methods(Static)

        function [A,B] = dynMat(x,y,z,r)            
            global mu

            A = [zeros(3)  eye(3);
                zeros(3) zeros(3)];
            A(4:6,1:3) = [3*mu*x^2/r^5 - mu/r^3    3*mu*x*y/r^5            3*mu*x*z;
                         3*mu*y*x/r^5              3*mu*y^2/r^5 - mu/r^3   3*mu*y*z/r^5;
                         3*mu*z*x/r^5              3*mu*y*z/r^5            3*mu*z^2/r^5 - mu/r^3];

            B = [zeros(3,3); eye(3,3)];
        end


        function [C,D] = measMat(r,l,ic,jc,kc)
            global fo

            C = zeros(2,6);
            C(2,1) = ((fo*kc(1)*(ic'*(l-r)))/((kc'*(l-r))^2))-(fo*ic(1)/(kc'*(l-r)));
            C(2,2) = ((fo*kc(2)*(ic'*(l-r)))/((kc'*(l-r))^2))-(fo*ic(2)/(kc'*(l-r)));
            C(2,3) = ((fo*kc(3)*(ic'*(l-r)))/((kc'*(l-r))^2))-(fo*ic(3)/(kc'*(l-r)));
            C(1,1) = ((fo*kc(1)*(jc'*(l-r)))/((kc'*(l-r))^2))-(fo*jc(1)/(kc'*(l-r)));
            C(1,2) = ((fo*kc(2)*(jc'*(l-r)))/((kc'*(l-r))^2))-(fo*jc(2)/(kc'*(l-r)));
            C(1,3) = ((fo*kc(3)*(jc'*(l-r)))/((kc'*(l-r))^2))-(fo*jc(3)/(kc'*(l-r)));

            D = eye(2);
        end

    end
end