classdef CTsys
    properties(Constant)
        B = [zeros(3,3); eye(3,3)];
        Gamma = [zeros(3) zeros(3); % Gamma should be size nxn
                 zeros(3) eye(3)];
        D = zeros(2);
    end

    methods(Static)

        function [AMat,BMat,GammaMat] = dynMat(x,y,z,r)
            AMat = CTsys.AEvaluated([x; y; z; 0; 0; 0]);
            BMat = CTsys.B;
            GammaMat = CTsys.Gamma;
        end

        function A = AEvaluated(X)
            global mu_A
            x = X(1);
            y = X(2);
            z = X(3);
            r = norm(X(1:3));

            A = [zeros(3)  eye(3);
                zeros(3) zeros(3)];
            A(4:6,1:3) = [3*mu_A*x^2/r^5 - mu_A/r^3    3*mu_A*x*y/r^5            3*mu_A*x*z;
                         3*mu_A*y*x/r^5              3*mu_A*y^2/r^5 - mu_A/r^3   3*mu_A*y*z/r^5;
                         3*mu_A*z*x/r^5              3*mu_A*y*z/r^5            3*mu_A*z^2/r^5 - mu_A/r^3];
        end

        function [C,DMat] = measMat(r,l,i_c,j_c,k_c)
            C = CTsys.CEvaluated(r, l, [i_c j_c k_c]);

            DMat = CTsys.D;
        end

        function C = CEvaluated(X, landmarkPosition, rotation_camera)
            assert(size(X, 2) == 1)
            assert(size(landmarkPosition, 1) == 3 && size(landmarkPosition, 2) == 1)
            assert(size(rotation_camera, 1) == 3 && size(rotation_camera, 2) == 3)

            global f_camera
            r = X(1:3);
            i_c = rotation_camera(:, 1);
            j_c = rotation_camera(:, 2);
            k_c = rotation_camera(:, 3);

            C = zeros(2,6);
            C(1,1) = ((f_camera*k_c(1)*(i_c'*(landmarkPosition-r)))/((k_c'*(landmarkPosition-r))^2))-(f_camera*i_c(1)/(k_c'*(landmarkPosition-r)));
            C(1,2) = ((f_camera*k_c(2)*(i_c'*(landmarkPosition-r)))/((k_c'*(landmarkPosition-r))^2))-(f_camera*i_c(2)/(k_c'*(landmarkPosition-r)));
            C(1,3) = ((f_camera*k_c(3)*(i_c'*(landmarkPosition-r)))/((k_c'*(landmarkPosition-r))^2))-(f_camera*i_c(3)/(k_c'*(landmarkPosition-r)));
            C(2,1) = ((f_camera*k_c(1)*(j_c'*(landmarkPosition-r)))/((k_c'*(landmarkPosition-r))^2))-(f_camera*j_c(1)/(k_c'*(landmarkPosition-r)));
            C(2,2) = ((f_camera*k_c(2)*(j_c'*(landmarkPosition-r)))/((k_c'*(landmarkPosition-r))^2))-(f_camera*j_c(2)/(k_c'*(landmarkPosition-r)));
            C(2,3) = ((f_camera*k_c(3)*(j_c'*(landmarkPosition-r)))/((k_c'*(landmarkPosition-r))^2))-(f_camera*j_c(3)/(k_c'*(landmarkPosition-r)));
        end
    end
end