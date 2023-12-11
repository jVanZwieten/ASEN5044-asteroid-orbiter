classdef DTsys
    methods(Static)

        function [F,G,Omega] = dynMat(A,B,Gamma,delT)            
                 
            Ahat = [A B;
                    zeros(length(B(1,:)),length(A(1,:))),zeros(length(B(1,:)),length(B(1,:)))];
            expmA = expm(Ahat*delT);
            
            F = expmA(1:length(A(:,1)),1:length(A(1,:)));
            G = expmA(1:length(A(:,1)),length(A(1,:))+1:length(A(1,:))+length(B(1,:)));
            Omega = delT*Gamma;


        end

        function [Q] = noiseMat(W,delT)     %(A,Gamma,W,delT)
            Q = W*[(delT^3)/3 0 0 (delT^2)/2 0 0;
                    0 (delT^3)/3 0 0 (delT^2)/2 0;
                    0 0 (delT^3)/3 0 0 (delT^2)/2;
                    (delT^2)/2 0 0 delT 0 0;
                    0 (delT^2)/2 0 0 delT 0;
                    0 0 (delT^2)/2 0 0 delT];
            
%             n = length(A);
%             Z = delT*[-A Gamma*W*Gamma';
%                       zeros(n,n) A'];
%             expmZ = expm(Z);
%                         
%             F = expmZ(n+1:end,n+1:end);
%             FinvQ = expmZ(1:n,n+1:end);
%             Q2 = F*FinvQ;
        end

    end
end