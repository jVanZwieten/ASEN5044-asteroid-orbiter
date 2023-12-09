classdef DTsys
    methods(Static)

        function [F,G] = dynMat(A,B,delT)            
                 
            Ahat = [A B;
                    zeros(length(B(1,:)),length(A(1,:))),zeros(length(B(1,:)),length(B(1,:)))];
            expmA = expm(Ahat*delT);
            
            F = expmA(1:length(A(:,1)),1:length(A(1,:)));
            G = expmA(1:length(A(:,1)),length(A(1,:))+1:length(A(1,:))+length(B(1,:)));
        end

    end
end