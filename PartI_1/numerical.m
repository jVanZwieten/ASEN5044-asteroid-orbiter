classdef numerical
    methods(Static)
        function X = propagate(X0, dT, t_end, gammaW)
            assert(size(X0, 2) == 1)
            if nargin < 4
                gammaW = zeros(6, 3);
            end
        
            X = zeros(size(X0, 1), length(1:dT:t_end) + 1);
            X(:, 1) = X0;
        
            for k=1:t_end/dT+1
                X(:, k+1) = numerical.rk4_state(X(:, k), dT, gammaW);
            end
        end

        function state_step = rk4_state(state, step, gammaW)
            if nargin < 3
                gammaW = zeros(6, 3);
            end

            k1 = f.twoBody(state,gammaW)+f.solarRadPress();
            
            statea = state+(.5*k1*step);
            k2 = f.twoBody(statea,gammaW)+f.solarRadPress();
            
            stateb = state+(.5*k2*step);
            k3 = f.twoBody(stateb,gammaW)+f.solarRadPress();
            
            statec = state+(k3*step);
            k4 = f.twoBody(statec,gammaW)+f.solarRadPress();

            state_step = state+(k1+2*k2+2*k3+k4)*step/6;
            
        end
    end
end