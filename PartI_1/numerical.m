classdef numerical
    methods(Static)
        function X = propagate(X0, dT, t_end)
            assert(size(X0, 2) == 1)
        
            X = zeros(size(X0, 1), length(1:dT:t_end) + 1);
            X(:, 1) = X0;
        
            for k=1:t_end/dT
                X(:, k+1) = numerical.rk4_state(X(:, k), dT);
            end
        end

        function state_step = rk4_state(state,step)

            k1 = f.twoBody(state)+f.solarRadPress();
            
            statea = state+(.5*k1*step);
            k2 = f.twoBody(statea)+f.solarRadPress();
            
            stateb = state+(.5*k2*step);
            k3 = f.twoBody(stateb)+f.solarRadPress();
            
            statec = state+(k3*step);
            k4 = f.twoBody(statec)+f.solarRadPress();

            state_step = state+(k1+2*k2+2*k3+k4)*step/6;
            
        end
    end
end