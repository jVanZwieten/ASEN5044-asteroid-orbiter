function [vis] = isVisible(u,l,r,kc)
    global umin umax

    vis = false;

    if(u(1) >= 0 && u(1) <= umax(1))
        if(u(2) >= 0 && u(2) <= umax(2))
            if((l-r)'*kc > 0)
                if(l'*kc<0)
                    vis = true;
                end
            end
        end
    end
end