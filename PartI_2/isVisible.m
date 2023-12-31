function [vis] = isVisible(u,l,r,kc)
    global u_min u_max v_min v_max

    vis = u(1) >= u_min && u(1) <= u_max ...
        && u(2) >= v_min && u(2) <= v_max ...
        && (l-r)'*kc > 0 ...
        && l'*kc<0;
end