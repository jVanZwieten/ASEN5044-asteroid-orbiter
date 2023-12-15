classdef LinearizedKalmanFilter
    methods(Static)



        function [NL_state,noisy_NL_state] = genNLState(dx0,gammaW)
            setGlobalVariables()
            
            NL_state = zeros(6,length(1:delT_integration:t_end)+2);
            NL_state(:,1) = X_0+dx0;
            noisy_NL_state = NL_state;

            for i=1:(t_end/delT_integration)+1
                wTilde = noiseMaker(zeros(3,1),(sigma_w^2)*eye(3),1)';

                NL_state(:,i+1) = numerical.rk4_state(NL_state(:,i),delT_integration);%+ gammaW*wTilde(:,i);
                noisy_NL_state(:,i+1) = numerical.rk4_state(noisy_NL_state(:,i),delT_integration) + gammaW*wTilde;
            end
        end



        function [y_noiseless,y_table] = genNLMeas(NL_state,data)
            setGlobalVariables()

            R = (sigma_u^2)*eye(2);

            y_table = [];
            y_noiseless=[];

            for j=1:(t_end/delT_observation)+1
                time = (j-1)*delT_observation;

                if(time <= 72*60*60)
                    NLind = ((j-1)*10)+1;

                    r = NL_state(1:3,NLind);

                    Rcn = data.R_CtoN(:,:,j);
                    ic = Rcn(:,1);
                    jc = Rcn(:,2);
                    kc = Rcn(:,3);
            
                    theta = w_A*time;
                    Rna = [cos(theta) -sin(theta) 0;
                           sin(theta) cos(theta) 0;
                           0 0 1];
            
                    for l=1:50       
                        lpos = data.pos_lmks_A(:,l);
                        lrot = Rna*lpos;

                        meas_noise = noiseMaker(zeros(2,1),R,1);

                        u_noiseless = ((f_camera*(lrot-r)'*ic)/((lrot-r)'*kc)) + u_0;
                        v_noiseless = ((f_camera*(lrot-r)'*jc)/((lrot-r)'*kc)) + v_0;

                        u_y = u_noiseless+meas_noise(1);
                        v_y = v_noiseless+meas_noise(1);
                        
                        if(isVisible([u_y; v_y],lrot,r,kc))
                            y_table(end+1,:) = [time l u_y v_y];
                            y_noiseless(end+1,:) = [time l u_noiseless v_noiseless];
                        end
                    end
                end
            end
        end


        function [xP,P_P,filt_total_state,NEES_hist] = LKF(NL_state,dx0,y_table,data,Qkf)
            setGlobalVariables()

            nMeas = length(1:delT_observation:t_end)+1;
            R = (sigma_u^2)*eye(2);

            xM = zeros(n,nMeas);
            P_M = zeros(n,n,nMeas);

            xP = zeros(n,nMeas+1);
            P_P = zeros(n,n,nMeas+1);
            xP(:,1) = dx0;
            P_P(:,:,1) = P_0;

            filt_total_state = zeros(n,nMeas);   
            NEES_hist = zeros(1,nMeas);

            for j=1:(t_end/delT_observation)+1
                time = (j-1)*delT_observation;
                NLind = ((j-1)*10)+1;
            
                x = NL_state(1,NLind);
                y = NL_state(2,NLind);
                z = NL_state(3,NLind);
                NL_r = norm(NL_state(1:3,NLind));
            
                [A,B,Gamma] = CTsys.dynMat(x,y,z,NL_r);
                [F,G,Omega] = DTsys.dynMat(A,B,Gamma,delT_observation);      
            
                [xM(:,j),P_M(:,:,j)] = LinearizedKalmanFilter.timeUpd(xP(:,j),[0 0 0]',P_P(:,:,j),F,G,Qkf,Omega);
                xP_temp = xM(:,j);
                P_P_temp = P_M(:,:,j);
            
                lmks = y_table(find(y_table(:,1)==time),2:4); 
            
                if time <= y_table(end,1) % stop processing measurements after 3days
                    Rcn = data.R_CtoN(:,:,j);
                    ic = Rcn(:,1);
                    jc = Rcn(:,2);
                    kc = Rcn(:,3);
                
                    theta = w_A*time;
                    Rna = [cos(theta) -sin(theta) 0;
                           sin(theta) cos(theta) 0;
                           0 0 1];
                    for l=1:50
                        lmk = lmks(find(lmks(:,1)==l),2:3);
                        
                
                        if(~isempty(lmk))
                            lpos = data.pos_lmks_A(:,l);
                            lrot = Rna*lpos;
                            [H,~] = CTsys.measMat(NL_state(1:3,NLind),lrot,ic,jc,kc);
                            lmk = [0 0] + noiseMaker([0 0],R,1);
                
                            [xP_temp,P_P_temp] = LinearizedKalmanFilter.measUpd(n,xP_temp,lmk',P_P_temp,H,R);
                            %linear_meas(end+1,:) = [time l (H*linear_state(:,j))'];
                        end
                
                    end
                end
                xP(:,j+1) = xP_temp;
                P_P(:,:,j+1) = P_P_temp;
            
                filt_total_state(:,j) = NL_state(:,NLind)+xP(:,j+1);
                NEES_hist(j) = (xP(:,j+1))'*pinv(P_P(:,:,j+1))*(xP(:,j+1));
            end
        end
        


        function [xM1,PM1] = timeUpd(xP0,u0,PP0,F,G,Q,Omega)
            xM1 = F*xP0 + G*u0;
            PM1 = F*PP0*F' + Omega*Q*Omega';
        end



        function [xP1,PP1] = measUpd(n,xM1,y1,PM1,H1,R)
            K1 = PM1*H1'*pinv(H1*PM1*H1'+R);
            xP1 = xM1 + K1*(y1 - H1*xM1);
            PP1 = (eye(n,n)-K1*H1)*PM1;
            % xM1
            H1*xM1;
            PP1;
            chol(PP1);  % Should crash if we do something dumb
        end



    end
end