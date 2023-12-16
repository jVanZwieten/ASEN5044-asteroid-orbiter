classdef utilities
    methods(Static)
        function stds = standardDeviationsFromCovarienceMatricies(P)
            dimentions=size(P);
            assert(dimentions(1) == dimentions(2))
            stds=zeros(dimentions(1), dimentions(3));
            
            for k=1:dimentions(3)
                stds(:, k)=utilities.standardDeviationsFromCovarienceMatrix(P(:, :, k));
            end
        end
        
        function stds = standardDeviationsFromCovarienceMatrix(P)
            stds = sqrt(diag(P));
        end
        
        function confidencePlots(X, Y, P, confidenceInterval, figureTitle, labels, stateUnits)
            assert(size(X, 2) == size(Y, 2) && size(X, 2) == size(P, 3))
            assert(size(Y, 1) == size(P, 1) && size(P, 1) == size(P, 2))
            
            sigmas = utilities.standardDeviationsFromCovarienceMatricies(P);
            n = size(Y, 1);
            
            figure
            sgtitle(figureTitle)
            for i=1:n
                utilities.setSubplotByColumns(i, 2, n)
                stateElement = Y(i, :);
                stateElementSigmas = sigmas(i, :);
                
                plot(X, stateElement, X, stateElement + stateElementSigmas*confidenceInterval, '--', X, stateElement - stateElementSigmas*confidenceInterval, '--')
                utilities.documentConfidencePlot(i, labels(i), stateUnits(i), confidenceInterval)
            end
        end

        function setSubplotByColumns(index, columns, plots)
            rows = ceil(plots/columns);
            row = mod((index - 1), rows) + 1;
            column = ceil(index/rows);
            subplot(rows, columns, column + (row - 1)*columns)
        end

        function documentConfidencePlot(i, xLabel, xUnit, confidenceInterval)
                title("x_" + i)
                ylabel(xLabel + " (" + xUnit + ")")
                legend("x_" + i + " (" + xUnit + ")", "+" + confidenceInterval + "sigma", "-" + confidenceInterval + "sigma");
        end
        
        function plot3D(X, bounds, plotTitle)
            figure
            view(3)
            hold on
            colors = jet(size(X, 2));
            for timeStep = 1:(size(X, 2) - 1)
                plot3(X(1, timeStep:(timeStep + 1)), X(2, timeStep:(timeStep + 1)), X(3, timeStep:(timeStep + 1)), ...
                    'color', colors(timeStep, :))
            end
            
            hold off
            xlim(bounds(1, :))
            ylim(bounds(2, :))
            zlim(bounds(3, :))
            title(plotTitle)
        end
        
        function errorPlots(T, X_true, X_estimate, P, confidenceInterval, figureTitle, xlabels, xUnits)
            n = size(X_true, 1);
            assert(size(X_estimate, 1) == n && size(P, 1) == n)
            assert(utilities.isSquare(P(:, :, 1)))
            Tsize = size(T, 2);
            assert(size(X_true, 2) == Tsize && size(X_estimate, 2) == Tsize && size(P, 3) == Tsize)
            
            error = X_true - X_estimate ;
            sigmas = utilities.standardDeviationsFromCovarienceMatricies(P);
            
            figure
            sgtitle(figureTitle)
            for i = 1:n
                utilities.setSubplotByColumns(i, 2, n)
                plot(T, error(i, :), T, confidenceInterval*sigmas(i, :), '--', T, -confidenceInterval*sigmas(i, :), '--')
                utilities.documentConfidencePlot(i, xlabels(i), xUnits(i), confidenceInterval)
            end
        end
        
        function square = isSquare(M)
            square = size(M, 1) == size(M, 2);
        end
        
        function NEES = calculateNEES(X_truth, X, P)
            steps = size(X_truth, 2);
            assert(size(X, 2) == steps && size(P, 3) == steps)
            assert(size(X_truth, 1) == size(X, 1) && size(X, 1) == size(P, 1) && size(P, 1) == size(P, 2))

            NEES = zeros(1, size(X_truth, 2));
            error = X_truth - X;
            for k = 1:steps
                NEES(k) = error(:, k)'*P(:, :, k)^-1*error(:, k);
            end
        end

        function NEESPlot(NEES_samples, figureTitle, boundToR1x)
            global n
            k = size(NEES_samples, 1);

            epsNEESbar = mean(NEES_samples,1);
            alphaNEES = 0.05; %%significance level
            Nnx = k*n;
            r1x = chi2inv(alphaNEES/2, Nnx )./ k;
            r2x = chi2inv(1-alphaNEES/2, Nnx )./ k;

            figure
            plot(epsNEESbar,'ro','MarkerSize',6,'LineWidth',2),hold on
            plot(r1x*ones(size(epsNEESbar)),'r--','LineWidth',2)
            plot(r2x*ones(size(epsNEESbar)),'r--','LineWidth',2)
            ylabel('NEES statistic, $\bar{\epsilon}_x$','Interpreter','latex', 'FontSize',14)
            xlabel('time step, k','FontSize',14)
            title('NEES Estimation Results, ' + figureTitle,'FontSize',14)
            legend('NEES @ time k', 'r_1 bound', 'r_2 bound'),grid on

            if nargin > 2 && boundToR1x
                ylim([r1x-2 r2x+2])
            end
        end

        function NISPlot(NIS_samples, figureTitle, boundToR1y)
            global p
            k = size(NIS_samples, 2);

            epsNISbar = mean(NIS_samples,1);
            alphaNIS = 0.05;
            Nny = k*p;
            r1y = chi2inv(alphaNIS/2, Nny )./ k;
            r2y = chi2inv(1-alphaNIS/2, Nny )./ k;

            figure
            plot(epsNISbar,'bo','MarkerSize',6,'LineWidth',2),hold on
            plot(r1y*ones(size(epsNISbar)),'b--','LineWidth',2)
            plot(r2y*ones(size(epsNISbar)),'b--','LineWidth',2)
            ylabel('NIS statistic, $\bar{\epsilon}_y$','Interpreter','latex','FontSize',14)
            xlabel('time step, k','FontSize',14)
            title('NIS Estimation Results, ' + figureTitle,'FontSize',14)
            legend('NIS @ time k', 'r_1 bound', 'r_2 bound'),grid on

            if nargin > 2 && boundToR1y
                ylim([r1y-2 r2y+2])
            end
        end

        function plotSimulatedYNoise(Y_noiseless, Y_noisy)
            figure()
            for ii = 1:3
                meastime = Y_noisy(1, Y_noisy(2, :)==ii);
                noisy_u = Y_noisy(3, Y_noisy(2, :)==ii);
                noisy_v = Y_noisy(4, Y_noisy(2, :)==ii);
                noiseless_u = Y_noiseless(3, Y_noiseless(2, :)==ii);
                noiseless_v = Y_noiseless(4, Y_noiseless(2, :)==ii);
                subplot(211)
                plot(meastime/3600,noisy_u,'x')
                hold on
                plot(meastime/3600,noiseless_u,'o')
                ylabel('u (pixels)')

                subplot(212)
                plot(meastime/3600,noisy_v,'x')
                hold on
                plot(meastime/3600,noiseless_v,'o')
                ylabel('v (pixels)')
            end
            sgtitle('Simulated NL Measurements (lmks 1-3)')
            legend('Noisy','Noiseless')

            xlabel('Time(hours)')
        end

        function  plotSimulatedTrajectories(T, X_noiseless, X_noisy)
            figure()
            subplot(611)
            sgtitle('Typical Noisy Truth Trajectory')
            plot(T, X_noiseless(1,:), T, X_noisy(1,:),'--')
            ylabel('x (km)')
            legend('Noiseless','Noisy')

            subplot(612)
            plot(T,X_noiseless(2,:), T,X_noisy(2,:),'--')
            ylabel('y (km)')

            subplot(613)
            plot(T,X_noiseless(3,:),T,X_noisy(3,:),'--')
            ylabel('z (km)')

            subplot(614)
            plot(T,X_noiseless(4,:),T,X_noisy(4,:),'--')
            ylabel('xdot (km/s)')

            subplot(615)
            plot(T,X_noiseless(5,:),T,X_noisy(5,:),'--')
            ylabel('ydot (km/s)')

            subplot(616)
            plot(T,X_noiseless(6,:),T,X_noisy(6,:),'--')
            ylabel('zdot (km/s)')
            xlabel('Time (hours)')
        end

        function plotStateEstimationErrors (T, X_truth, X_estimate, P)
            error = X_estimate - X_truth;
            sigmas = utilities.standardDeviationsFromCovarienceMatricies(P);

            figure()
            sgtitle('Position state estimation error')
            subplot(311)
            plot(T, error(1,:))
            hold on
            plot(T,2*sigmas(1, :),'black --')
            plot(T,-2*sigmas(1, :),'black --')
            ylabel('x (km)')
            legend('Error','\pm2\sigma bounds')

            subplot(312)
            plot(T,error(2,:))
            hold on
            plot(T,2*sigmas(2, :),'black --')
            plot(T,-2*sigmas(2, :),'black --')
            ylabel('y (km)')

            subplot(313)
            plot(T,error(3,:))
            hold on
            plot(T,2*sigmas(3, :),'black --')
            plot(T,-2*sigmas(3, :),'black --')
            ylabel('z (km)')
            xlabel('Time (hours)')
            %xlim([0 10])
            %ylim([-1e-3 1e-3])

            figure()
            sgtitle('Velocity state estimation error')
            subplot(311)
            plot(T,error(4,:))
            hold on
            plot(T,2*sigmas(4, :),'black --')
            plot(T,-2*sigmas(4, :),'black --')
            ylabel('xdot (km)')
            legend('Error','\pm2\sigma bounds')
            %xlim([0 10])
            %ylim([-1e-4 1e-4])

            subplot(312)
            plot(T,error(5,:))
            hold on
            plot(T,2*sigmas(5, :),'black --')
            plot(T,-2*sigmas(5, :),'black --')
            ylabel('ydot (km)')
            %xlim([0 10])
            %ylim([-1e-3 1e-3])

            subplot(313)
            plot(T,error(6,:))
            hold on
            plot(T,2*sigmas(6, :),'black --')
            plot(T,-2*sigmas(6, :),'black --')
            ylabel('zdot (km)')
            xlabel('Time (hours)')
        end
    end
end