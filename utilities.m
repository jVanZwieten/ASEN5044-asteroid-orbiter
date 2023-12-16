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
        

    end
end