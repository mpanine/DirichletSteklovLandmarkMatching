function [ plot_axes ] = plot( results, resultField, plotType, varargin )
    %PLOT Plots the values given the plotType.
    p = inputParser;
    addParameter(p,'metricType','percentage_error', @ischar);
    addParameter(p,'name',resultField, @ischar);
    addParameter(p, 'xlim', [0 1], @isnumeric);
    addParameter(p, 'step', 0.005, @isnumeric);
    addParameter(p,'color_map',brewermap(3,'Accent'), @isnumeric);
    addParameter(p, 'LineStyle', '-', @ischar);
    addParameter(p, 'xOffset', 0, @isnumeric);
    addParameter(p, 'appearance', 'none', @ischar);
    addParameter(p, 'parameterList', []);
    addParameter(p, 'parameterName', 'Parameter');
    addParameter(p, 'shapeList', []);
    addParameter(p,'showAll',false, @islogical);
    parse(p,varargin{:});
    inputs = p.Results;
    
    metricType = inputs.metricType;
    name = inputs.name;
    xlimits = inputs.xlim;
    step = inputs.step;
    color_map = inputs.color_map;
    xOffset = inputs.xOffset;
    appearance = inputs.appearance;
    parameterList = inputs.parameterList;
    parameterName = inputs.parameterName;
    showAll = inputs.showAll;
    LineStyle = inputs.LineStyle;
    
    
    plotArgs_mean.LineWidth = 2;
    plotArgs_mean.LineStyle = LineStyle;
    plotArgs_mean.Color = color_map(1,:);
    plotArgs_mean.Marker = 'none';
    
    plotArgs_all.LineWidth = 0.5;
    plotArgs_all.LineStyle = '-.';
    plotArgs_all.Color = color_map(2,:);
    plotArgs_all.Marker = 'none';
    
    plot_axes = struct;
    
    
    switch plotType
        case 'percentage_error'
            thresholds = 0:step:1;

            curves = zeros(size(results, 2), length(thresholds));
            for i=1:size(results, 2)
                curves(i, :) = compute_curve(thresholds,...
                                results{i}.(resultField).values);

                if showAll
                    plot_axes.values = plot_values(thresholds,...
                                curves(i, :), metricType, xlimits,...
                                plotArgs_all, name);hold on;
                end
            end

            mean_values = mean(curves, 1);
            plot_axes.mean = plot_values(thresholds, mean_values,...
                                metricType, xlimits, plotArgs_mean, name);
        case 'parameterWise_error'
            curve_mean = zeros(size(results, 2), 1);
            curve_worst = zeros(size(results, 2), 1);
            curve_variance = zeros(size(results, 2), 1);
            for i=1:size(results, 2)
                curve_mean(i) = results{i}.(resultField).mean;
                curve_worst(i) = results{i}.(resultField).worstCase;
                curve_variance(i) = results{i}.(resultField).standardDeviation;
            end
            
            curve_stdP = curve_mean + curve_variance;
            curve_stdM = curve_mean - curve_variance;
            
            mkSz_mean = 7;
            mkSz_std = 7;
            lnW_mean = 1;
            lnW_std = 1;
            
            xVals = 1:size(results, 2);
            xVals = xVals + xOffset;
            
            set(groot,'defaultAxesTickLabelInterpreter','none');  
            
            switch appearance
                case 'continuous'
                    hold on;plot_axes.mean = plot(xVals, curve_mean,'LineStyle',LineStyle,...
                        'MarkerSize',mkSz_mean,'LineWidth',lnW_mean,...
                        'MarkerEdgeColor','k','MarkerFaceColor',color_map(1,:),...
                        'Color', color_map(1,:));
                otherwise
                    for i=1:size(results, 2)
                        hold on;plot_axes.standardDeviation = plot([i+xOffset i+xOffset],...
                            [curve_stdM(i) curve_stdP(i)],'Color','k',...
                            'LineWidth',2*lnW_std);
                    end;
                    hold on;plot_axes.mean = plot(xVals, curve_mean,'s',...
                        'MarkerSize',mkSz_mean,'LineWidth',lnW_mean,...
                        'MarkerEdgeColor','k','MarkerFaceColor',color_map(1,:));
                    hold on;plot(xVals, curve_stdP,'v','MarkerSize',...
                        mkSz_std,'LineWidth',lnW_std,'MarkerEdgeColor','k',...
                        'MarkerFaceColor',color_map(2,:));
                    hold on;plot(xVals, curve_stdM,'^','MarkerSize',...
                        mkSz_std,'LineWidth',lnW_std,'MarkerEdgeColor','k',...
                        'MarkerFaceColor',color_map(2,:));
            end
            xticks(1:size(results, 2));
            if isempty(parameterList)
                xticklabels(1:size(results, 2));
            else
                if length(parameterList{1}) >= 10
                    xtickangle(90);
                end
                xticklabels(parameterList);
            end
            xlabel(parameterName);
            ylabel(sprintf('%s Error',metricType));
            title(name);
            
        case 'matrix_error'
            if ~isfield(inputs,'shapeList')
                throw(...
                MException('VISUALIZESINGLEBASISEVALUATIONERROR:noShapeList',...
                'No shapeList keyword argument, specifing the list of unique shapes in the result was given to the function.'));
            end
            
            shapeList = inputs.shapeList;
            numShapes = size(shapeList, 2);
            
            if sqrt(size(results, 2)) ~= numShapes
                throw(...
                MException('VISUALIZESINGLEBASISEVALUATIONERROR:incoherentNumberOfShapes',...
                'The size of results does not match size(shapeList, 2)^2.'));
            end
            
            mean_error_matrix = zeros(numShapes);
            variance_error_matrix = zeros(numShapes);
            worstCase_error_matrix = zeros(numShapes);
            for i=1:numShapes
                for j=1:numShapes
                   mean_error_matrix(i, j) = ...
                       results{i+(j-1)*numShapes}.(resultField).mean;
                   variance_error_matrix(i, j) = ...
                       results{i+(j-1)*numShapes}.(resultField).standardDeviation;
                   worstCase_error_matrix(i, j) = ...
                       results{i+(j-1)*numShapes}.(resultField).worstCase;
                end
            end
            
            figure;
            plot_axes.mean = axes;
            figure;
            plot_axes.standardDeviation = axes;
            figure;
            plot_axes.worstCase = axes;
            
            imagesc(plot_axes.mean, mean_error_matrix);
            imagesc(plot_axes.standardDeviation, variance_error_matrix);
            imagesc(plot_axes.worstCase, worstCase_error_matrix);
            
            xticks(plot_axes.mean, 1:numShapes);
            yticks(plot_axes.mean, 1:numShapes);
            xtickangle(plot_axes.mean, 90);
            xticklabels(plot_axes.mean, shapeList);
            yticklabels(plot_axes.mean, shapeList);
            title(plot_axes.mean, ['Mean ', name, ' Error Matrix']);
            colorbar(plot_axes.mean);
            colormap(plot_axes.mean, color_map);
            
            xticks(plot_axes.standardDeviation, 1:numShapes);
            yticks(plot_axes.standardDeviation, 1:numShapes);
            xtickangle(plot_axes.standardDeviation, 90);
            xticklabels(plot_axes.standardDeviation, shapeList);
            yticklabels(plot_axes.standardDeviation, shapeList);
            title(plot_axes.standardDeviation, ['Standard Deviation ', name, ' Error Matrix']);
            colorbar(plot_axes.standardDeviation);
            colormap(plot_axes.standardDeviation, color_map);
            
            xticks(plot_axes.worstCase, 1:numShapes);
            yticks(plot_axes.worstCase, 1:numShapes);
            xtickangle(plot_axes.worstCase, 90);
            xticklabels(plot_axes.worstCase, shapeList);
            yticklabels(plot_axes.worstCase, shapeList);
            title(plot_axes.worstCase, ['Worst Case ', name, ' Error Matrix']);
            colorbar(plot_axes.worstCase);
            colormap(plot_axes.worstCase, color_map);
        case 'partiality_error'
            partiality = zeros(size(results, 2), 1);
            for i=1:size(results.entries, 2)
                Atotal = results.entries{i}.area_Src;
                Apartial = results.entries{i}.area_Tar;
                partiality(i) = 100-100*Apartial/Atotal;
            end
            
            part_thresh = 100:step:0;
            curves = zeros(size(results, 2), length(part_thresh));
            for i=1:size(results.entries, 2)
                curves(i, :) = compute_partial_curve(part_thresh, partiality,...
                                results.entries{i}.(resultField).values);

                if showAll
                    plot_axes.values = plot_partial_values(part_thresh,...
                                curves(i, :), metricType, xlimits,...
                                plotArgs_all, name);hold on;
                end
            end

            mean_values = mean(curves, 1);
            plot_axes.mean = plot_partial_values(part_thresh, mean_values,...
                                metricType, xlimits, plotArgs_mean, name);
    end
end

function [ curve ] = compute_curve(thresholds, values)
    
    curve = zeros(1, length(thresholds));
    for i=1:length(thresholds)
        curve(i) = 100*sum(values <= thresholds(i)) / length(values);
    end
end

function [ curve ] = compute_partial_curve(thresholds, partiality, values)
    
    curve = zeros(1, length(thresholds));
    for i=1:length(thresholds)
        curve(i) = mean(values(partiality <= thresholds(i)));
    end
end

function [ ax ] = plot_partial_values(thresholds, curve, metricType, xlimits, plotArgs, name)
    ax = plot(thresholds, curve, plotArgs);
    xlim(xlimits);
    ylabel(sprintf('%s Error', metricType));
    xlabel('Partiality (%)');
    title(name);
end

function [ ax ] = plot_values(thresholds, curve, metricType, xlimits, plotArgs, name)
    ax = plot(thresholds, curve, plotArgs);
    xlim(xlimits);
    xlabel(sprintf('%s Error', metricType));
    ylabel('% Correspondences');
    title(name);
end
