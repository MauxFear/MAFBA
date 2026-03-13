% Supplementary Figure S4
% Transport-reaction acetate threshold sensitivity.

enable_gamma_line_plot = true;
enable_contour_plot = true;

scriptDir = fileparts(mfilename('fullpath'));
repoRoot = fullfile(scriptDir, '..', '..');
dataFile = fullfile(repoRoot, 'data', 'outputs', ...
    'Acetate_Threshold_Transport_Sensitivity', ...
    'acetate_threshold_transport_sensitivity.csv');
outputDir = fullfile(scriptDir, 'output');

if ~isfile(dataFile)
    error('S4 data file not found: %s', dataFile);
end

if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

aceLData = readtable(dataFile);
paramCol = 'Gamma';
cfCol = 'Change_Factor';
paramName = '\gamma';

paramValues = unique(aceLData.(paramCol));
changeFactors = unique(aceLData.(cfCol));
heatmapData = create_heatmap_data(aceLData, paramValues, changeFactors);

if enable_gamma_line_plot
    create_gamma_line_plot_visualization(aceLData, paramCol, paramName, cfCol, outputDir);
end

if enable_contour_plot
    create_contour_plot_visualization(heatmapData, paramValues, changeFactors, paramName, outputDir);
end

fprintf('Saved Supplementary Figure S4 outputs to: %s\n', outputDir);

function heatmapData = create_heatmap_data(aceLData, paramValues, changeFactors)
    heatmapData = NaN(length(paramValues), length(changeFactors));
    for i = 1:length(paramValues)
        for j = 1:length(changeFactors)
            idx = aceLData.Gamma == paramValues(i) & aceLData.Change_Factor == changeFactors(j);
            if any(idx)
                heatmapData(i, j) = aceLData.Consensus_Threshold(find(idx, 1, 'first'));
            end
        end
    end
end

function create_gamma_line_plot_visualization(aceLData, paramCol, paramName, cfCol, outputDir)
    changeFactors = unique(aceLData.(cfCol));
    selectedChangeFactors = [];
    stdValues = [1, 2, 5, 10, 0.5, 0.2, 0.1];
    for val = stdValues
        if any(abs(changeFactors - val) < 0.01)
            selectedChangeFactors(end + 1) = val; %#ok<AGROW>
        end
    end
    if isempty(selectedChangeFactors)
        selectedChangeFactors = changeFactors';
    else
        selectedChangeFactors = sort(selectedChangeFactors);
    end

    thresholdTypes = {'Acetate_Threshold_Main'};
    thresholdLabels = {'Main Threshold (flux > 0.1)'};
    color = {sscanf('65AF65','%2x%2x%2x',[1 3])/255, ...
             sscanf('C74848','%2x%2x%2x',[1 3])/255, ...
             sscanf('9065D1','%2x%2x%2x',[1 3])/255, ...
             sscanf('266DD7','%2x%2x%2x',[1 3])/255, ...
             sscanf('F2A1E4','%2x%2x%2x',[1 3])/255, ...
             sscanf('D2C725','%2x%2x%2x',[1 3])/255, ...
             sscanf('B2B2B2','%2x%2x%2x',[1 3])/255, ...
             sscanf('E40303','%2x%2x%2x',[1 3])/255, ...
             [0.07450980392156863, 0.06666666666666667, 0.0196078431372549], ...
             sscanf('FFFFFF','%2x%2x%2x',[1 3])/255};
    markerShapes = {'o', 'd', 's', '^', 'd', 'v', '>', '<'};

    for t = 1:length(thresholdTypes)
        thresholdType = thresholdTypes{t};
        if ~ismember(thresholdType, aceLData.Properties.VariableNames)
            continue;
        end

        fig = figure('Units', 'centimeters', 'Position', [1, 1, 22, 15]);
        title(['Transport Reactions: Acetate Threshold (' thresholdLabels{t} ') vs. ' paramName], ...
            'FontSize', 16);
        xlabel([paramName, ' (\phi^{M}_{max} / \phi_{max})'], 'FontSize', 14);
        ylabel('Acetate Threshold (h^{-1})', 'FontSize', 14);
        set(gca, 'FontSize', 14);
        hold on;

        lineData = cell(length(selectedChangeFactors), 3);
        legendHandles = [];

        for i = 1:length(selectedChangeFactors)
            cf = selectedChangeFactors(i);
            subset = aceLData(aceLData.(cfCol) == cf, :);
            if height(subset) < 2
                continue;
            end

            [sortedParam, idx] = sort(subset.(paramCol));
            sortedThreshold = subset.(thresholdType)(idx);
            lineData{i, 1} = sortedParam;
            lineData{i, 2} = sortedThreshold;
            lineData{i, 3} = ['CF = ' num2str(cf)];
            plot(sortedParam, sortedThreshold, '-', 'Color', 'k', 'LineWidth', 4, ...
                'HandleVisibility', 'off');
        end

        for i = 1:length(selectedChangeFactors)
            if isempty(lineData{i, 1})
                continue;
            end
            colorIdx = mod(i - 1, 8) + 1;
            plot(lineData{i, 1}, lineData{i, 2}, '-', 'Color', color{colorIdx}, ...
                'LineWidth', 4, 'HandleVisibility', 'off');
        end

        y = ylim;
        if y(1) >= 0
            ylim([0, y(2)]);
        end

        for i = 1:length(selectedChangeFactors)
            if isempty(lineData{i, 1})
                continue;
            end
            sortedParam = lineData{i, 1};
            sortedThreshold = lineData{i, 2};
            if length(sortedParam) <= 10
                indicesToShow = (1:length(sortedParam))';
            else
                lowIdx = find(sortedParam <= 0.3);
                highIdx = find(sortedParam > 0.3);
                lowShow = downsample_indices(lowIdx, 8);
                highShow = downsample_indices(highIdx, 8);
                indicesToShow = sort([lowShow; highShow]);
            end

            colorIdx = mod(i - 1, 8) + 1;
            markerIdx = mod(i - 1, length(markerShapes)) + 1;
            if colorIdx <= 4
                h = plot(sortedParam(indicesToShow), sortedThreshold(indicesToShow), ...
                    markerShapes{markerIdx}, 'Color', color{colorIdx}, 'LineWidth', 1.5, ...
                    'MarkerEdgeColor', color{9}, 'MarkerFaceColor', color{colorIdx}, ...
                    'MarkerSize', 10, 'DisplayName', lineData{i, 3});
            else
                h = plot(sortedParam(indicesToShow), sortedThreshold(indicesToShow), ...
                    markerShapes{markerIdx}, 'Color', color{colorIdx}, 'LineWidth', 1.5, ...
                    'MarkerEdgeColor', color{colorIdx}, 'MarkerFaceColor', color{10}, ...
                    'MarkerSize', 10, 'DisplayName', lineData{i, 3});
            end
            legendHandles = [legendHandles; h]; %#ok<AGROW>
        end

        grid on;
        lgd = legend(legendHandles, 'Location', 'southoutside', ...
            'NumColumns', min(9, length(legendHandles)));
        set(lgd, 'Box', 'off', 'FontSize', 14);
        title(lgd, 'Change Factors', 'FontSize', 14);

        filenameBase = sprintf('acetate_threshold_%s_vs_param', ...
            strrep(lower(thresholdType), 'acetate_threshold_', ''));
        saveas(fig, fullfile(outputDir, [filenameBase '.png']), 'png');
        saveas(fig, fullfile(outputDir, [filenameBase '.svg']), 'svg');
        close(fig);
    end
end

function create_contour_plot_visualization(heatmapData, paramValues, changeFactors, paramName, outputDir)
    fig = figure('Units', 'centimeters', 'Position', [1, 1, 25, 20]);
    [X, Y] = meshgrid(changeFactors, paramValues);
    [C, h] = contourf(X, Y, heatmapData, 20);
    colormap(jet);
    c = colorbar;
    ylabel(c, 'Acetate Threshold (h^{-1})', 'FontSize', 12);
    clabel(C, h, 'FontSize', 8);
    set(gca, 'XScale', 'log');
    xlabel('Change Factor (log scale)', 'FontSize', 14);
    ylabel([paramName, ' (Respiratory Cost Factor)'], 'FontSize', 14);
    title('Transport Reactions: Contour Plot of Acetate Threshold', 'FontSize', 16);
    set(gca, 'FontSize', 10);
    saveas(fig, fullfile(outputDir, 'acetate_threshold_contour_plot.png'), 'png');
    saveas(fig, fullfile(outputDir, 'acetate_threshold_contour_plot.svg'), 'svg');
    close(fig);
end

function indicesToShow = downsample_indices(indices, targetCount)
    if isempty(indices)
        indicesToShow = [];
        return;
    end
    step = max(1, floor(length(indices) / targetCount));
    indicesToShow = indices(1:step:end);
end
