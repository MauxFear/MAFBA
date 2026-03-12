% Supplementary Figure S3
% Acetate threshold sensitivity to proteome cost scaling (all reactions).

enable_gamma_line_plot = true;
enable_contour_plot = true;

scriptDir = fileparts(mfilename('fullpath'));
repoRoot = fullfile(scriptDir, '..', '..');
dataFile = fullfile(repoRoot, 'data', 'outputs', 'Acetate_Threshold_Sensitivity', ...
    'acetate_threshold_sensitivity_all_reactions.csv');
outputDir = fullfile(scriptDir, 'output');

if ~isfile(dataFile)
    error('S3 data file not found: %s', dataFile);
end

if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

thresholdData = readtable(dataFile);

if ~ismember('Consensus_Threshold', thresholdData.Properties.VariableNames)
    if ismember('Acetate_threshold', thresholdData.Properties.VariableNames)
        thresholdData.Consensus_Threshold = thresholdData.Acetate_threshold;
    else
        error('Expected Consensus_Threshold or Acetate_threshold column in %s', dataFile);
    end
end

gammaValues = unique(thresholdData.Gamma);
changeFactors = unique(thresholdData.Change_Factor);

selectedChangeFactors = [1, 2, 5, 10, 0.5, 0.2, 0.1];
heatmapData = create_heatmap_data(thresholdData, sort(gammaValues), sort(changeFactors));

if enable_gamma_line_plot
    create_gamma_line_plot_visualization(thresholdData, selectedChangeFactors, outputDir);
end

if enable_contour_plot
    create_contour_plot_visualization(heatmapData, sort(gammaValues), sort(changeFactors), outputDir);
end

fprintf('Saved Supplementary Figure S3 outputs to: %s\n', outputDir);

function heatmapData = create_heatmap_data(thresholdData, gammaValues, changeFactors)
    heatmapData = NaN(length(gammaValues), length(changeFactors));
    for i = 1:length(gammaValues)
        for j = 1:length(changeFactors)
            idx = thresholdData.Gamma == gammaValues(i) & ...
                  thresholdData.Change_Factor == changeFactors(j);
            if any(idx)
                heatmapData(i, j) = thresholdData.Consensus_Threshold(find(idx, 1, 'first'));
            end
        end
    end
end

function create_gamma_line_plot_visualization(thresholdData, selectedChangeFactors, outputDir)
    fig = figure('Units', 'centimeters', 'Position', [1, 1, 22, 15]);

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

    title('Acetate Threshold vs. \gamma for Different Change Factors', 'FontSize', 18);
    xlabel('\gamma (\phi^{M}_{max} / \phi^{P}_{max})', 'FontSize', 14);
    ylabel('Acetate Threshold (h^{-1})', 'FontSize', 14);
    set(gca, 'FontSize', 14);
    hold on;

    lineData = cell(length(selectedChangeFactors), 3);
    legendHandles = [];

    for i = 1:length(selectedChangeFactors)
        cf = selectedChangeFactors(i);
        subset = thresholdData(thresholdData.Change_Factor == cf, :);
        if height(subset) < 2
            continue;
        end

        [sortedGamma, idx] = sort(subset.Gamma);
        sortedThreshold = subset.Consensus_Threshold(idx);

        lineData{i, 1} = sortedGamma;
        lineData{i, 2} = sortedThreshold;
        lineData{i, 3} = ['CF = ' num2str(cf)];

        plot(sortedGamma, sortedThreshold, '-', 'Color', 'k', 'LineWidth', 4, ...
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

        sortedGamma = lineData{i, 1};
        sortedThreshold = lineData{i, 2};
        lowIdx = find(sortedGamma <= 0.3);
        highIdx = find(sortedGamma > 0.3);

        lowShow = downsample_indices(lowIdx, 8);
        highShow = downsample_indices(highIdx, 8);
        indicesToShow = sort([lowShow; highShow]);
        if length(indicesToShow) < 5 && ~isempty(sortedGamma)
            indicesToShow = (1:length(sortedGamma))';
        end

        colorIdx = mod(i - 1, 8) + 1;
        markerIdx = mod(i - 1, length(markerShapes)) + 1;

        if colorIdx <= 4
            h = plot(sortedGamma(indicesToShow), sortedThreshold(indicesToShow), ...
                markerShapes{markerIdx}, 'Color', color{colorIdx}, 'LineWidth', 1.5, ...
                'MarkerEdgeColor', color{9}, 'MarkerFaceColor', color{colorIdx}, ...
                'MarkerSize', 10, 'DisplayName', lineData{i, 3});
        else
            h = plot(sortedGamma(indicesToShow), sortedThreshold(indicesToShow), ...
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

    saveas(fig, fullfile(outputDir, 'acetate_threshold_vs_gamma.png'), 'png');
    saveas(fig, fullfile(outputDir, 'acetate_threshold_vs_gamma.svg'), 'svg');
    close(fig);
end

function create_contour_plot_visualization(heatmapData, gammaValues, changeFactors, outputDir)
    fig = figure('Units', 'centimeters', 'Position', [1, 1, 25, 20]);
    [X, Y] = meshgrid(changeFactors, gammaValues);
    [C, h] = contourf(X, Y, heatmapData, 20);
    colormap(jet);
    c = colorbar;
    ylabel(c, 'Acetate Threshold (h^{-1})', 'FontSize', 12);
    clabel(C, h, 'FontSize', 8);
    set(gca, 'XScale', 'log');
    xlabel('Change Factor (log scale)', 'FontSize', 14);
    ylabel('\gamma (Respiratory Cost Factor)', 'FontSize', 14);
    title('Contour Plot of Acetate Threshold', 'FontSize', 16);
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
