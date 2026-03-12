% Basic homogeneous sector analysis.
% 1. Stage the homogeneous sector bundle under data/outputs/Sectors/
% 2. Generate stacked sector-allocation plots similar to Figure S1

baseDir = fileparts(mfilename('fullpath'));
repoRoot = fullfile(baseDir, '..', '..');
addpath(fullfile(repoRoot, 'utils'));

inputFolder = fullfile(repoRoot, 'data', 'outputs', 'iML1515_Homogeneous');
stagedFolder = fullfile(repoRoot, 'data', 'outputs', 'Sectors');
outputFolder = fullfile(baseDir, 'output');
ensure_dir(stagedFolder);
ensure_dir(outputFolder);

fluxFile = fullfile(inputFolder, 'dataFlux_Homogeneous.csv');
sectorsFile = fullfile(inputFolder, 'sectors_Homogeneous.csv');
aceFile = fullfile(inputFolder, 'aceT_Homogeneous.csv');

if ~isfile(fluxFile) || ~isfile(sectorsFile)
    error('Homogeneous outputs are missing. Run 1_Glucose_Sweep/run_glucose_sweep_homogeneous.m first.');
end

fprintf('Input folder: %s\n', inputFolder);
fprintf('Staged folder: %s\n', stagedFolder);
fprintf('Output folder: %s\n', outputFolder);

fluxTable = readtable(fluxFile);
sectorsTable = readtable(sectorsFile);

fluxTable = fill_missing_numeric_with_zero(fluxTable);
sectorsTable = fill_missing_numeric_with_zero(sectorsTable);

writetable(fluxTable, fullfile(stagedFolder, 'data_modelV0_kcatscomb_Homog.csv'));
writetable(sectorsTable, fullfile(stagedFolder, 'sectors_modelV0_kcatscomb_Homog.csv'));

if isfile(aceFile)
    aceTable = readtable(aceFile);
    writetable(aceTable, fullfile(stagedFolder, 'aceTData_modelV0_kcatscomb_Homog.csv'));
end

fprintf('Staged sectors bundle written.\n');

gammaValues = unique(sectorsTable.Gamma);
gammaValues = sort(gammaValues(:)');
fprintf('Gamma values to plot: %s\n', mat2str(gammaValues));

color = {sscanf('65AF65','%2x%2x%2x',[1 3])/255, ...
         sscanf('C74848','%2x%2x%2x',[1 3])/255, ...
         sscanf('9065D1','%2x%2x%2x',[1 3])/255, ...
         sscanf('266DD7','%2x%2x%2x',[1 3])/255, ...
         sscanf('F2A1E4','%2x%2x%2x',[1 3])/255, ...
         sscanf('D2C725','%2x%2x%2x',[1 3])/255, ...
         sscanf('B2B2B2','%2x%2x%2x',[1 3])/255, ...
         sscanf('EC8327','%2x%2x%2x',[1 3])/255, ...
         sscanf('35B2B2','%2x%2x%2x',[1 3])/255, ...
         sscanf('E40303','%2x%2x%2x',[1 3])/255};

generatedCount = 0;

for idx = 1:numel(gammaValues)
    gammaVal = gammaValues(idx);
    gammaLabel = ['\gamma = ', num2str(gammaVal)];
    fprintf('\nProcessing gamma = %.3f...\n', gammaVal);

    fluxSlice = fluxTable(fluxTable.Gamma == gammaVal, :);
    if height(fluxSlice) < 2
        fprintf('  Skipping gamma %.3f: not enough flux data\n', gammaVal);
        continue;
    end
    fluxSlice = sortrows(fluxSlice, 'BIOMASS_Ec_iML1515_WT_75p37M');
    grAc = detect_acetate_threshold(fluxSlice);

    sectorsSlice = sectorsTable(sectorsTable.Gamma == gammaVal, :);
    if height(sectorsSlice) < 2
        fprintf('  Skipping gamma %.3f: not enough sectors data\n', gammaVal);
        continue;
    end
    sectorsSlice = sortrows(sectorsSlice, 'Growth');

    phiPmax = max(sectorsSlice.phiPmax);
    phiMmax = max(sectorsSlice.phiMmax);
    [fig, h] = get_plot(sectorsSlice, gammaLabel, color, phiPmax, phiMmax, grAc);

    filenameSvg = fullfile(outputFolder, sprintf('Homogeneous_sectors_g%03d.svg', round(gammaVal * 100)));
    saveas(fig, filenameSvg, 'svg');
    fprintf('  Saved: %s\n', filenameSvg);

    filenamePng = fullfile(outputFolder, sprintf('Homogeneous_sectors_g%03d.png', round(gammaVal * 100)));
    saveas(fig, filenamePng, 'png');
    fprintf('  Saved: %s\n', filenamePng);

    hatchMask = [1 1 1 1 0 1];
    hatchAngles = [-45, 45, -45, 90, 0, 45];
    hold on;
    for hIdx = 1:length(h)
        if hatchMask(hIdx) == 1
            hatchfill2(h(hIdx), 'single', 'HatchAngle', hatchAngles(hIdx), ...
                'HatchDensity', 40, 'HatchColor', 'black');
        end
    end
    hold off;

    filenameHatch = fullfile(outputFolder, sprintf('Homogeneous_sectors_g%03d_hatch.png', round(gammaVal * 100)));
    saveas(fig, filenameHatch, 'png');
    fprintf('  Saved: %s\n', filenameHatch);

    close(fig);
    generatedCount = generatedCount + 1;
end

fprintf('\nSector analysis complete.\n');
fprintf('Generated %d gamma plots in %s\n', generatedCount, outputFolder);


function tableData = fill_missing_numeric_with_zero(tableData)
    variableNames = tableData.Properties.VariableNames;
    for idx = 1:numel(variableNames)
        columnName = variableNames{idx};
        columnData = tableData.(columnName);
        if isnumeric(columnData)
            nanMask = isnan(columnData);
            if any(nanMask)
                columnData(nanMask) = 0;
                tableData.(columnName) = columnData;
            end
        end
    end
end


function grAc = detect_acetate_threshold(fluxTable)
    difV = diff(fluxTable.ATPS4rpp_f);
    difG = diff(fluxTable.BIOMASS_Ec_iML1515_WT_75p37M);
    slopeV = difV ./ difG;
    validSlopes = slopeV(isfinite(slopeV));

    if isempty(validSlopes)
        acIdx = find(fluxTable.EX_ac_e_f >= 0.01, 1, 'first');
        if isempty(acIdx)
            grAc = fluxTable.BIOMASS_Ec_iML1515_WT_75p37M(1);
        else
            grAc = fluxTable.BIOMASS_Ec_iML1515_WT_75p37M(acIdx);
        end
        return;
    end

    maxSlope = max(validSlopes);
    thresholds = [0.01, 0.05, 0.2];
    acetateCutoffs = [0.01, 0, 0];
    acIdx = [];

    for idx = 1:numel(thresholds)
        membraneSat = find(slopeV <= maxSlope * thresholds(idx));
        if isempty(membraneSat)
            continue;
        end
        thresholdIdx = find(fluxTable.EX_ac_e_f(membraneSat) >= acetateCutoffs(idx), 1, 'first');
        if ~isempty(thresholdIdx)
            acIdx = membraneSat(thresholdIdx);
            break;
        end
    end

    if isempty(acIdx)
        acIdx = 1;
    end

    grAc = fluxTable.BIOMASS_Ec_iML1515_WT_75p37M(acIdx);
end


function [fig, h] = get_plot(sectorsTable, gammaLabel, color, phiPmax, phiMmax, grAc)
    fig = figure('Units', 'centimeters', 'Position', [1, 1, 18, 14]);
    axes('Position', [0.1, 0.1, 0.7, 0.8]);
    hold on;

    title(gammaLabel);
    xlabel('Growth rate (h^{-1})', 'FontSize', 12);
    ylabel('Cumulative proteome fraction', 'FontSize', 12);

    gr = sectorsTable.Growth;
    gdata = [sectorsTable.phiCm, sectorsTable.phiEm, sectorsTable.phiEr, ...
             sectorsTable.phiR, sectorsTable.phiEc, sectorsTable.phiCc];
    h = area(gr, gdata);

    set(h(1), 'FaceColor', color{1}, 'LineWidth', 1, 'DisplayName', '\phi C_{m}');
    set(h(2), 'FaceColor', color{6}, 'LineWidth', 1, 'DisplayName', '\phi E_{m}');
    set(h(3), 'FaceColor', color{4}, 'LineWidth', 1, 'DisplayName', '\phi E_{r}');
    set(h(4), 'FaceColor', color{3}, 'LineWidth', 1, 'DisplayName', '\phi R');
    set(h(5), 'FaceColor', color{2}, 'LineWidth', 1, 'DisplayName', '\phi E_{c}');
    set(h(6), 'FaceColor', color{8}, 'LineWidth', 1, 'DisplayName', '\phi C_{c}');

    xlim([0, 1]);
    y = ylim;
    hAc = plot([grAc, grAc], [0, y(2)], '--', 'Color', color{10}, 'LineWidth', 2, ...
        'DisplayName', '\lambda_{ac}');
    ylim(y);

    lgd = legend([h(1), h(2), h(3), h(4), h(5), h(6), hAc]);
    set(lgd, 'Box', 'off', 'Location', 'southoutside', 'NumColumns', 4, 'FontSize', 20);

    hLine1 = line([gr(1), gr(end)], [phiPmax, phiPmax], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 2.5);
    hLine2 = line([gr(1), gr(end)], [phiMmax, phiMmax], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 2.5);

    threshold = 0.05;
    if abs(phiPmax - phiMmax) <= threshold
        text(gr(end) * 1.015, phiPmax, ['\phi P_{max} = ', num2str(phiPmax, '%.2f')], ...
            'VerticalAlignment', 'bottom');
        text(gr(end) * 1.015, phiMmax, ['\phi M_{max} = ', num2str(phiMmax, '%.2f')], ...
            'VerticalAlignment', 'top');
    else
        text(gr(end) * 1.015, phiPmax, ['\phi P_{max} = ', num2str(phiPmax, '%.2f')], ...
            'VerticalAlignment', 'middle');
        text(gr(end) * 1.015, phiMmax, ['\phi M_{max} = ', num2str(phiMmax, '%.2f')], ...
            'VerticalAlignment', 'middle');
    end

    axis([0 gr(end) 0 phiPmax]);
    set(gca, 'FontSize', 12);
    grid on;

    set(get(get(hLine1, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
    set(get(get(hLine2, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
end
