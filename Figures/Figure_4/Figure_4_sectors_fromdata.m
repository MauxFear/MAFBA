%% FIGURE 4 - PARAMETERIZED PROTEOME SECTOR ALLOCATION
% Uses the same baseline sensitivity dataset as Figure 5, but plots
% proteome sectors instead of reaction fluxes.

%% CONFIGURATION
baseDir = fileparts(mfilename('fullpath'));
repoRoot = fullfile(baseDir, '..', '..');
addpath(fullfile(repoRoot, 'utils'));

datafolder = fullfile(repoRoot, 'data', 'outputs', 'Parametrized');
fprintf('Data folder: %s\n', datafolder);

output_folder = fullfile(baseDir, 'output');
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end
fprintf('Output folder: %s\n', output_folder);

figurePrefix = 'iML1515_MAFBA_parametrized';
sectorsfilename = fullfile(datafolder, 'sectors_modelV0_kcatscomb_Parametrized.csv');
datafilename = fullfile(datafolder, 'data_modelV0_kcatscomb_Parametrized.csv');

%% LOAD DATA
fprintf('\nLoading data...\n');
tableSectors = readtable(sectorsfilename);
tableData = readtable(datafilename);
fprintf('Sectors data: %d rows\n', height(tableSectors));
fprintf('Flux data: %d rows\n', height(tableData));

gamma_values = [1.0, 0.3, 0.25];
change_factor = 1.0;
fprintf('Gamma values to plot: %s\n', mat2str(gamma_values));
if ismember('Change_Factor', tableData.Properties.VariableNames)
    fprintf('Using Change_Factor = %.1f\n', change_factor);
else
    fprintf('No Change_Factor column found. Using the full dataset.\n');
end

growth_col = resolve_column(tableData.Properties.VariableNames, {'rx_BIOMASS_Ec_iML1515_WT_75p37M', 'BIOMASS_Ec_iML1515_WT_75p37M'});
atp_col = resolve_column(tableData.Properties.VariableNames, {'rx_ATPS4rpp_f', 'ATPS4rpp_f'});
acetate_col = resolve_column(tableData.Properties.VariableNames, {'rx_EX_ac_e_f', 'EX_ac_e_f'});

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

%% GENERATE FIGURES
fprintf('\nGenerating figures for %d gamma values...\n', length(gamma_values));
generated_count = 0;

for f = 1:length(gamma_values)
    gamma_val = gamma_values(f);
    gamma_string = ['\gamma = ', num2str(gamma_val)];
    fprintf('\n  Processing gamma = %.2f...\n', gamma_val);
    
    if ismember('Change_Factor', tableData.Properties.VariableNames)
        T = tableData(abs(tableData.Gamma - gamma_val) < 1e-9 & abs(tableData.Change_Factor - change_factor) < 1e-9, :);
    else
        T = tableData(abs(tableData.Gamma - gamma_val) < 1e-9, :);
    end
    if height(T) < 2
        fprintf('    Skipping gamma %.2f: not enough baseline flux data\n', gamma_val);
        continue;
    end
    T = sortrows(T, growth_col);
    gr_ac = detect_acetate_threshold(T, growth_col, atp_col, acetate_col);
    
    if ismember('Change_Factor', tableSectors.Properties.VariableNames)
        S = tableSectors(abs(tableSectors.Gamma - gamma_val) < 1e-9 & abs(tableSectors.Change_Factor - change_factor) < 1e-9, :);
    else
        S = tableSectors(abs(tableSectors.Gamma - gamma_val) < 1e-9, :);
    end
    if height(S) < 2
        fprintf('    Skipping gamma %.2f: not enough baseline sector data\n', gamma_val);
        continue;
    end
    S = sortrows(S, 'Growth');
    
    phiPmax = max(S.phiPmax);
    phiMmax = max(S.phiMmax);
    [fig, h] = get_plot(S, gamma_string, color, phiPmax, phiMmax, gr_ac);

    filename_svg = sprintf('%s/%s_sectors_g%03d_cf%03d.svg', output_folder, figurePrefix, round(gamma_val * 100), round(change_factor * 100));
    saveas(fig, filename_svg, 'svg');
    fprintf('    Saved: %s\n', filename_svg);

    filename_png = sprintf('%s/%s_sectors_g%03d_cf%03d.png', output_folder, figurePrefix, round(gamma_val * 100), round(change_factor * 100));
    saveas(fig, filename_png, 'png');
    fprintf('    Saved: %s\n', filename_png);
    
    hatch_mask = [1 1 1 1 0 1];
    hatch_angles = [-45, 45, -45, 90, 0, 45];
    hold on;
    for i = 1:length(h)
        if hatch_mask(i) == 1
            hatchfill2(h(i), 'single', 'HatchAngle', hatch_angles(i), 'HatchDensity', 40, 'HatchColor', 'black');
        end
    end
    hold off;

    filename_hatch = sprintf('%s/%s_sectors_g%03d_cf%03d_hatch.png', output_folder, figurePrefix, round(gamma_val * 100), round(change_factor * 100));
    saveas(fig, filename_hatch, 'png');
    fprintf('    Saved: %s\n', filename_hatch);
    close(fig);
    generated_count = generated_count + 1;
end

fprintf('\n✓ Figure 4 generation complete!\n');
fprintf('  Total figures: %d gamma values × 3 formats = %d files\n', generated_count, generated_count * 3);


function gr_ac = detect_acetate_threshold(T, growth_col, atp_col, acetate_col)
    dif_v = diff(T.(atp_col));
    dif_g = diff(T.(growth_col));
    slope_v = dif_v ./ dif_g;
    valid_slopes = slope_v(isfinite(slope_v));
    if isempty(valid_slopes)
        ac_idx = find(T.(acetate_col) >= 0.01, 1, 'first');
        if isempty(ac_idx)
            gr_ac = T.(growth_col)(1);
        else
            gr_ac = T.(growth_col)(ac_idx);
        end
        return;
    end

    m_s = max(valid_slopes);
    thresholds = [0.01, 0.05, 0.2];
    acetate_cutoffs = [0.01, 0, 0];
    ac_t = [];

    for i = 1:numel(thresholds)
        memb_s = find(slope_v <= m_s * thresholds(i));
        if isempty(memb_s)
            continue;
        end
        t_i = find(T.(acetate_col)(memb_s) >= acetate_cutoffs(i), 1, 'first');
        if ~isempty(t_i)
            ac_t = memb_s(t_i);
            break;
        end
    end

    if isempty(ac_t)
        ac_t = 1;
    end

    gr_ac = T.(growth_col)(ac_t);
end


function column_name = resolve_column(variable_names, candidates)
    column_name = '';
    for idx = 1:numel(candidates)
        if ismember(candidates{idx}, variable_names)
            column_name = candidates{idx};
            return;
        end
    end
    error('None of the candidate columns were found: %s', strjoin(candidates, ', '));
end


function [fig, h] = get_plot(T, gamma_string, color, phiPmax, phiMmax, gr_ac)
    fig = figure('Units', 'centimeters', 'Position', [1, 1, 18, 14]);
    axes('Position', [0.1, 0.1, 0.7, 0.8]);
    hold on;

    title(gamma_string);
    xlabel('Growth rate (h^{-1})', 'FontSize', 12);
    ylabel('Cumulative proteome fraction', 'FontSize', 12);

    gr = T.Growth;
    gdata = [T.phiCm, T.phiEm, T.phiEr, T.phiR, T.phiEc, T.phiCc];
    h = area(gr, gdata);

    set(h(1), 'FaceColor', color{1}, 'LineWidth', 1, 'DisplayName', '\phi C_{m}');
    set(h(2), 'FaceColor', color{6}, 'LineWidth', 1, 'DisplayName', '\phi E_{m}');
    set(h(3), 'FaceColor', color{4}, 'LineWidth', 1, 'DisplayName', '\phi E_{r}');
    set(h(4), 'FaceColor', color{3}, 'LineWidth', 1, 'DisplayName', '\phi R');
    set(h(5), 'FaceColor', color{2}, 'LineWidth', 1, 'DisplayName', '\phi E_{c}');
    set(h(6), 'FaceColor', color{8}, 'LineWidth', 1, 'DisplayName', '\phi C_{c}');

    xlim([0, 1]);
    y = ylim;
    h8 = plot([gr_ac, gr_ac], [0, y(2)], '--', 'Color', color{10}, 'LineWidth', 2, 'DisplayName', '\lambda_{ac}');
    ylim(y);

    lgd = legend([h(1), h(2), h(3), h(4), h(5), h(6), h8]);
    set(lgd, 'Box', 'off', 'Location', 'southoutside', 'NumColumns', 4, 'FontSize', 20);

    hLine1 = line([gr(1), gr(end)], [phiPmax, phiPmax], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 2.5);
    hLine2 = line([gr(1), gr(end)], [phiMmax, phiMmax], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 2.5);

    threshold = 0.05;
    if abs(phiPmax - phiMmax) <= threshold
        text(gr(end) * 1.015, phiPmax, ['\phi P_{max} = ', num2str(phiPmax, '%.2f')], 'VerticalAlignment', 'bottom');
        text(gr(end) * 1.015, phiMmax, ['\phi M_{max} = ', num2str(phiMmax, '%.2f')], 'VerticalAlignment', 'top');
    else
        text(gr(end) * 1.015, phiPmax, ['\phi P_{max} = ', num2str(phiPmax, '%.2f')], 'VerticalAlignment', 'middle');
        text(gr(end) * 1.015, phiMmax, ['\phi M_{max} = ', num2str(phiMmax, '%.2f')], 'VerticalAlignment', 'middle');
    end

    axis([0 gr(end) 0 phiPmax]);
    set(gca, 'FontSize', 12);
    grid on;

    set(get(get(hLine1, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
    set(get(get(hLine2, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
end
