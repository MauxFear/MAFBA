%% FIGURE 7 - PROTEIN EXPRESSION LEVELS: ACETATE EXCRETION PATTERNS
% Shows how cytosolic (Y) and membrane (Z) protein expression levels
% affect acetate excretion patterns across growth rates
%
% Panel A: Cytosolic protein expression (Y levels)
% Panel B: Membrane protein expression (Z levels)
%
% Data source: Combined protein expression sensitivity analysis
% Outputs: Two acetate excretion plots
%% CONFIGURATION
% Set file paths
baseDir = fileparts(mfilename('fullpath'));
repoRoot = fullfile(baseDir, '..', '..');

% Data folder
datafolder = fullfile(repoRoot, 'data', 'outputs', 'Protein_Expression');
fprintf('Data folder: %s\n', datafolder);

% Output folder
output_folder = fullfile(baseDir, 'output');
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end
fprintf('Output folder: %s\n', output_folder);

%% SETUP
% Output naming
figurePrefix = 'iML1515_MAFBA_protExp';

% Analysis parameters
selected_model = 'MAFBA';  % MAFBA parametrized model
gamma_value = 0.25;            % Membrane constraint parameter

% Protein expression levels
y_levels = [0, 0.01, 0.02, 0.04, 0.06, 0.08];  % Cytosolic (Panel A)
z_levels = [0, 0.02, 0.04, 0.10, 0.15, 0.20];  % Membrane (Panel B)

% Color palettes
% Panel A (Y, cytosolic): Blue/green palette
colors_y = {
    sscanf('6B8E23','%2x%2x%2x',[1 3])/255,... % Olive Green
    sscanf('008080','%2x%2x%2x',[1 3])/255,... % Teal
    sscanf('40E0D0','%2x%2x%2x',[1 3])/255,... % Turquoise
    sscanf('87CEEB','%2x%2x%2x',[1 3])/255,... % Sky Blue
    sscanf('6495ED','%2x%2x%2x',[1 3])/255,... % Cornflower Blue
    sscanf('1F56A8','%2x%2x%2x',[1 3])/255     % Blue
};

% Panel B (Z, membrane): Yellow/red palette
colors_z = {
    sscanf('FFD700','%2x%2x%2x',[1 3])/255,... % Sunny Yellow
    sscanf('FFA500','%2x%2x%2x',[1 3])/255,... % Mango Orange
    sscanf('FF6F61','%2x%2x%2x',[1 3])/255,... % Coral
    sscanf('FF3E3E','%2x%2x%2x',[1 3])/255,... % Tomato Red
    sscanf('DC143C','%2x%2x%2x',[1 3])/255,... % Crimson
    sscanf('C00217','%2x%2x%2x',[1 3])/255     % Dark Red
};

% Markers
markers = {'o', 's', '>', '^', 'd', 'v'};

%% LOAD DATA
fprintf('\n=== LOADING PROTEIN EXPRESSION DATA ===\n');

% Use filtered data file (0.1 MB vs 2.3 GB original)
datafilename = fullfile(datafolder, 'data_protein_expression_figure7.csv');

if ~isfile(datafilename)
    error('Data file not found: %s', datafilename);
end

fprintf('Loading filtered data...\n');
data = readtable(datafilename);
fprintf('Loaded %d rows × %d columns\n', height(data), width(data));

% Note: Data is pre-filtered for MAFBA model, gamma=0.25

%% PANEL A - CYTOSOLIC PROTEIN EXPRESSION (Y)
fprintf('\n=== GENERATING PANEL A: CYTOSOLIC PROTEIN (Y) ===\n');

fig_a = figure('Units', 'centimeters', 'Position', [1, 1, 13, 12.5]);
hold on;

for u = 1:length(y_levels)
    y_level = y_levels(u);
    fprintf('  Processing Y level = %.2f (%d/%d)\n', y_level, u, length(y_levels));
    
    % Filter data: gamma, Y level, Z=0
    mask = data.Gamma == gamma_value & ...
           data.Y_Level == y_level & ...
           data.Z_Level == 0;
    
    if sum(mask) < 2
        fprintf('    Skipping: insufficient data points\n');
        continue;
    end
    
    % Extract data
    T_data = data(mask, :);
    T_data = sortrows(T_data, 'rx_BIOMASS_Ec_iML1515_WT_75p37M');
    
    gr = T_data.rx_BIOMASS_Ec_iML1515_WT_75p37M;
    ac = T_data.rx_EX_ac_e_f;
    
    % Detect acetate threshold
    gr_ac = detect_acetate_threshold(T_data);
    
    % Plot line
    plot(gr, ac, '-', 'Color', colors_y{u}, 'LineWidth', 4, 'HandleVisibility', 'off');
    
    % Plot markers (downsampled)
    downsample_factor = 5;
    indices = 1:downsample_factor:length(gr);
    label_text = sprintf('\\phi_Y %d%%', round(y_level * 100));
    plot(gr(indices), ac(indices), markers{u}, ...
        'Color', colors_y{u}, 'LineWidth', 1.5, ...
        'MarkerEdgeColor', colors_y{u}, 'MarkerFaceColor', 'w', ...
        'MarkerSize', 10, 'DisplayName', label_text);
    
    % Plot acetate threshold
    ymax = 45;
    plot([gr_ac, gr_ac], [0, ymax], '--', 'Color', colors_y{u}, ...
        'LineWidth', 2, 'HandleVisibility', 'off');
    
    fprintf('    λ_ac = %.3f h^{-1}\n', gr_ac);
end

% Customize Panel A
title(sprintf('Cytosolic Protein Expression Level\n\\gamma = %.2f', gamma_value), ...
    'FontSize', 16, 'FontWeight', 'bold');
xlabel('Growth rate (h^{-1})', 'FontSize', 14);
ylabel('Acetate flux (mmol gDW^{-1} h^{-1})', 'FontSize', 14);
xlim([0, 1]);
ylim([0, 45]);
yticks(0:5:45);
legend('Location', 'southoutside', 'NumColumns', 2, 'Box', 'off', 'FontSize', 14);
grid on;
set(gca, 'FontSize', 12);

% Save Panel A
saveas(fig_a, fullfile(output_folder, sprintf('%s_panel_A.svg', figurePrefix)), 'svg');
saveas(fig_a, fullfile(output_folder, sprintf('%s_panel_A.png', figurePrefix)), 'png');
fprintf('  Saved: %s_panel_A.svg/.png\n', figurePrefix);
close(fig_a);

%% PANEL B - MEMBRANE PROTEIN EXPRESSION (Z)
fprintf('\n=== GENERATING PANEL B: MEMBRANE PROTEIN (Z) ===\n');

fig_b = figure('Units', 'centimeters', 'Position', [1, 1, 13, 12.5]);
hold on;

for u = 1:length(z_levels)
    z_level = z_levels(u);
    fprintf('  Processing Z level = %.2f (%d/%d)\n', z_level, u, length(z_levels));
    
    % Filter data: gamma, Z level, Y=0
    mask = data.Gamma == gamma_value & ...
           data.Z_Level == z_level & ...
           data.Y_Level == 0;
    
    if sum(mask) < 2
        fprintf('    Skipping: insufficient data points\n');
        continue;
    end
    
    % Extract data
    T_data = data(mask, :);
    T_data = sortrows(T_data, 'rx_BIOMASS_Ec_iML1515_WT_75p37M');
    
    gr = T_data.rx_BIOMASS_Ec_iML1515_WT_75p37M;
    ac = T_data.rx_EX_ac_e_f;
    
    % Detect acetate threshold
    gr_ac = detect_acetate_threshold(T_data);
    
    % Plot line
    plot(gr, ac, '-', 'Color', colors_z{u}, 'LineWidth', 4, 'HandleVisibility', 'off');
    
    % Plot markers (downsampled)
    downsample_factor = 5;
    indices = 1:downsample_factor:length(gr);
    label_text = sprintf('\\phi_Z %d%%', round(z_level * 100));
    plot(gr(indices), ac(indices), markers{u}, ...
        'Color', colors_z{u}, 'LineWidth', 1.5, ...
        'MarkerEdgeColor', colors_z{u}, 'MarkerFaceColor', 'w', ...
        'MarkerSize', 10, 'DisplayName', label_text);
    
    % Plot acetate threshold
    ymax = 45;
    plot([gr_ac, gr_ac], [0, ymax], '--', 'Color', colors_z{u}, ...
        'LineWidth', 2, 'HandleVisibility', 'off');
    
    fprintf('    λ_ac = %.3f h^{-1}\n', gr_ac);
end

% Customize Panel B
title(sprintf('Membrane Protein Expression Level\n\\gamma = %.2f', gamma_value), ...
    'FontSize', 16, 'FontWeight', 'bold');
xlabel('Growth rate (h^{-1})', 'FontSize', 14);
ylabel('Acetate flux (mmol gDW^{-1} h^{-1})', 'FontSize', 14);
xlim([0, 1]);
ylim([0, 45]);
yticks(0:5:45);
legend('Location', 'southoutside', 'NumColumns', 2, 'Box', 'off', 'FontSize', 14);
grid on;
set(gca, 'FontSize', 12);

% Save Panel B
saveas(fig_b, fullfile(output_folder, sprintf('%s_panel_B.svg', figurePrefix)), 'svg');
saveas(fig_b, fullfile(output_folder, sprintf('%s_panel_B.png', figurePrefix)), 'png');
fprintf('  Saved: %s_panel_B.svg/.png\n', figurePrefix);
close(fig_b);

%% COMPLETION
fprintf('\n=== FIGURE 7 COMPLETE ===\n');
fprintf('Generated figures:\n');
fprintf('  Model: MAFBA (pre-filtered data)\n');
fprintf('  Gamma: %.2f\n', gamma_value);
fprintf('  Panel A: Y protein levels (cytosolic) - %d levels\n', length(y_levels));
fprintf('  Panel B: Z protein levels (membrane) - %d levels\n', length(z_levels));
fprintf('  Total panels: 2\n');
fprintf('  Total files: 2 panels × 2 formats = 4 files\n');
fprintf('  Output location: %s\n', output_folder);

%% HELPER FUNCTION - ACETATE THRESHOLD DETECTION
function gr_ac = detect_acetate_threshold(T_data)
    % Detects acetate threshold (λ_ac) from membrane flux slope change
    % Uses ATPS4rpp_f as membrane flux proxy
    
    memb = T_data.rx_ATPS4rpp_f;
    gr = T_data.rx_BIOMASS_Ec_iML1515_WT_75p37M;
    ac = T_data.rx_EX_ac_e_f;
    
    % Calculate slope of membrane flux
    dif_memb = diff(memb);
    dif_gr = diff(gr);
    slope = dif_memb ./ dif_gr;
    max_slope = max(slope);
    
    % Find where slope drops significantly (membrane saturation)
    memb_sat = find(slope <= max_slope * 0.01);
    t_i = find(ac(memb_sat) >= 0.01, 1, 'first');
    ac_threshold_idx = memb_sat(t_i);
    
    % Fallback strategies if not found
    if isempty(ac_threshold_idx)
        memb_sat = find(slope <= max_slope * 0.05);
        t_i = find(ac(memb_sat) >= 0, 1, 'first');
        ac_threshold_idx = memb_sat(t_i);
    end
    
    if isempty(ac_threshold_idx)
        memb_sat = find(slope <= max_slope * 0.2);
        t_i = find(ac(memb_sat) >= 0, 1, 'first');
        ac_threshold_idx = memb_sat(t_i);
    end
    
    if isempty(ac_threshold_idx)
        ac_threshold_idx = 1;
    end
    
    gr_ac = gr(ac_threshold_idx);
end
