%% HETEROGENEOUS FLUX PROFILE — ALL GAMMA VALUES
% Average metabolic fluxes vs growth rate using the heterogeneous approach.
%
% Each panel shows average fluxes from 1000 simulations with random
% protein cost per reaction at a specific gamma value.
%
% Data source: data/outputs/Heterogeneous/wc_iML_rand_table_g_*.mat
% Output: PNG per gamma saved to ./output/

%% PATHS
baseDir  = fileparts(mfilename('fullpath'));
repoRoot = fullfile(baseDir, '..', '..');
dataDir  = fullfile(repoRoot, 'data', 'outputs', 'Heterogeneous');
outputDir = fullfile(baseDir, 'output');
if ~exist(outputDir, 'dir'), mkdir(outputDir); end

%% LOAD DATA
T_g_1   = load(fullfile(dataDir, 'wc_iML_rand_table_g_1.mat'));
T_g_025 = load(fullfile(dataDir, 'wc_iML_rand_table_g_025.mat'));
T_g_023 = load(fullfile(dataDir, 'wc_iML_rand_table_g_023.mat'));
T_g_02  = load(fullfile(dataDir, 'wc_iML_rand_table_g_02.mat'));
T_g_015 = load(fullfile(dataDir, 'wc_iML_rand_table_g_015.mat'));

gamma_v = [1, 0.25, 0.23, 0.2, 0.15];
T_v     = {T_g_1.T, T_g_025.T, T_g_023.T, T_g_02.T, T_g_015.T};

%% COLORS
color = {
    sscanf('65AF65','%2x%2x%2x',[1 3])/255, ...  % green  — glucose
    sscanf('C74848','%2x%2x%2x',[1 3])/255, ...  % red    — acetate
    sscanf('9065D1','%2x%2x%2x',[1 3])/255, ...  % purple — AKGDH
    sscanf('266DD7','%2x%2x%2x',[1 3])/255, ...  % blue   — membrane
    sscanf('F2A1E4','%2x%2x%2x',[1 3])/255, ...  % pink   — glyoxylate shunt
    sscanf('D2C725','%2x%2x%2x',[1 3])/255, ...  % yellow — ED pathway
    sscanf('B2B2B2','%2x%2x%2x',[1 3])/255, ...  % gray   — CO2
    sscanf('EC8327','%2x%2x%2x',[1 3])/255, ...  % orange — wCc
};

%% PLOT EACH GAMMA
for g = 1:length(gamma_v)
    gammaVal     = gamma_v(g);
    gamma_string = ['\gamma = ', num2str(gammaVal)];

    T = T_v{g};
    T = sortrows(T, 'gr', 'ascend');

    % Acetate threshold detection
    ac_t = find(T.ac >= 0.01, 1, 'first');
    if isempty(ac_t), ac_t = 1; end
    gr_ac = T.gr(ac_t);
    fprintf('gamma = %.2f  ->  lambda_ac = %.4f h^-1\n', gammaVal, gr_ac);

    % Figure
    fig = figure('Position', [10 10 800 600]);
    hold on;
    set(fig, 'defaultAxesColorOrder', [color{4}; color{8}]);

    l1 = plot(T.gr, T.glc,      '-o', 'Color', color{1}, 'LineWidth', 1, 'MarkerFaceColor', color{1}, 'MarkerSize', 8);
    l3 = plot(T.gr, T.akgdh,    '-^', 'Color', color{3}, 'LineWidth', 1, 'MarkerFaceColor', color{3}, 'MarkerSize', 8);
    l4 = plot(T.gr, 0.5*T.memb, '-d', 'Color', color{4}, 'LineWidth', 1, 'MarkerFaceColor', color{4}, 'MarkerSize', 8);
    l5 = plot(T.gr, T.mals,     '-v', 'Color', color{5}, 'LineWidth', 1, 'MarkerFaceColor', color{5}, 'MarkerSize', 8);
    l6 = plot(T.gr, T.edd,      '->', 'Color', color{6}, 'LineWidth', 1, 'MarkerFaceColor', color{6}, 'MarkerSize', 8);
    l7 = plot(T.gr, T.co2,      '-<', 'Color', color{7}, 'LineWidth', 1, 'MarkerFaceColor', color{7}, 'MarkerSize', 8);
    l2 = plot(T.gr, T.ac,       '-s', 'Color', color{2}, 'LineWidth', 1, 'MarkerFaceColor', color{2}, 'MarkerSize', 8);

    y   = ylim;
    lac = plot([gr_ac, gr_ac], [0, y(2)], '--', 'Color', color{2}, 'LineWidth', 1);

    title({'MAFBA: Flux profile vs growth rate', 'Heterogeneous weights', ...
           ['varying w_{Cc}   ', gamma_string]});
    xlabel('Growth rate (h^{-1})', 'FontSize', 18);
    ylabel('Flux (mmol/g_{DW}h)', 'FontSize', 18);
    set(gca, 'FontSize', 16);
    axis([0 inf 0 inf]);

    lgd = legend([l1; l2; l3; l4; l5; l6; l7; lac], ...
        'Glucose uptake', 'Acetate excretion', 'AKGDH (TCA)', ...
        'ATPS4r flux 0.5\times (Membrane)', 'Glyoxylate shunt', ...
        'ED pathway', 'CO_{2} excretion', '\lambda_{ac}');
    set(lgd, 'Box', 'off', 'Location', 'northwest', 'FontSize', 16);
    hold off;

    % Save
    str_g    = strrep(num2str(gammaVal), '.', '');
    filename = fullfile(outputDir, ['heterogeneous_flux_g', str_g, '.png']);
    saveas(fig, filename, 'png');
    fprintf('Saved: %s\n', filename);
    close(fig);
end

fprintf('\nAll gamma plots saved to: %s\n', outputDir);
