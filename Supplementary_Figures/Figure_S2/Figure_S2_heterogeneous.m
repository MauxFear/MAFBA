%% SUPPLEMENTARY FIGURE S2
% Flux simulations using the heterogeneous approach.
%
% Simulation at gamma = 0.25 using glucose as carbon source.
% Shows average fluxes from 1000 simulations using random protein cost
% per reaction to represent cellular heterogeneity.
%
% Data source: pre-computed averaged flux tables (MATLAB .mat files)
%   data/outputs/Heterogeneous/wc_iML_rand_table_g_025.mat
%
% Output: PNG saved to ./output/Figure_S2_heterogeneous_g025.png
%
% See also: Heterogeneous_plot.m (plots all gamma values)

%% PATHS
baseDir  = fileparts(mfilename('fullpath'));
repoRoot = fullfile(baseDir, '..', '..');
dataDir  = fullfile(repoRoot, 'data', 'outputs', 'Heterogeneous');
outputDir = fullfile(baseDir, 'output');
if ~exist(outputDir, 'dir'), mkdir(outputDir); end

%% LOAD DATA — gamma = 0.25
dataFile = fullfile(dataDir, 'wc_iML_rand_table_g_025.mat');
if ~isfile(dataFile)
    error('Data file not found: %s', dataFile);
end
loaded = load(dataFile);
T = loaded.T;                       % MATLAB table: gr, ac, glc, memb, akgdh, mals, edd, co2, ...
T = sortrows(T, 'gr', 'ascend');    % ensure sorted by growth rate

fprintf('Loaded %d data points for gamma = 0.25\n', height(T));

%% ACETATE THRESHOLD DETECTION
ac_t = find(T.ac >= 0.01, 1, 'first');
if isempty(ac_t)
    ac_t = 1;
end
gr_ac = T.gr(ac_t);
fprintf('Acetate threshold (lambda_ac) = %.4f h^-1\n', gr_ac);

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

%% PLOT
fig = figure('Position', [10 10 800 600]);
hold on;

l1 = plot(T.gr, T.glc,       '-o', 'Color', color{1}, 'LineWidth', 1, ...
          'MarkerFaceColor', color{1}, 'MarkerSize', 8);
l3 = plot(T.gr, T.akgdh,     '-^', 'Color', color{3}, 'LineWidth', 1, ...
          'MarkerFaceColor', color{3}, 'MarkerSize', 8);
l4 = plot(T.gr, 0.5*T.memb,  '-d', 'Color', color{4}, 'LineWidth', 1, ...
          'MarkerFaceColor', color{4}, 'MarkerSize', 8);
l5 = plot(T.gr, T.mals,      '-v', 'Color', color{5}, 'LineWidth', 1, ...
          'MarkerFaceColor', color{5}, 'MarkerSize', 8);
l6 = plot(T.gr, T.edd,       '->', 'Color', color{6}, 'LineWidth', 1, ...
          'MarkerFaceColor', color{6}, 'MarkerSize', 8);
l7 = plot(T.gr, T.co2,       '-<', 'Color', color{7}, 'LineWidth', 1, ...
          'MarkerFaceColor', color{7}, 'MarkerSize', 8);
l2 = plot(T.gr, T.ac,        '-s', 'Color', color{2}, 'LineWidth', 1, ...
          'MarkerFaceColor', color{2}, 'MarkerSize', 8);

y = ylim;
lac = plot([gr_ac, gr_ac], [0, y(2)], '--', 'Color', color{2}, 'LineWidth', 1);

title({'MAFBA: Flux profile vs growth rate', ...
       'Heterogeneous weights (1000 simulations)', ...
       'varying w_{Cc},  \gamma = 0.25'});
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

%% SAVE
outFile = fullfile(outputDir, 'Figure_S2_heterogeneous_g025.png');
saveas(fig, outFile, 'png');
fprintf('Saved: %s\n', outFile);
