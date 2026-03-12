%% SUPPLEMENTARY FIGURE S10
% Correlation between acetate threshold and normalized heterologous protein level.
%
% Panels:
%   A - Cytosolic protein (phi_Z) normalized by total proteome capacity (phi_P_max)
%   B - Membrane-associated protein (phi_Y) normalized by phi_P_max
%   C - Membrane-associated protein (phi_Y) normalized by membrane capacity (phi_M_max)
%   D - Annotated contour for cytosolic protein showing correlation regions
%
% Data source: pre-computed acetate thresholds (aceL_protein_expression_combined.csv)
% Output: SVG files saved to ./output/
%
%% PATHS
baseDir  = fileparts(mfilename('fullpath'));
repoRoot = fullfile(baseDir, '..', '..');

dataFile = fullfile(repoRoot, 'data', 'outputs', 'Protein_Expression', ...
                   'aceL_protein_expression_combined.csv');
if ~isfile(dataFile)
    error('Data file not found: %s', dataFile);
end

outputDir = fullfile(baseDir, 'output');
if ~exist(outputDir, 'dir'), mkdir(outputDir); end

%% LOAD & PRE-PROCESS DATA
data = readtable(dataFile);
fprintf('Loaded %d rows from aceL data.\n', height(data));

% Use only baseline simulations (Change_Factor == 1)
data = data(data.Change_Factor == 1, :);

phi_max = 0.484;  % total proteome capacity (phi_P_max)

% ---- Z data (cytosolic protein): Y_Level == 0 ----------------------------
z_raw = data(data.Y_Level == 0, {'Gamma', 'Z_Level', 'Acetate_threshold'});
z_raw.Properties.VariableNames = {'gamma', 'phi_Z', 't_ac'};
z_raw.Z_divided_by_max = z_raw.phi_Z / phi_max;   % phi_Z / phi_P_max

% ---- Y data (membrane protein): Z_Level == 0 -----------------------------
y_raw = data(data.Z_Level == 0, {'Gamma', 'Y_Level', 'Acetate_threshold'});
y_raw.Properties.VariableNames = {'gamma', 'phi_Y', 't_ac'};
y_raw.Y_divided_by_max  = y_raw.phi_Y / phi_max;               % phi_Y / phi_P_max
y_raw.Y_divided_by_Mmax = y_raw.phi_Y ./ (y_raw.gamma * phi_max); % phi_Y / phi_M_max

fprintf('Z data: %d points   Y data: %d points\n', height(z_raw), height(y_raw));

% Build concatenated Y table filtered to the feasible phi_Y/phi_M_max <= 1 region
% and clip to xlim [0,1].
concatenatedTable = table();
ygammas = unique(y_raw.gamma);
for i = 1:numel(ygammas)
    gammaVal = ygammas(i);
    sub = y_raw(y_raw.gamma == gammaVal, :);
    sub = sortrows(sub, 'Y_divided_by_Mmax');
    % Keep only feasible normalised range (phi_Y / phi_M_max <= 1)
    sub = sub(sub.Y_divided_by_Mmax <= 1, :);
    if height(sub) == 0, continue; end
    concatenatedTable = vertcat(concatenatedTable, sub);
end
fprintf('Concatenated Y/Mmax data: %d points\n', height(concatenatedTable));

%% COLORMAPS
% --- Blue-green gradient (Y panels, n=20 per segment) — Panels A-like Y ---
white_c  = [1 1 1];
green_c  = sscanf('6B8E23','%2x%2x%2x',[1 3])/255;
aqua_c   = sscanf('40E0D0','%2x%2x%2x',[1 3])/255;
blue_c   = sscanf('6495ED','%2x%2x%2x',[1 3])/255;
dblue_c  = sscanf('1F56A8','%2x%2x%2x',[1 3])/255;

n20 = 20;
Y_cmap20 = [
    interp1Colors(white_c, green_c, n20);
    interp1Colors(green_c, aqua_c,  n20);
    interp1Colors(aqua_c,  blue_c,  n20);
    dblue_c
];  % 61 entries for Panels B

n30 = 30;
Y_cmap30 = [
    interp1Colors(white_c, green_c, n30);
    interp1Colors(green_c, aqua_c,  n30);
    interp1Colors(aqua_c,  blue_c,  n30);
    dblue_c
];  % 91 entries for Panel C (NormContour, n=30)

% --- Red gradient (Z panel, n=20 per segment) ---
yellow_c = sscanf('FFD700','%2x%2x%2x',[1 3])/255;
orange_c = sscanf('FF6F61','%2x%2x%2x',[1 3])/255;
red_c    = sscanf('DC143C','%2x%2x%2x',[1 3])/255;
dred_c   = sscanf('C00217','%2x%2x%2x',[1 3])/255;

Z_cmap = [
    interp1Colors(white_c,  yellow_c, n20);
    interp1Colors(yellow_c, orange_c, n20);
    interp1Colors(orange_c, red_c,    n20);
    dred_c
];  % 61 entries for Panel A

%% PANEL A — Cytosolic protein / phi_P_max  (Contour_uProtZvsAceT)
% x=Z/phiMax, y=gamma, 100x100 grid, 20 contour levels
fig_A = contourPanel( ...
    z_raw.Z_divided_by_max, z_raw.gamma, z_raw.t_ac, ...
    '\phi Z /\phi_{max}', ...
    '\gamma ratio (\phi^{m}_{max}/\phi_{max})', ...
    Z_cmap, [0 1]);
saveas(fig_A, fullfile(outputDir, 'panel_A_cytosolic_phiPmax.svg'), 'svg');
fprintf('Saved panel A\n');

%% PANEL B — Membrane protein / phi_P_max  (Contour_uProtYvsAceT)
% x=Y/phiMax, y=gamma, 100x100 grid, 20 contour levels
fig_B = contourPanel( ...
    y_raw.Y_divided_by_max, y_raw.gamma, y_raw.t_ac, ...
    '\phi Y /\phi_{max}', ...
    '\gamma ratio (\phi^{m}_{max}/\phi_{max})', ...
    Y_cmap20, [0 1]);
saveas(fig_B, fullfile(outputDir, 'panel_B_membrane_phiPmax.svg'), 'svg');
fprintf('Saved panel B\n');

%% PANEL C — Membrane protein / phi_M_max  (NormContour_uProtYvsAceT)
% x=Y/phiMmax (filtered, <=1), y=gamma, 100x100 grid, 20 contour levels, n=30 cmap
fig_C = contourPanel( ...
    concatenatedTable.Y_divided_by_Mmax, concatenatedTable.gamma, concatenatedTable.t_ac, ...
    '\phi Y /\phi^{m}_{max}', ...
    '\gamma ratio (\phi^{m}_{max}/\phi_{max})', ...
    Y_cmap30, [0 1]);
saveas(fig_C, fullfile(outputDir, 'panel_C_membrane_phiMmax.svg'), 'svg');
fprintf('Saved panel C\n');

%% PANEL D — Annotated cytosolic contour  (same data as A with region labels)
fig_D = contourPanel( ...
    z_raw.Z_divided_by_max, z_raw.gamma, z_raw.t_ac, ...
    '\phi Z /\phi_{max}', ...
    '\gamma ratio (\phi^{m}_{max}/\phi_{max})', ...
    Z_cmap, [0 1]);
ax_D = gca;

% Horizontal dashed line at gamma = 0.25 (transition boundary)
yline(ax_D, 0.25, '--w', 'LineWidth', 1.5);

% Region 1 (gamma > 0.25): threshold depends only on phi_Z, not on gamma
text(ax_D, 0.35, 0.62, ...
    {'Linear correlation with \phi_Z', 'No effect of \gamma'}, ...
    'HorizontalAlignment', 'center', 'FontSize', 9, ...
    'Color', 'w', 'FontWeight', 'bold');

% Region 2 (gamma < 0.25): threshold depends on both phi_Z and gamma
text(ax_D, 0.18, 0.13, ...
    {'Linear correlation', 'with \phi_Z and \gamma'}, ...
    'HorizontalAlignment', 'center', 'FontSize', 9, ...
    'Color', 'w', 'FontWeight', 'bold');

% Infeasible region (upper right: high phi_Z forces the model to be infeasible)
text(ax_D, 0.82, 0.50, ...
    {'Infeasible', 'region'}, ...
    'HorizontalAlignment', 'center', 'FontSize', 9, ...
    'Color', [0.3 0.3 0.3], 'FontWeight', 'bold');

saveas(fig_D, fullfile(outputDir, 'panel_D_annotated.svg'), 'svg');
fprintf('Saved panel D\n');

fprintf('\nFigure S10 complete. SVGs saved to: %s\n', outputDir);

%% LOCAL FUNCTIONS

function fig = contourPanel(xVec, yVec, cVec, xlbl, ylbl, cmap, xlims)
    % Filled contour using griddata
    fig = figure();

    x_grid = linspace(min(xVec), max(xVec), 100);   % 100x100 grid
    y_grid = linspace(min(yVec), max(yVec), 100);
    [X, Y]  = meshgrid(x_grid, y_grid);
    Z_color = griddata(xVec, yVec, cVec, X, Y);      % default (linear) interpolation

    colormap(cmap);
    contourf(X, Y, Z_color, 20, 'LineStyle', ':', 'LineWidth', 1);  % 20 levels
    colormap(cmap);   % re-apply after contourf
    cb = colorbar;
    ylabel(cb, 'Acetate Threshold (\lambda_{ac}) (h^{-1})', 'FontSize', 14);
    set(gca, 'FontSize', 12);
    xlim(xlims);
    xlabel(xlbl, 'FontSize', 14);
    ylabel(ylbl, 'FontSize', 14);
end

function rgb = interp1Colors(c1, c2, n)
    rgb = [linspace(c1(1), c2(1), n);
           linspace(c1(2), c2(2), n);
           linspace(c1(3), c2(3), n)].';
end
