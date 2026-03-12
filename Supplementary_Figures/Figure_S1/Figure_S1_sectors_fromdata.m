%% SUPPLEMENTARY FIGURE S1 - HOMOGENEOUS PROTEOME SECTOR ALLOCATION
% Homogeneous sector-allocation workflow.

%% CONFIGURATION
% Set file paths
baseDir = fileparts(mfilename('fullpath'));
repoRoot = fullfile(baseDir, '..', '..');
addpath(fullfile(repoRoot, 'utils'));

% Data folder
datafolder = fullfile(repoRoot, 'data', 'outputs', 'Sectors');
fprintf('Data folder: %s\n', datafolder);

% Output folder
output_folder = fullfile(baseDir, 'output');
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end
fprintf('Output folder: %s\n', output_folder);

% Data files
modelV = 0;
dateVersion = 'Homogeneous';
sectorsfilename = sprintf('%s/sectors_modelV%i_kcatscomb_Homog.csv', datafolder, modelV);
datafilename = sprintf('%s/data_modelV%i_kcatscomb_Homog.csv', datafolder, modelV);
%% LOAD DATA
fprintf('\nLoading data...\n');
tableSectors = readtable(sectorsfilename);
tableData = readtable(datafilename);
fprintf('Sectors data: %d rows\n', height(tableSectors));
fprintf('Flux data: %d rows\n', height(tableData));

% Get unique gamma values
gamma_values = unique(tableSectors.Gamma);
fprintf('Gamma values to plot: %s\n', mat2str(gamma_values'));
% gamma_values = [1]; % Uncomment to plot single gamma

% Define colors for the plot
color = {sscanf('65AF65','%2x%2x%2x',[1 3])/255, ...
         sscanf('C74848','%2x%2x%2x',[1 3])/255, ...
         sscanf('9065D1','%2x%2x%2x',[1 3])/255, ...
         sscanf('266DD7','%2x%2x%2x',[1 3])/255, ...
         sscanf('F2A1E4','%2x%2x%2x',[1 3])/255, ...
         sscanf('D2C725','%2x%2x%2x',[1 3])/255, ...
         sscanf('B2B2B2','%2x%2x%2x',[1 3])/255, ...
         sscanf('EC8327','%2x%2x%2x',[1 3])/255, ...
         sscanf('35B2B2','%2x%2x%2x',[1 3])/255,...
         sscanf('E40303','%2x%2x%2x',[1 3])/255};%sligthly bright red;

%% GENERATE FIGURES
fprintf('\nGenerating figures for %d gamma values...\n', length(gamma_values));

for f = 1:length(gamma_values)
    gamma_val = gamma_values(f);
    gamma_string = ['\gamma = ', num2str(gamma_val)];
    fprintf('\n  Processing gamma = %.2f...\n', gamma_val);
    
    % Filter the table data for the current gamma value
    T = tableData(tableData.Gamma == gamma_val, :);
    
    % Sort rows by growth rate
    T = sortrows(T, 'BIOMASS_Ec_iML1515_WT_75p37M');
    
    % Detect acetate threshold >=0.1
    dif_v = diff(T.ATPS4rpp_f);
    dif_g = diff(T.BIOMASS_Ec_iML1515_WT_75p37M);
    slope_v = dif_v ./ dif_g;
    m_s = max(slope_v);
    
    memb_s = find(slope_v <= m_s * 0.01);
    t_i = find(T.EX_ac_e_f(memb_s) >= 0.01, 1, 'first');
    ac_t = memb_s(t_i);
    
    if isempty(ac_t)
        memb_s = find(slope_v <= m_s * 0.05);
        t_i = find(T.EX_ac_e_f(memb_s) >= 0, 1, 'first');
        ac_t = memb_s(t_i);
    end
    if isempty(ac_t)
        memb_s = find(slope_v <= m_s * 0.2);
        t_i = find(T.EX_ac_e_f(memb_s) >= 0, 1, 'first');
        ac_t = memb_s(t_i);
    end
    if isempty(ac_t)
        ac_t = 1;
    end
    
    gr_ac = T.BIOMASS_Ec_iML1515_WT_75p37M(ac_t);
    
    % Filter the table data for the current gamma value
    S = tableSectors(tableSectors.Gamma == gamma_val, :);
    
    % Sort rows by growth rate
    S = sortrows(S, 'Growth');
    
    % Calculate phiPmax and phiMmax
    phiPmax = max(S.phiR + S.phiEc + S.phiEr + S.phiCm + S.phiEm);
    phiMmax = gamma_val*phiPmax;
    
    % Plot the data
    % fig = get_plot(S, gamma_string, color, phiPmax, phiMmax, 
    % gr_ac);
    [fig, h] = get_plot(S, gamma_string, color, phiPmax, phiMmax, gr_ac);
    
    
    
    % Save the figure as SVG
    filename_svg = sprintf('%s/%s_sectors_g%03d.svg', output_folder, dateVersion, round(gamma_val * 100));
    saveas(fig, filename_svg, 'svg');
    fprintf('    Saved: %s\n', filename_svg);
    
    % Save the figure as PNG
    filename_png = sprintf('%s/%s_sectors_g%03d.png', output_folder, dateVersion, round(gamma_val * 100));
    saveas(fig, filename_png, 'png');
    fprintf('    Saved: %s\n', filename_png);
    
    % Apply hatching and save the hatched image
    % colors = [color{1}; color{6}; color{4}; color{3}; color{2}; 
    % color{8}];
    % im_hatch = applyhatch_pluscolor(fig, '|\x+c/', 1, [1 1 1 1 0 
    % 1], colors, 300, 5, 5);
    hatch_mask = [1 1 1 1 0 1];
    hatch_angles = [-45, 45, -45, 90, 0, 45]; % Use different angles for different patterns
    
    hold on;
    for i = 1:length(h)
        if hatch_mask(i) == 1
            hatchfill2(h(i), 'single', 'HatchAngle', hatch_angles(i), 'HatchDensity', 40, 'HatchColor', 'black');
        end
    end
    hold off;
    
    filename_hatch = sprintf('%s/%s_sectors_g%03d_hatch.png', output_folder, dateVersion, round(gamma_val * 100));
    saveas(fig, filename_hatch, 'png');
    fprintf('    Saved: %s\n', filename_hatch);
    
end

fprintf('\n✓ Supplementary Figure S1 generation complete!\n');
fprintf('  Total figures: %d gamma values × 3 formats = %d files\n', length(gamma_values), length(gamma_values)*3);


function [fig, h] = get_plot(T, gamma_string, color, phiPmax, phiMmax, gr_ac)
% Increase figure size (adjust these values as needed)
    fig = figure('Units', 'centimeters', 'Position', [1, 1, 18, 14]);
    % Create axes with specific position (adjust these values as needed)
    ax = axes('Position', [0.1, 0.1, 0.7, 0.8]);
    hold on;

%     Set the title, xlabel, and ylabel properties
    title(gamma_string);
    xlabel('Growth rate (h^{-1})', 'FontSize', 12);
    ylabel('Cumulative proteome fraction', 'FontSize', 12);

%     Extract relevant columns
    gr = T.Growth;
%     Calculate phiCc as phiPmax minus the sum of all other sectors
    phiCc = phiPmax - (T.phiR + T.phiEc + T.phiEr + T.phiCm + T.phiEm);
    gdata = [T.phiCm, T.phiEm, T.phiEr, T.phiR, T.phiEc, phiCc];

    % Plot the area plot
    h = area(gr, gdata);

    % Set colors and DisplayName for each sector
    set(h(1), 'FaceColor', color{1}, 'LineWidth', 1, 'DisplayName', '\phi C_{m}');
    set(h(2), 'FaceColor', color{6}, 'LineWidth', 1, 'DisplayName', '\phi E_{m}');
    set(h(3), 'FaceColor', color{4}, 'LineWidth', 1, 'DisplayName', '\phi E_{r}');
    set(h(4), 'FaceColor', color{3}, 'LineWidth', 1, 'DisplayName', '\phi R');
    set(h(5), 'FaceColor', color{2}, 'LineWidth', 1, 'DisplayName', '\phi E_{c}');
    set(h(6), 'FaceColor', color{8}, 'LineWidth', 1, 'DisplayName', '\phi C_{c}');

    xlim([0, 1]);
    y = ylim;

    % Plot the additional line
    h8 = plot([gr_ac, gr_ac], [0, y(2)], '--', 'Color', color{10}, 'LineWidth', 2, 'MarkerFaceColor', color{2}, 'DisplayName', '\lambda_{ac}');
    ylim(y);

    % Set legend
%     legendHandles = [h(1), h(2), h(3), h(4), h(5), h(6), h8];
%     legendLabels = {'\phi C_{m}', '\phi E_{m}', '\phi E_{r}', '\phi R', '\phi E_{c}', '\phi C_{c}', '\lambda_{ac}'};
%     lgd = legend(legendHandles, legendLabels);
    lgd = legend([h(1), h(2), h(3), h(4), h(5), h(6), h8]);
    set(lgd, 'Box', 'off', 'Location', 'southoutside', 'NumColumns', 4, 'FontSize', 20);    set(lgd, 'Box', 'off', 'Location', 'southoutside', 'NumColumns', 4, 'FontSize', 20);

    % Add horizontal dashed lines for phiPmax and phiMmax
    hLine1 = line([gr(1), gr(end)], [phiPmax, phiPmax], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 2.5);
    hLine2 = line([gr(1), gr(end)], [phiMmax, phiMmax], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 2.5);

    % Add annotations
    threshold = 0.05; % Adjust this value to define "close"
    if abs(phiPmax - phiMmax) <= threshold
        % Values are close or equal, adjust text placement
        text(gr(end)*1.015, phiPmax, ['\phi P_{max} = ', num2str(phiPmax, '%.2f')], 'VerticalAlignment', 'bottom');
        text(gr(end)*1.015, phiMmax, ['\phi M_{max} = ', num2str(phiMmax, '%.2f')], 'VerticalAlignment', 'top');
    else
        % Values are not close, use original placement
        text(gr(end)*1.015, phiPmax, ['\phi P_{max} = ', num2str(phiPmax, '%.2f')], 'VerticalAlignment', 'middle');
        text(gr(end)*1.015, phiMmax, ['\phi M_{max} = ', num2str(phiMmax, '%.2f')], 'VerticalAlignment', 'middle');
    end
    
    
%     Adjust axis and grid
    axis([0 gr(end) 0 phiPmax]);
    set(gca, 'FontSize', 12);
    grid on;
    
%     Remove horizontal lines from legend
    set(get(get(hLine1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    set(get(get(hLine2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
end
