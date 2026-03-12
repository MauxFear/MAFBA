% FIGURE 8 - NORMALIZED PROTEIN EXPRESSION LEVELS
%
% Creates four plots:
% - Panel B: Membrane (Z) normalized by phi_max^P
% - Panel C: Cytosolic (Y) normalized by phi_max^P
% - Panel D: Cytosolic (Y) normalized by phi_max^M
% - Combined: Both Y and Z on same plot (no panel letter)

% Enable validation warnings for acetate threshold detection
global ENABLE_VALIDATION_WARNINGS;
ENABLE_VALIDATION_WARNINGS = false;  % Set to true for debugging

% Set up paths
baseDir = fileparts(mfilename('fullpath'));
repoRoot = fullfile(baseDir, '..', '..');
datafolder = fullfile(repoRoot, 'data', 'outputs', 'Protein_Expression');

% Load the filtered data
data_file = fullfile(datafolder, 'data_protein_expression_figure8.csv');
if ~isfile(data_file)
    error('Filtered data file not found at: %s\nRun filter_data_for_figure8.py to create it', data_file);
end

% Read the data
data = readtable(data_file);
fprintf('Loaded data with %d rows\n', height(data));

% Map reaction names to expected column names with rx_ prefix (same as Figure_4)
reaction_map = struct();
if ismember('rx_EX_glc__D_e_b', data.Properties.VariableNames)
    reaction_map.glc = 'rx_EX_glc__D_e_b';
elseif ismember('glc', data.Properties.VariableNames)
    reaction_map.glc = 'glc';
else
    error('Glucose uptake column not found');
end

if ismember('rx_EX_ac_e_f', data.Properties.VariableNames)
    reaction_map.ac = 'rx_EX_ac_e_f';
elseif ismember('ac', data.Properties.VariableNames)
    reaction_map.ac = 'ac';
else
    error('Acetate excretion column not found');
end

if ismember('rx_AKGDH', data.Properties.VariableNames)
    reaction_map.akgdh = 'rx_AKGDH';
elseif ismember('akgdh', data.Properties.VariableNames)
    reaction_map.akgdh = 'akgdh';
else
    error('AKGDH column not found');
end

if ismember('rx_ATPS4rpp_f', data.Properties.VariableNames)
    reaction_map.memb = 'rx_ATPS4rpp_f';
elseif ismember('memb', data.Properties.VariableNames)
    reaction_map.memb = 'memb';
else
    error('Membrane/ATPS4rpp column not found');
end

if ismember('rx_MALS', data.Properties.VariableNames)
    reaction_map.mals = 'rx_MALS';
elseif ismember('mals', data.Properties.VariableNames)
    reaction_map.mals = 'mals';
else
    error('MALS column not found');
end

if ismember('rx_EDD', data.Properties.VariableNames)
    reaction_map.edd = 'rx_EDD';
elseif ismember('edd', data.Properties.VariableNames)
    reaction_map.edd = 'edd';
else
    error('EDD column not found');
end

if ismember('rx_EX_co2_e_f', data.Properties.VariableNames)
    reaction_map.co2 = 'rx_EX_co2_e_f';
elseif ismember('co2', data.Properties.VariableNames)
    reaction_map.co2 = 'co2';
else
    error('CO2 excretion column not found');
end

% Find growth rate column
if ismember('rx_BIOMASS_Ec_iML1515_WT_75p37M', data.Properties.VariableNames)
    reaction_map.gr = 'rx_BIOMASS_Ec_iML1515_WT_75p37M';
elseif ismember('rx_BIOMASS_Ec_iML1515_core_75p37M', data.Properties.VariableNames)
    reaction_map.gr = 'rx_BIOMASS_Ec_iML1515_core_75p37M';
elseif ismember('gr', data.Properties.VariableNames)
    reaction_map.gr = 'gr';
else
    error('Growth rate column not found');
end

% Get unique values for filtering
gamma_values = unique(data.Gamma);
z_levels = unique(data.Z_Level);
y_levels = unique(data.Y_Level);

fprintf('Found %d gamma values, %d Z levels, and %d Y levels\n', ...
    length(gamma_values), length(z_levels), length(y_levels));

% Process data to create arrays similar to original Figure_5
phi_max = 0.484;
gamma_select = [1, 0.35, 0.3, 0.25, 0.22, 0.2];

% Initialize arrays for Z and Y data
triArray_Z = [];
triArray_Y = [];

% Process Z data (cytosolic proteins)
fprintf('Processing Z (cytosolic) data...\n');
for g = 1:length(gamma_select)
    gamma_val = gamma_select(g);
    
    % Get unique Z levels for this gamma
    mask_gamma = data.Gamma == gamma_val & ...
                 data.Y_Level == 0;
    
    gamma_data = data(mask_gamma, :);
    unique_z_levels = unique(gamma_data.Z_Level);
    unique_z_levels = unique_z_levels(unique_z_levels >= 0); % Include zero level (no cytosolic proteins)
    
    for z = 1:length(unique_z_levels)
        z_level = unique_z_levels(z);
        
        % Filter data for current Z level and gamma
        mask = data.Gamma == gamma_val & ...
               data.Z_Level == z_level & ...
               data.Y_Level == 0;
        
        if sum(mask) < 2
            continue;
        end
        
        % Get data subset
        T_data = data(mask, :);
        T_data = sortrows(T_data, reaction_map.gr);
        
        % Create table with same column names as Figure_4
        T = table();
        T.gr = T_data.(reaction_map.gr);
        T.ac = T_data.(reaction_map.ac);
        T.glc = T_data.(reaction_map.glc);
        T.memb = T_data.(reaction_map.memb);
        T.akgdh = T_data.(reaction_map.akgdh);
        T.mals = T_data.(reaction_map.mals);
        T.edd = T_data.(reaction_map.edd);
        T.co2 = T_data.(reaction_map.co2);
        
        % Get acetate threshold using same function as Figure_4
        gr_ac = get_acetateT(T);
        
        % Add to triArray_Z: [gamma, phi_Z, t_ac, Z_divided_by_max]
        triArray_Z = [triArray_Z; gamma_val, z_level, gr_ac, z_level/phi_max];
    end
end

% Process Y data (membrane proteins)
fprintf('Processing Y (membrane) data...\n');
for g = 1:length(gamma_select)
    gamma_val = gamma_select(g);
    
    % Get unique Y levels for this gamma
    mask_gamma = data.Gamma == gamma_val & ...
                 data.Z_Level == 0;
    
    gamma_data = data(mask_gamma, :);
    unique_y_levels = unique(gamma_data.Y_Level);
    unique_y_levels = unique_y_levels(unique_y_levels >= 0); % Include zero level (no membrane proteins)
    
    for y = 1:length(unique_y_levels)
        y_level = unique_y_levels(y);
        
        % Filter data for current Y level and gamma
        mask = data.Gamma == gamma_val & ...
               data.Y_Level == y_level & ...
               data.Z_Level == 0;
        
        if sum(mask) < 2
            continue;
        end
        
        % Get data subset
        T_data = data(mask, :);
        T_data = sortrows(T_data, reaction_map.gr);
        
        % Create table with same column names as Figure_4
        T = table();
        T.gr = T_data.(reaction_map.gr);
        T.ac = T_data.(reaction_map.ac);
        T.glc = T_data.(reaction_map.glc);
        T.memb = T_data.(reaction_map.memb);
        T.akgdh = T_data.(reaction_map.akgdh);
        T.mals = T_data.(reaction_map.mals);
        T.edd = T_data.(reaction_map.edd);
        T.co2 = T_data.(reaction_map.co2);
        
        % Get acetate threshold using same function as Figure_4
        gr_ac = get_acetateT(T);
        
        % Add to triArray_Y: [gamma, phi_Y, t_ac, Y_divided_by_max, Y_divided_by_Mmax]
        phi_m_max = gamma_val * phi_max;
        triArray_Y = [triArray_Y; gamma_val, y_level, gr_ac, y_level/phi_max, y_level/phi_m_max];
    end
end

% Create tables to manipulate the data
z_df = array2table(triArray_Z, 'VariableNames', {'gamma', 'phi_Z', 't_ac', 'Z_divided_by_max'});
z_df = sortrows(z_df, {'gamma','Z_divided_by_max'});
disp('Z data:');
disp(z_df(1:min(10,height(z_df)),:))

y_df = array2table(triArray_Y, 'VariableNames', {'gamma', 'phi_Y', 't_ac', 'Y_divided_by_max', 'Y_divided_by_Mmax'});
y_df = sortrows(y_df, {'gamma','Y_divided_by_max'});
disp('Y data:');
disp(y_df(1:min(10,height(y_df)),:))

fprintf('\nNote: Using enhanced acetate threshold detection with 10x more stringent criteria\n');
fprintf('(slope threshold 0.001 vs 0.01) for improved accuracy in transition point detection.\n');

% Add data range diagnostics
fprintf('\nData Range Diagnostics:\n');
fprintf('Z data (Cytosolic): %d total points\n', height(z_df));
if height(z_df) > 0
    fprintf('  - Acetate threshold range: %.4f to %.4f h^-1\n', min(z_df.t_ac), max(z_df.t_ac));
    fprintf('  - phi_Z range: %.4f to %.4f\n', min(z_df.Z_divided_by_max), max(z_df.Z_divided_by_max));
    fprintf('  - Points with t_ac ≤ 0.1: %d\n', sum(z_df.t_ac <= 0.1));
end

fprintf('Y data (Membrane): %d total points\n', height(y_df));
if height(y_df) > 0
    fprintf('  - Acetate threshold range: %.4f to %.4f h^-1\n', min(y_df.t_ac), max(y_df.t_ac));
    fprintf('  - phi_Y/phi_max range: %.4f to %.4f\n', min(y_df.Y_divided_by_max), max(y_df.Y_divided_by_max));
    fprintf('  - phi_Y/phi_Mmax range: %.4f to %.4f\n', min(y_df.Y_divided_by_Mmax), max(y_df.Y_divided_by_Mmax));
    fprintf('  - Points with t_ac ≤ 0.1: %d\n', sum(y_df.t_ac <= 0.1));
    fprintf('  - Points with t_ac = 0: %d (first zero preserved, rest removed for clean plots)\n', sum(y_df.t_ac == 0));
end

%%
% PLOT 1: Z-only plot (Cytosolic proteins)
fig1 = figure('Units', 'centimeters', 'Position', [1, 1, 13, 12.5]);
hold on;

markers = {'o', '^', 's', 'd', '>', 'v', '<'};
colors = {
        sscanf('FFA500','%2x%2x%2x',[1 3])/255,... % Mango Orange 
        sscanf('FF6F61','%2x%2x%2x',[1 3])/255,... % Coral 
        sscanf('FF3E3E','%2x%2x%2x',[1 3])/255,... % Tomato Red
        sscanf('DC143C','%2x%2x%2x',[1 3])/255,... % Crimson
        sscanf('C00217','%2x%2x%2x',[1 3])/255,... % dark red
        [0.07450980392156863, 0.06666666666666667, 0.0196078431372549],... %black
        sscanf('FFFFFF','%2x%2x%2x',[1 3])/255,... %white
             };
         
xlabel('\phi Z /\phi^{P}_{max}', 'FontSize', 14)
ylabel('Acetate Threshold (\lambda_{ac}) (h^{-1})', 'FontSize', 14)
set(gca, 'FontSize', 12)

% Plot Cytosolic Protein for specified gamma values
gamma_select_reversed = gamma_select(end:-1:1);

for i=1:numel(gamma_select_reversed)
    gammaVal = gamma_select_reversed(i);
    sub_z_df = z_df(z_df.gamma == gammaVal,:); 
    sub_z_df = sortrows(sub_z_df, {'gamma','Z_divided_by_max'});
    
    if height(sub_z_df) == 0
        continue;
    end
    
    % FIXED: Preserve only the first zero acetate threshold (transition point) for cleaner plots
    % Keep all data up to and including first zero, remove everything after
    indexZero = find(sub_z_df.t_ac == 0, 1);
    if ~isempty(indexZero)
        sub_z_df = sub_z_df(1:indexZero, :);  % Keep up to and including first zero
    else
        % If no zeros, just remove negative values
        validPoints = sub_z_df.t_ac >= 0;
        sub_z_df = sub_z_df(validPoints, :);
    end
    
    if height(sub_z_df) == 0
        continue;
    end
    
    colorAdjusted = [colors{i}, 1];    
    plot(sub_z_df.Z_divided_by_max, sub_z_df.t_ac, '-', 'Color', colorAdjusted, 'LineWidth', 4, 'HandleVisibility', 'off');
    
    % Downsample the data to only keep points with markers
    downsampleFactor = max(1, floor(height(sub_z_df) / 10));
    indicesToShow = 1 : downsampleFactor : height(sub_z_df);
    
    if i <= 2 && height(sub_z_df) > 20
        indicesToShow = sort([indicesToShow, height(sub_z_df)-1, height(sub_z_df)]);
    elseif height(sub_z_df) > 10
        indicesToShow = sort([indicesToShow, height(sub_z_df)-1, height(sub_z_df)]);
    end
    
    indicesToShow = unique(indicesToShow);
    indicesToShow = indicesToShow(indicesToShow <= height(sub_z_df));
    
    gammaLabel = sprintf('Cytosolic (\\phi Z) \\gamma =  %s', num2str(gammaVal)); 
    plot(sub_z_df.Z_divided_by_max(indicesToShow), sub_z_df.t_ac(indicesToShow),...
    markers{i}, 'Color', colorAdjusted, 'LineWidth', 1.5,'MarkerEdgeColor', colors{i},...
    'MarkerFaceColor', colors{end}, 'MarkerSize', 8, 'DisplayName', gammaLabel);
end

xlim([0, 1]);
% Set the legend
lgd = legend('Location', 'southoutside','NumColumns', 2);
set(lgd, 'Box', 'off', 'FontSize', 12);
% Show the grid
grid on;

% Save the figure as SVG
output_dir = fullfile(baseDir, 'output');
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end
figurePrefix = 'iML1515_MAFBA_protNorm';
filename = fullfile(output_dir, sprintf('%s_panel_B.svg', figurePrefix));
saveas(fig1, filename, 'svg');
filename = fullfile(output_dir, sprintf('%s_panel_B.png', figurePrefix));
saveas(fig1, filename, 'png');

%%
% PLOT 2: Y-only plot (Membrane proteins) - PANEL C
fig2 = figure('Units', 'centimeters', 'Position', [1, 1, 13, 12.5]);
hold on;

markers = {'o', '^', 's', 'd', '>', 'v', '<'};
colors = {
        sscanf('6B8E23','%2x%2x%2x',[1 3])/255,... % Olive Green
        sscanf('008080','%2x%2x%2x',[1 3])/255,... % Teal
        sscanf('40E0D0','%2x%2x%2x',[1 3])/255,... % Turquoise
        sscanf('87CEEB','%2x%2x%2x',[1 3])/255,... % Sky Blue
        sscanf('6495ED','%2x%2x%2x',[1 3])/255,... % Cornflower Blue
        sscanf('1F56A8','%2x%2x%2x',[1 3])/255,... % blue
        [0.07450980392156863, 0.06666666666666667, 0.0196078431372549],... %black
        sscanf('FFFFFF','%2x%2x%2x',[1 3])/255,... %white
             };
         
xlabel('\phi Y /\phi^{P}_{max}', 'FontSize', 14)
ylabel('Acetate Threshold (\lambda_{ac}) (h^{-1})', 'FontSize', 14)
set(gca, 'FontSize', 12)

% Plot Membrane proteins for specified gamma values
gamma_select_reversed = gamma_select(end:-1:1);

for i=1:numel(gamma_select_reversed)
    gammaVal = gamma_select_reversed(i);
    sub_y_df = y_df(y_df.gamma == gammaVal,:); 
    sub_y_df = sortrows(sub_y_df, {'gamma','Y_divided_by_max'});
    
    if height(sub_y_df) == 0
        continue;
    end
    
    % FIXED: Preserve only the first zero acetate threshold (transition point) for cleaner plots
    % Keep all data up to and including first zero, remove everything after
    indexZero = find(sub_y_df.t_ac == 0, 1);
    if ~isempty(indexZero)
        sub_y_df = sub_y_df(1:indexZero, :);  % Keep up to and including first zero
    else
        % If no zeros, just remove negative values
        validPoints = sub_y_df.t_ac >= 0;
        sub_y_df = sub_y_df(validPoints, :);
    end
    
    if height(sub_y_df) == 0
        continue;
    end
        
    plot(sub_y_df.Y_divided_by_max, sub_y_df.t_ac, '-', 'Color', colors{i}, 'LineWidth', 4, 'HandleVisibility', 'off');
           
    % Downsample the data appropriately
    if i == numel(gamma_select_reversed)
        downsampleFactor = max(1, floor(height(sub_y_df) / 8));
        indicesToShow = 1 : downsampleFactor : height(sub_y_df);
        if height(sub_y_df) > 20
            indicesToShow = sort([indicesToShow, height(sub_y_df)-1, height(sub_y_df)]);
        end
    elseif i == numel(gamma_select_reversed)-1
        downsampleFactor = max(1, floor(height(sub_y_df) / 8));
        indicesToShow = 1 : downsampleFactor : height(sub_y_df);
    else
        downsampleFactor = max(1, floor(height(sub_y_df) / 10));
        indicesToShow = 1 : downsampleFactor : height(sub_y_df);
    end
    
    indicesToShow = unique(indicesToShow);
    indicesToShow = indicesToShow(indicesToShow <= height(sub_y_df));
    
    gammaLabel = sprintf('Membrane (\\phi Y) \\gamma =  %s', num2str(gammaVal)); 
    plot(sub_y_df.Y_divided_by_max(indicesToShow), sub_y_df.t_ac(indicesToShow),...
    markers{i+1}, 'Color', colors{i}, 'LineWidth', 1.5,'MarkerEdgeColor', colors{i},...
    'MarkerFaceColor', colors{end}, 'MarkerSize', 8, 'DisplayName', gammaLabel);
end

xlim([0, 1]);
% Set the legend
lgd = legend('Location', 'southoutside','NumColumns', 2);
set(lgd, 'Box', 'off', 'FontSize', 12);
% Show the grid
grid on;

% Save the figure as SVG
filename = fullfile(output_dir, sprintf('%s_panel_C.svg', figurePrefix));
saveas(fig2, filename, 'svg');
filename = fullfile(output_dir, sprintf('%s_panel_C.png', figurePrefix));
saveas(fig2, filename, 'png');

%%
% PLOT 3: Combined plot (both Z and Y with multiple gamma values) - COMBINED (no panel letter)
fig3 = figure('Units', 'centimeters', 'Position', [1, 1, 13, 12.5]);
hold on;

markers = {'-o', '-^', '-s', '-d', '->', '-v', '-<'};
colors_combined = {
        sscanf('C00217','%2x%2x%2x',[1 3])/255,... % dark red for Z
        sscanf('6B8E23','%2x%2x%2x',[1 3])/255,... % Olive Green for Y gamma=0.2
        sscanf('008080','%2x%2x%2x',[1 3])/255,... % Teal for Y gamma=0.23
        sscanf('40E0D0','%2x%2x%2x',[1 3])/255,... % Turquoise for Y gamma=0.25
        sscanf('87CEEB','%2x%2x%2x',[1 3])/255,... % Sky Blue
        sscanf('6495ED','%2x%2x%2x',[1 3])/255,... % Cornflower Blue
        sscanf('1F56A8','%2x%2x%2x',[1 3])/255,... % blue for Y gamma=1
        [0.07450980392156863, 0.06666666666666667, 0.0196078431372549],... %black
        sscanf('FFFFFF','%2x%2x%2x',[1 3])/255,... %white
             };

xlabel('\phi Z /\phi^{P}_{max} or \phi Y /\phi^{P}_{max}', 'FontSize', 14)
ylabel('Acetate Threshold (\lambda_{ac}) (h^{-1})', 'FontSize', 14)
set(gca, 'FontSize', 12)

% Plot Membrane proteins
gamma_select_plot = [1, 0.35,0.3,0.25, 0.22, 0.2];
gamma_select_reversed = gamma_select_plot(end:-1:1);

for i=1:numel(gamma_select_reversed)
    gammaVal = gamma_select_reversed(i);
    sub_y_df = y_df(y_df.gamma == gammaVal,:); 
    sub_y_df = sortrows(sub_y_df, {'gamma','Y_divided_by_max'});
    
    if height(sub_y_df) == 0
        continue;
    end
    
    % FIXED: Preserve only the first zero acetate threshold (transition point) for cleaner plots
    % Keep all data up to and including first zero, remove everything after
    indexZero = find(sub_y_df.t_ac == 0, 1);
    if ~isempty(indexZero)
        sub_y_df = sub_y_df(1:indexZero, :);  % Keep up to and including first zero
    else
        % If no zeros, just remove negative values
        validPoints = sub_y_df.t_ac >= 0;
        sub_y_df = sub_y_df(validPoints, :);
    end
    
    if height(sub_y_df) == 0
        continue;
    end
        
    plot(sub_y_df.Y_divided_by_max, sub_y_df.t_ac, '-', 'Color', colors_combined{1+i}, 'LineWidth', 4, 'HandleVisibility', 'off');
    
    % Downsample the data to only keep points with markers
    downsampleFactor = max(1, floor(height(sub_y_df) / 10));
    indicesToShow = 1 : downsampleFactor : height(sub_y_df);
       
    if gammaVal == 1 && height(sub_y_df) > 20
        indicesToShow = sort([indicesToShow, height(sub_y_df)-1, height(sub_y_df)]);
    end
    
    indicesToShow = unique(indicesToShow);
    indicesToShow = indicesToShow(indicesToShow <= height(sub_y_df));
    
    gammaLabel = sprintf('Membrane (\\phi Y) \\gamma =  %s', num2str(gammaVal)); 
    plot(sub_y_df.Y_divided_by_max(indicesToShow), sub_y_df.t_ac(indicesToShow),...
    markers{1+i}, 'Color', colors_combined{1+i}, 'LineWidth', 1.5,'MarkerEdgeColor', colors_combined{1+i},...
    'MarkerFaceColor', colors_combined{end}, 'MarkerSize', 8, 'DisplayName', gammaLabel);
end

% Plot Cytosolic Protein (multiple gamma values for better comparison)
gamma_select_z_plot = [1, 0.35, 0.3, 0.25, 0.22, 0.2];  % Show key gamma values for Z (same as Plot 4)
gamma_select_z_reversed = gamma_select_z_plot(end:-1:1);  % Reverse order like Plot 1
z_markers = {'-o', '-s', '-d', '-^', '->', '-v'};  % Different markers for each Z gamma (expanded to match gamma values)

% Define Plot 1 colors for Z data consistency
z_plot1_colors = {
    sscanf('FFA500','%2x%2x%2x',[1 3])/255,... % Mango Orange 
    sscanf('FF6F61','%2x%2x%2x',[1 3])/255,... % Coral 
    sscanf('FF3E3E','%2x%2x%2x',[1 3])/255,... % Tomato Red
    sscanf('DC143C','%2x%2x%2x',[1 3])/255,... % Crimson
    sscanf('C00217','%2x%2x%2x',[1 3])/255,... % dark red
    [0.07450980392156863, 0.06666666666666667, 0.0196078431372549]... %black
};

for i=1:numel(gamma_select_z_reversed)
    gammaVal = gamma_select_z_reversed(i);
    sub_z_df = z_df(z_df.gamma == gammaVal,:);
    
    if height(sub_z_df) == 0
        continue;
    end
    
    % FIXED: Preserve only the first zero acetate threshold (transition point) for cleaner plots
    % Keep all data up to and including first zero, remove everything after
    indexZero = find(sub_z_df.t_ac == 0, 1);
    if ~isempty(indexZero)
        sub_z_df = sub_z_df(1:indexZero, :);  % Keep up to and including first zero
    else
        % If no zeros, just remove negative values
        validPoints = sub_z_df.t_ac >= 0;
        sub_z_df = sub_z_df(validPoints, :);
    end
    
    if height(sub_z_df) == 0
        continue;
    end
    
    % Use same colors as Plot 1 for Z data (now using reversed order like Plot 1)
    z_color = z_plot1_colors{i};
    
    plot(sub_z_df.Z_divided_by_max, sub_z_df.t_ac, '-', 'Color', z_color, 'LineWidth', 4, 'HandleVisibility', 'off');
    
    % Downsample the data to only keep points with markers
    downsampleFactor = max(1, floor(height(sub_z_df) / 8));
    indicesToShow = 1 : downsampleFactor : height(sub_z_df);
    indicesToShow = sort([indicesToShow, height(sub_z_df)-1, height(sub_z_df)]);
    
    indicesToShow = unique(indicesToShow);
    indicesToShow = indicesToShow(indicesToShow <= height(sub_z_df));
    
    gammaLabel = sprintf('Cytosolic (\\phi Z) \\gamma =  %s', num2str(gammaVal));
    plot(sub_z_df.Z_divided_by_max(indicesToShow), sub_z_df.t_ac(indicesToShow),...
        z_markers{i}, 'Color', z_color, 'LineWidth', 1.5,'MarkerEdgeColor', z_color,...
        'MarkerFaceColor', colors_combined{end}, 'MarkerSize', 8, 'DisplayName', gammaLabel);
end

xlim([0, 1]);
% Set the legend with more columns to accommodate additional gamma values
lgd = legend('Location', 'southoutside','NumColumns', 2);
set(lgd, 'Box', 'off', 'FontSize', 10);
% Show the grid
grid on;

% Save the figure as SVG
filename = fullfile(output_dir, sprintf('%s_combined.svg', figurePrefix));
saveas(fig3, filename, 'svg');
filename = fullfile(output_dir, sprintf('%s_combined.png', figurePrefix));
saveas(fig3, filename, 'png');

%%
% PLOT 4: Combined plot with Y data normalized by membrane constraint (phi_Mmax) - PANEL D
fig4 = figure('Units', 'centimeters', 'Position', [1, 1, 13, 12.5]);
hold on;

markers = {'-o', '-^', '-s', '-d', '->', '-v', '-<'};
colors_combined = {
        sscanf('C00217','%2x%2x%2x',[1 3])/255,... % dark red for Z
        sscanf('6B8E23','%2x%2x%2x',[1 3])/255,... % Olive Green for Y gamma=0.2
        sscanf('008080','%2x%2x%2x',[1 3])/255,... % Teal for Y gamma=0.23
        sscanf('40E0D0','%2x%2x%2x',[1 3])/255,... % Turquoise for Y gamma=0.25
        sscanf('87CEEB','%2x%2x%2x',[1 3])/255,... % Sky Blue
        sscanf('6495ED','%2x%2x%2x',[1 3])/255,... % Cornflower Blue
        sscanf('1F56A8','%2x%2x%2x',[1 3])/255,... % blue for Y gamma=1
        [0.07450980392156863, 0.06666666666666667, 0.0196078431372549],... %black
        sscanf('FFFFFF','%2x%2x%2x',[1 3])/255,... %white
             };

% xlabel('\phi Z /\phi^{P}_{max} or \phi Y /\phi^{M}_{max}', 'FontSize', 14)
xlabel('\phi Y /\phi^{M}_{max}', 'FontSize', 14)

ylabel('Acetate Threshold (\lambda_{ac}) (h^{-1})', 'FontSize', 14)
set(gca, 'FontSize', 12)

% Plot Membrane proteins using Y_divided_by_Mmax normalization
gamma_select_plot = [1, 0.35,0.3,0.25, 0.22, 0.2];
gamma_select_reversed = gamma_select_plot(end:-1:1);

for i=1:numel(gamma_select_reversed)
    gammaVal = gamma_select_reversed(i);
    sub_y_df = y_df(y_df.gamma == gammaVal,:); 
    sub_y_df = sortrows(sub_y_df, {'gamma','Y_divided_by_Mmax'});  % Sort by Mmax normalization
    
    if height(sub_y_df) == 0
        continue;
    end
    
    % FIXED: Preserve only the first zero acetate threshold (transition point) for cleaner plots
    % Keep all data up to and including first zero, remove everything after
    indexZero = find(sub_y_df.t_ac == 0, 1);
    if ~isempty(indexZero)
        sub_y_df = sub_y_df(1:indexZero, :);  % Keep up to and including first zero
    else
        % If no zeros, just remove negative values
        validPoints = sub_y_df.t_ac >= 0;
        sub_y_df = sub_y_df(validPoints, :);
    end
    
    if height(sub_y_df) == 0
        continue;
    end
        
    plot(sub_y_df.Y_divided_by_Mmax, sub_y_df.t_ac, '-', 'Color', colors_combined{1+i}, 'LineWidth', 4, 'HandleVisibility', 'off');
    
    % Downsample the data to only keep points with markers
    downsampleFactor = max(1, floor(height(sub_y_df) / 10));
    indicesToShow = 1 : downsampleFactor : height(sub_y_df);
       
    if gammaVal == 1 && height(sub_y_df) > 20
        indicesToShow = sort([indicesToShow, height(sub_y_df)-1, height(sub_y_df)]);
    end
    
    indicesToShow = unique(indicesToShow);
    indicesToShow = indicesToShow(indicesToShow <= height(sub_y_df));
    
    gammaLabel = sprintf('Membrane (\\phi Y) \\gamma =  %s', num2str(gammaVal)); 
    plot(sub_y_df.Y_divided_by_Mmax(indicesToShow), sub_y_df.t_ac(indicesToShow),...
    markers{1+i}, 'Color', colors_combined{1+i}, 'LineWidth', 1.5,'MarkerEdgeColor', colors_combined{1+i},...
    'MarkerFaceColor', colors_combined{end}, 'MarkerSize', 8, 'DisplayName', gammaLabel);
end

% % Plot Cytosolic Protein (same as Plot 3 - normalized by phi_max)
% gamma_select_z_plot = [1, 0.35, 0.3, 0.25, 0.22, 0.2];  % Show key gamma values for Z
% gamma_select_z_reversed = gamma_select_z_plot(end:-1:1);  % Reverse order like Plot 1
% z_markers = {'-o', '-s', '-d', '-^', '->', '-v'};  % Different markers for each Z gamma (expanded to match gamma values)
% 
% % Define Plot 1 colors for Z data consistency (same as Plot 3)
% z_plot1_colors = {
%     sscanf('FFA500','%2x%2x%2x',[1 3])/255,... % Mango Orange 
%     sscanf('FF6F61','%2x%2x%2x',[1 3])/255,... % Coral 
%     sscanf('FF3E3E','%2x%2x%2x',[1 3])/255,... % Tomato Red
%     sscanf('DC143C','%2x%2x%2x',[1 3])/255,... % Crimson
%     sscanf('C00217','%2x%2x%2x',[1 3])/255,... % dark red
%     [0.07450980392156863, 0.06666666666666667, 0.0196078431372549]... %black
% };
% 
% for i=1:numel(gamma_select_z_reversed)
%     gammaVal = gamma_select_z_reversed(i);
%     sub_z_df = z_df(z_df.gamma == gammaVal,:);
% 
%     if height(sub_z_df) == 0
%         continue;
%     end
% 
%     % FIXED: Preserve only the first zero acetate threshold (transition point) for cleaner plots
%     % Keep all data up to and including first zero, remove everything after
%     indexZero = find(sub_z_df.t_ac == 0, 1);
%     if ~isempty(indexZero)
%         sub_z_df = sub_z_df(1:indexZero, :);  % Keep up to and including first zero
%     else
%         % If no zeros, just remove negative values
%         validPoints = sub_z_df.t_ac >= 0;
%         sub_z_df = sub_z_df(validPoints, :);
%     end
% 
%     if height(sub_z_df) == 0
%         continue;
%     end
% 
%     % Use same colors as Plot 1 for Z data (now using reversed order like Plot 1)
%     z_color = z_plot1_colors{i};
% 
%     plot(sub_z_df.Z_divided_by_max, sub_z_df.t_ac, '-', 'Color', z_color, 'LineWidth', 4, 'HandleVisibility', 'off');
% 
%     % Downsample the data to only keep points with markers
%     downsampleFactor = max(1, floor(height(sub_z_df) / 8));
%     indicesToShow = 1 : downsampleFactor : height(sub_z_df);
%     indicesToShow = sort([indicesToShow, height(sub_z_df)-1, height(sub_z_df)]);
% 
%     indicesToShow = unique(indicesToShow);
%     indicesToShow = indicesToShow(indicesToShow <= height(sub_z_df));
% 
%     gammaLabel = sprintf('Cytosolic (\\phi Z) \\gamma =  %s', num2str(gammaVal));
%     plot(sub_z_df.Z_divided_by_max(indicesToShow), sub_z_df.t_ac(indicesToShow),...
%         z_markers{i}, 'Color', z_color, 'LineWidth', 1.5,'MarkerEdgeColor', z_color,...
%         'MarkerFaceColor', colors_combined{end}, 'MarkerSize', 8, 'DisplayName', gammaLabel);
% end

xlim([0, 1]);
% Set the legend with more columns to accommodate additional gamma values
lgd = legend('Location', 'southoutside','NumColumns', 2);
set(lgd, 'Box', 'off', 'FontSize', 10);
% Show the grid
grid on;

% Save the figure as SVG
filename = fullfile(output_dir, sprintf('%s_panel_D.svg', figurePrefix));
saveas(fig4, filename, 'svg');
filename = fullfile(output_dir, sprintf('%s_panel_D.png', figurePrefix));
saveas(fig4, filename, 'png');

fprintf('\n=== FIGURE 8 COMPLETE ===\n');
fprintf('Generated normalized protein expression figures:\n');
fprintf('  Panel B (Membrane Z, phi_max^P): %s\n', fullfile(output_dir, sprintf('%s_panel_B.svg', figurePrefix)));
fprintf('  Panel C (Cytosolic Y, phi_max^P): %s\n', fullfile(output_dir, sprintf('%s_panel_C.svg', figurePrefix)));
fprintf('  Panel D (Cytosolic Y, phi_max^M): %s\n', fullfile(output_dir, sprintf('%s_panel_D.svg', figurePrefix)));
fprintf('  Combined (Y+Z): %s\n', fullfile(output_dir, sprintf('%s_combined.svg', figurePrefix)));
fprintf('\nTotal: 4 figures × 2 formats (SVG + PNG) = 8 files\n');
fprintf('Output location: %s\n', output_dir);

%%
function gr_ac = get_acetateT(T)
    % Advanced acetate threshold detection method
    % This method uses a more stringent initial slope threshold for better accuracy
    
    % Access the global validation flag
    global ENABLE_VALIDATION_WARNINGS;
    if isempty(ENABLE_VALIDATION_WARNINGS)
        ENABLE_VALIDATION_WARNINGS = false;
    end
    
    % Ensure the table is sorted by growth rate
    T = sortrows(T, 'gr');
    
    % Calculate the slope of membrane/respiration flux vs growth rate
    dif_v = diff(T.memb);  % Change in membrane flux (respiration)
    dif_g = diff(T.gr);    % Change in growth rate
    
    % Handle edge case where growth rate doesn't change
    if any(dif_g == 0)
        if ENABLE_VALIDATION_WARNINGS
            fprintf('Warning: Zero growth rate differences detected. Using fallback method.\n');
        end
        % Find first point with meaningful acetate production
        ac_idx = find(T.ac >= 0.01, 1, 'first');
        if isempty(ac_idx)
            ac_idx = find(T.ac >= 0.005, 1, 'first');
        end
        if isempty(ac_idx)
            ac_idx = 1;
        end
        gr_ac = T.gr(ac_idx);
        return;
    end
    
    slope_v = dif_v./dif_g;
    m_s = max(slope_v);

    % Track which method was used for debugging
    method_used = 'none';
    
    % Advanced method: Start with much more stringent threshold (0.001 vs 0.01)
    % This detects the transition point more accurately
    memb_s = find(slope_v <= m_s*0.001);
    t_i = find(T.ac(memb_s) >= 0.01, 1, 'first');
    ac_t = memb_s(t_i);
    if ~isempty(ac_t)
        method_used = 'stringent_0.001';
    end

    % First fallback: Use 5% threshold
    if isempty(ac_t)
        memb_s = find(slope_v <= m_s*0.05);
        t_i = find(T.ac(memb_s) >= 0, 1, 'first');
        ac_t = memb_s(t_i);
        if ~isempty(ac_t)
            method_used = 'fallback_0.05';
        end
    end
    
    % Second fallback: Use 20% threshold  
    if isempty(ac_t)
        memb_s = find(slope_v <= m_s*0.2);
        t_i = find(T.ac(memb_s) >= 0, 1, 'first');
        ac_t = memb_s(t_i);
        if ~isempty(ac_t)
            method_used = 'fallback_0.2';
        end
    end

    % Final fallback: Use first data point if no threshold found
    if isempty(ac_t)
        ac_t = 1;
        method_used = 'final_fallback';
        if ENABLE_VALIDATION_WARNINGS
            fprintf('Warning: No acetate threshold found, using first data point.\n');
        end
    end
    
    % Return the growth rate at the detected acetate threshold
    gr_ac = T.gr(ac_t);
    
    % Add validation to check for abrupt changes and data quality
    if length(T.gr) > 2
        % Check if the detected threshold makes sense relative to the data
        growth_range = max(T.gr) - min(T.gr);
        
        % Validate that the threshold isn't too early (within first 10% of growth range)
        if gr_ac < min(T.gr) + 0.1 * growth_range && strcmp(method_used, 'stringent_0.001')
            if ENABLE_VALIDATION_WARNINGS
                fprintf('Info: Early acetate threshold detected with stringent method (gr=%.4f, method=%s)\n', gr_ac, method_used);
            end
            
            % Additional validation: check for consistent acetate production
            consistent_ac_idx = find(T.ac >= 0.005);
            if length(consistent_ac_idx) >= 2
                % Ensure we're not picking a spurious early point
                ac_values_at_threshold = T.ac(max(1, ac_t-1):min(length(T.ac), ac_t+1));
                if any(ac_values_at_threshold < 0.001) && ac_t > 1
                    % Look for a more stable threshold
                    stable_idx = find(T.ac >= 0.005 & [0; diff(T.ac)] >= 0, 1, 'first');
                    if ~isempty(stable_idx) && stable_idx > ac_t
                        ac_t = stable_idx;
                        gr_ac = T.gr(ac_t);
                        if ENABLE_VALIDATION_WARNINGS
                            fprintf('Adjusted to more stable threshold at gr=%.4f\n', gr_ac);
                        end
                    end
                end
            end
        end
        
        % Check for potential data quality issues
        if length(unique(T.ac)) < 3 && ENABLE_VALIDATION_WARNINGS
            fprintf('Warning: Low acetate variability detected - results may be unreliable.\n');
        end
    end
    
    % Optional debug output
    if ENABLE_VALIDATION_WARNINGS && ~strcmp(method_used, 'stringent_0.001')
        fprintf('Acetate threshold detection used method: %s (gr=%.4f)\n', method_used, gr_ac);
    end
end
