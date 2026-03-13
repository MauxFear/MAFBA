% LOAD DATA FROM CSV

% Set file path for the CSV file
baseDir = fileparts(mfilename('fullpath'));
repoRoot = fullfile(baseDir, '..', '..');
datafolder = fullfile(repoRoot, 'data', 'outputs', 'iML1515_Homogeneous');
fprintf('Data folder: %s\n', datafolder);

% Set the data file name for homogeneous analysis
figTag = 'Homogeneous';
datafilename = fullfile(datafolder, 'dataFlux_Homogeneous.csv');

% If the specific file doesn't exist, try to find any CSV file in the directory
if ~isfile(datafilename)
    fprintf('Specified data file not found. Looking for matching CSV files...\n');
    
    % Get a list of all CSV files in the directory
    files = dir(fullfile(datafolder, 'dataFlux_*.csv'));
    
    if isempty(files)
        error('No CSV files found in %s. Please check the directory.', datafolder);
    else
        % Use the newest CSV file found
        [~, newestIdx] = max([files.datenum]);
        datafilename = fullfile(datafolder, files(newestIdx).name);
        fprintf('Using file: %s\n', files(newestIdx).name);
    end
end

fprintf('Data file found: %s\n', datafilename);

% Read the table data from CSV file
tableData = readtable(datafilename);

% Display column names for debugging
fprintf('Data columns: %s\n', strjoin(tableData.Properties.VariableNames, ', '));

% Check for expected column names and adapt if needed
if ismember('Gamma', tableData.Properties.VariableNames)
    gamma_col = 'Gamma';
elseif ismember('gamaVal', tableData.Properties.VariableNames)
    gamma_col = 'gamaVal';
else
    % If no gamma column found, assume it's 1
    fprintf('No gamma column found. Assuming gamma = 1\n');
    tableData.Gamma = ones(height(tableData), 1);
    gamma_col = 'Gamma';
end

if ismember('Change_Factor', tableData.Properties.VariableNames)
    cf_col = 'Change_Factor';
elseif ismember('changeFactor', tableData.Properties.VariableNames)
    cf_col = 'changeFactor';
else
    % If no change factor column found, assume it's 1
    fprintf('No change factor column found. Assuming change factor = 1\n');
    tableData.Change_Factor = ones(height(tableData), 1);
    cf_col = 'Change_Factor';
end

% Map reaction names to expected column names
reaction_map = struct();
if ismember('rx_EX_glc__D_e_b', tableData.Properties.VariableNames)
    reaction_map.glc = 'rx_EX_glc__D_e_b';
elseif ismember('EX_glc__D_e_b', tableData.Properties.VariableNames)
    reaction_map.glc = 'EX_glc__D_e_b';
elseif ismember('glc', tableData.Properties.VariableNames)
    reaction_map.glc = 'glc';
else
    error('Glucose uptake column not found');
end

if ismember('rx_EX_ac_e_f', tableData.Properties.VariableNames)
    reaction_map.ac = 'rx_EX_ac_e_f';
elseif ismember('EX_ac_e_f', tableData.Properties.VariableNames)
    reaction_map.ac = 'EX_ac_e_f';
elseif ismember('ac', tableData.Properties.VariableNames)
    reaction_map.ac = 'ac';
else
    error('Acetate excretion column not found');
end

if ismember('rx_AKGDH', tableData.Properties.VariableNames)
    reaction_map.akgdh = 'rx_AKGDH';
elseif ismember('AKGDH', tableData.Properties.VariableNames)
    reaction_map.akgdh = 'AKGDH';
elseif ismember('akgdh', tableData.Properties.VariableNames)
    reaction_map.akgdh = 'akgdh';
else
    error('AKGDH column not found');
end

if ismember('rx_ATPS4rpp_f', tableData.Properties.VariableNames)
    reaction_map.memb = 'rx_ATPS4rpp_f';
elseif ismember('ATPS4rpp_f', tableData.Properties.VariableNames)
    reaction_map.memb = 'ATPS4rpp_f';
elseif ismember('memb', tableData.Properties.VariableNames)
    reaction_map.memb = 'memb';
else
    error('Membrane/ATPS4rpp column not found');
end

if ismember('rx_MALS', tableData.Properties.VariableNames)
    reaction_map.mals = 'rx_MALS';
elseif ismember('MALS', tableData.Properties.VariableNames)
    reaction_map.mals = 'MALS';
elseif ismember('mals', tableData.Properties.VariableNames)
    reaction_map.mals = 'mals';
else
    error('MALS column not found');
end

if ismember('rx_EDD', tableData.Properties.VariableNames)
    reaction_map.edd = 'rx_EDD';
elseif ismember('EDD', tableData.Properties.VariableNames)
    reaction_map.edd = 'EDD';
elseif ismember('edd', tableData.Properties.VariableNames)
    reaction_map.edd = 'edd';
else
    error('EDD column not found');
end

if ismember('rx_EX_co2_e_f', tableData.Properties.VariableNames)
    reaction_map.co2 = 'rx_EX_co2_e_f';
elseif ismember('EX_co2_e_f', tableData.Properties.VariableNames)
    reaction_map.co2 = 'EX_co2_e_f';
elseif ismember('co2', tableData.Properties.VariableNames)
    reaction_map.co2 = 'co2';
else
    error('CO2 excretion column not found');
end

% Find growth rate column
if ismember('rx_BIOMASS_Ec_iML1515_WT_75p37M', tableData.Properties.VariableNames)
    reaction_map.gr = 'rx_BIOMASS_Ec_iML1515_WT_75p37M';
elseif ismember('BIOMASS_Ec_iML1515_WT_75p37M', tableData.Properties.VariableNames)
    reaction_map.gr = 'BIOMASS_Ec_iML1515_WT_75p37M';
elseif ismember('rx_BIOMASS_Ec_iML1515_core_75p37M', tableData.Properties.VariableNames)
    reaction_map.gr = 'rx_BIOMASS_Ec_iML1515_core_75p37M';
elseif ismember('BIOMASS_Ec_iML1515_core_75p37M', tableData.Properties.VariableNames)
    reaction_map.gr = 'BIOMASS_Ec_iML1515_core_75p37M';
elseif ismember('gr', tableData.Properties.VariableNames)
    reaction_map.gr = 'gr';
else
    error('Growth rate column not found');
end

% Get unique gamma values
gamma_values = unique(tableData.(gamma_col));
fprintf('Found %d unique gamma values\n', length(gamma_values));

% Get unique change factors
change_factors = unique(tableData.(cf_col));
fprintf('Found %d unique change factors\n', length(change_factors));

% Create a timestamped output folder
output_folder = fullfile(baseDir, 'output');

% Create the folder if it doesn't exist
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
    fprintf('Created output folder: %s\n', output_folder);
end

% Display information about the analysis
fprintf('Creating plots for %d gamma values and %d change factors\n', length(gamma_values), length(change_factors));
fprintf('Total plots to create: %d\n', length(gamma_values) * length(change_factors));

% Loop through each gamma value
for g = 1:length(gamma_values)
    gamma_val = gamma_values(g);
    gamma_string = ['\gamma = ', num2str(gamma_val)];
    fprintf('Processing gamma = %f (%d of %d)\n', gamma_val, g, length(gamma_values));
    
    % Loop through each change factor
    for c = 1:length(change_factors)
        change_factor = change_factors(c);
        fprintf('  Processing change factor = %f (%d of %d)\n', change_factor, c, length(change_factors));
        
        % Filter the table data for the current gamma value and change factor
        T = tableData(tableData.(gamma_col) == gamma_val & tableData.(cf_col) == change_factor, :);
        
        % Skip if no data for this combination
        if height(T) < 2
            fprintf('  Skipping: Not enough data points for gamma=%f, change factor=%f\n', gamma_val, change_factor);
            continue;
        end
        
        % Sort rows by growth rate
        T = sortrows(T, reaction_map.gr);
        
        % Detect acetate threshold based on acetate production and respiration slope change
        try
            % Step 1: First check if acetate production exists
            ac_positive_indices = find(T.(reaction_map.ac) >= 0.01);
            
            if isempty(ac_positive_indices)
                % No acetate production above threshold, use first point
                fprintf('  No acetate production above 0.01 found. Using first point.\n');
                ac_t = 1;
            else
                % Step 2: Calculate respiration slope changes
                dif_memb = diff(T.(reaction_map.memb));
                dif_gr = diff(T.(reaction_map.gr));
                
                % Calculate slopes with respect to growth rate
                slope_memb = dif_memb ./ dif_gr;
                
                % Handle case where slopes might be empty or all NaN
                if isempty(slope_memb) || all(isnan(slope_memb)) || all(isinf(slope_memb))
                    % Can't calculate slopes, use first acetate point
                    ac_t = ac_positive_indices(1);
                    fprintf('  Cannot calculate valid respiration slopes. Using first acetate point at growth rate %.4f\n', T.(reaction_map.gr)(ac_t));
                else
                    % Find the maximum respiration slope
                    valid_slopes = slope_memb(~isnan(slope_memb) & ~isinf(slope_memb));
                    if isempty(valid_slopes)
                        max_slope = 0;
                    else
                        max_slope = max(valid_slopes);
                    end
                    
                    % Find where respiration slope significantly decreases
                    slope_threshold = 0.1; % Threshold for significant slope change
                    slope_drop_indices = find(slope_memb <= max_slope * slope_threshold);
                    
                    if isempty(slope_drop_indices)
                        % No significant slope change, use first acetate point
                        ac_t = ac_positive_indices(1);
                        fprintf('  No significant respiration slope change. Using first acetate point at growth rate %.4f\n', T.(reaction_map.gr)(ac_t));
                    else
                        % Find the first point where both conditions are met:
                        % 1. Respiration slope has dropped
                        % 2. Acetate is above threshold
                        
                        % Find the first slope drop index that occurs at or before the first acetate point
                        first_ac_idx = ac_positive_indices(1);
                        valid_slope_drops = slope_drop_indices(slope_drop_indices <= first_ac_idx);
                        
                        if ~isempty(valid_slope_drops)
                            % Use the last slope drop before acetate appears
                            ac_t = valid_slope_drops(end);
                            fprintf('  Found threshold where respiration slope changes before acetate production at growth rate %.4f\n', T.(reaction_map.gr)(ac_t));
                        else
                            % No slope drop before acetate, use first acetate point
                            ac_t = first_ac_idx;
                            fprintf('  No respiration slope change before acetate. Using first acetate point at growth rate %.4f\n', T.(reaction_map.gr)(ac_t));
                        end
                    end
                end
            end
        catch e
            % If any error occurs, use a safe default
            fprintf('  Error detecting acetate threshold: %s\n', e.message);
            fprintf('  Using default threshold.\n');
            ac_t = 1;
        end
        
        % Ensure the index is valid
        if isempty(ac_t) || ac_t < 1
            ac_t = 1;
        end
        
        if ac_t > height(T)
            ac_t = height(T);
        end
        
        % Get the growth rate at the acetate threshold
        gr_ac = T.(reaction_map.gr)(ac_t);
        fprintf('  Acetate threshold detected at growth rate = %.4f\n', gr_ac);
        
        % Plot the data
        fig = get_plot(T, [gamma_string, ', CF = ', num2str(change_factor)], gr_ac, reaction_map);
        
        if numel(change_factors) > 1
            filename = fullfile(output_folder, sprintf('%s_flux_g%03d_cf%03d.svg', figTag, round(gamma_val * 100), round(change_factor * 100)));
            filename_png = fullfile(output_folder, sprintf('%s_flux_g%03d_cf%03d.png', figTag, round(gamma_val * 100), round(change_factor * 100)));
        else
            filename = fullfile(output_folder, sprintf('%s_flux_g%03d.svg', figTag, round(gamma_val * 100)));
            filename_png = fullfile(output_folder, sprintf('%s_flux_g%03d.png', figTag, round(gamma_val * 100)));
        end

        % Save the figure
        saveas(fig, filename, 'svg');
        saveas(fig, filename_png, 'png');
        fprintf('  Saved flux plot for gamma=%f, change factor=%f.\n', gamma_val, change_factor);
        
        % Close the figure to free up memory
        close(fig);
    end
end

fprintf('All plots have been saved to: %s\n', output_folder);


function fig = get_plot(T, gamma_string, gr_ac, reaction_map)
    % Create a figure with the specified size and position in centimeters
    fig = figure('Units', 'centimeters', 'Position', [1, 1, 20, 13]);
    hold on;
    
    % Define colors for each reaction - using the same colors as in gradientPlot
    color = {sscanf('65AF65','%2x%2x%2x',[1 3])/255, ... % green - Glucose uptake
             sscanf('C74848','%2x%2x%2x',[1 3])/255, ... % red - Acetate excretion
             sscanf('9065D1','%2x%2x%2x',[1 3])/255, ... % purple - TCA flux
             sscanf('266DD7','%2x%2x%2x',[1 3])/255, ... % blue - Respiration
             sscanf('F2A1E4','%2x%2x%2x',[1 3])/255, ... % pink - Glyoxylate shunt
             sscanf('D2C725','%2x%2x%2x',[1 3])/255, ... % yellow - ED pathway
             sscanf('B2B2B2','%2x%2x%2x',[1 3])/255, ... % gray - CO2 excretion
             sscanf('E40303','%2x%2x%2x',[1 3])/255, ... % bright red
             [0.07450980392156863, 0.06666666666666667, 0.0196078431372549], ... % black
             sscanf('FFFFFF','%2x%2x%2x',[1 3])/255, ... % white
             };
                 
    title(['Sensitivity Analysis: ', gamma_string], 'FontSize', 14)
    xlabel('Growth rate (h^{-1})', 'FontSize', 12)
    ylabel('Flux (mmol/g_{DW}h)', 'FontSize', 12)

    % Set font size
    set(gca, 'FontSize', 10)
    
    % Get growth rate values
    growth_rate = T.(reaction_map.gr);
    
    % Downsample the data to only keep every 5th point with markers
    downsampleFactor = max(1, floor(height(T) / 20));  % Aim for about 20 points
    numPoints = numel(growth_rate);
    indicesToShow = 1 : downsampleFactor : numPoints;
    
    % Plot each reaction with both lines and markers
    % First plot black lines for all reactions (for visual effect)
    plot(growth_rate, T.(reaction_map.edd), '-', 'Color', 'k', 'LineWidth', 4, 'HandleVisibility', 'off');
    plot(growth_rate, T.(reaction_map.akgdh), '-', 'Color', 'k', 'LineWidth', 4, 'HandleVisibility', 'off');
    plot(growth_rate, T.(reaction_map.mals), '-', 'Color', 'k', 'LineWidth', 4, 'HandleVisibility', 'off');
    plot(growth_rate, T.(reaction_map.co2), '-', 'Color', 'k', 'LineWidth', 4, 'HandleVisibility', 'off');
    plot(growth_rate, T.(reaction_map.glc), '-', 'Color', 'k', 'LineWidth', 4, 'HandleVisibility', 'off');
    plot(growth_rate, 0.5*T.(reaction_map.memb), '-', 'Color', 'k', 'LineWidth', 4, 'HandleVisibility', 'off');
    plot(growth_rate, T.(reaction_map.ac), '-', 'Color', 'k', 'LineWidth', 4, 'HandleVisibility', 'off');
    
    % Then plot colored lines
    plot(growth_rate, T.(reaction_map.edd), '-', 'Color', color{6}, 'LineWidth', 4, 'HandleVisibility', 'off');
    plot(growth_rate, T.(reaction_map.akgdh), '-', 'Color', color{3}, 'LineWidth', 4, 'HandleVisibility', 'off');
    plot(growth_rate, T.(reaction_map.mals), '-', 'Color', color{5}, 'LineWidth', 4, 'HandleVisibility', 'off');
    plot(growth_rate, T.(reaction_map.co2), '-', 'Color', color{7}, 'LineWidth', 4, 'HandleVisibility', 'off');
    plot(growth_rate, T.(reaction_map.glc), '-', 'Color', color{1}, 'LineWidth', 4, 'HandleVisibility', 'off');
    plot(growth_rate, 0.5*T.(reaction_map.memb), '-', 'Color', color{4}, 'LineWidth', 4, 'HandleVisibility', 'off');
    plot(growth_rate, T.(reaction_map.ac), '-', 'Color', color{2}, 'LineWidth', 4, 'HandleVisibility', 'off');
    
    % Set x-axis limits
    if max(growth_rate) > 0
        xlim([0, max(growth_rate) * 1.05]);
    else
        xlim([0, 1]);
    end
    
    % Get y-axis limits
    y = ylim;
    
    % Plot acetate threshold line
    h8 = plot([gr_ac, gr_ac], [0, y(2)], '--', 'Color', color{8}, 'LineWidth', 2, 'MarkerFaceColor', color{2}, 'DisplayName', '\lambda_{ac}');
    ylim(y);
    
    % Plot markers on top of lines for better visibility
    h5 = plot(growth_rate(indicesToShow), T.(reaction_map.edd)(indicesToShow), 's', 'Color', color{6}, 'LineWidth', 1.5, ...
        'MarkerEdgeColor', color{6}, 'MarkerFaceColor', color{10}, 'MarkerSize', 10, 'DisplayName', 'ED pathway');
    h4 = plot(growth_rate(indicesToShow), T.(reaction_map.akgdh)(indicesToShow), '^', 'Color', color{3}, 'LineWidth', 0.5, ...
        'MarkerEdgeColor', color{9}, 'MarkerFaceColor', color{3}, 'MarkerSize', 10, 'DisplayName', 'TCA flux');
    h6 = plot(growth_rate(indicesToShow), T.(reaction_map.mals)(indicesToShow), '^', 'Color', color{5}, 'LineWidth', 1.5, ...
        'MarkerEdgeColor', color{5}, 'MarkerFaceColor', color{10}, 'MarkerSize', 10, 'DisplayName', 'Glyoxylate shunt');
    h7 = plot(growth_rate(indicesToShow), T.(reaction_map.co2)(indicesToShow), 'd', 'Color', color{7}, 'LineWidth', 1.5, ...
        'MarkerEdgeColor', color{7}, 'MarkerFaceColor', color{10}, 'MarkerSize', 10, 'DisplayName', 'CO_2 excretion');
    h1 = plot(growth_rate(indicesToShow), T.(reaction_map.glc)(indicesToShow), 'o', 'Color', color{1}, 'LineWidth', 0.5, ...
        'MarkerEdgeColor', color{9}, 'MarkerFaceColor', color{1}, 'MarkerSize', 10, 'DisplayName', 'Glucose uptake');
    h2 = plot(growth_rate(indicesToShow), T.(reaction_map.memb)(indicesToShow), 'd', 'Color', color{4}, 'LineWidth', 0.5, ...
        'MarkerEdgeColor', color{9}, 'MarkerFaceColor', color{4}, 'MarkerSize', 10, 'DisplayName', 'Respiration 0.5x');    
    h3 = plot(growth_rate(indicesToShow), T.(reaction_map.ac)(indicesToShow), 's', 'Color', color{2}, 'LineWidth', 0.5, ...
        'MarkerEdgeColor', color{9}, 'MarkerFaceColor', color{2}, 'MarkerSize', 10, 'DisplayName', 'Acetate excretion');
    
    % Set the legend
    lgd = legend([h1, h2, h3, h4, h5, h6, h7, h8], ...
        'Location', 'southoutside', 'NumColumns', 4);
    set(lgd, 'Box', 'off', 'FontSize', 12);
    title(lgd, 'Reaction/Pathway', 'FontSize', 12); % Set the legend title
    
    % Show the grid
    grid on;
    
    % Add annotation with the value of gr_ac
    text(gr_ac + 0.02, y(2) * 0.9, sprintf('\\lambda_{ac} = %.2f', gr_ac), 'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'top', 'FontSize', 12, 'FontWeight', 'bold');

    hold off;
end
