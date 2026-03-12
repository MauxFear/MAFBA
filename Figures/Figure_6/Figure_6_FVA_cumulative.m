%% FIGURE 6 - FLUX VARIABILITY ANALYSIS (FVA): CUMULATIVE DISTRIBUTIONS
% Shows cumulative distributions of flux variability comparing:
%   - MAFBA (parametrized model) at multiple gamma values
%   - irreversible iML1515 (no enzymatic constraints)
%
% Data source: FVA analysis from combined_FVA_results
% Outputs: Three cumulative distribution plots (all, membrane, non-membrane reactions)
%% CONFIGURATION
% Set file paths
baseDir = fileparts(mfilename('fullpath'));
repoRoot = fullfile(baseDir, '..', '..');

% Data folder
datafolder = fullfile(repoRoot, 'data', 'outputs', 'FVA_Results');
fprintf('Data folder: %s\n', datafolder);

% Output folder
output_folder = fullfile(baseDir, 'output');
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end
fprintf('Output folder: %s\n', output_folder);

%% SETUP
% Output naming
figurePrefix = 'iML1515_MAFBA_FVA';

% Models to compare
selected_models = {'MAFBA', 'irriML1515'};  % MAFBA = parametrized model with both constraints
models_labels = {'MAFBA', 'irreversible iML1515'};

% Analysis parameters
uptake_bound = 1000;       % High uptake bound
FVA_percentage = 0;        % 0% FVA (exact optimum)
analysis_types = {'all_reactions', 'membrane', 'nonmembrane'};  % Three panels

% Define colors for different gamma values and models
base_colors = {
    sscanf('40E0D0','%2x%2x%2x',[1 3])/255,... % Turquoise (γ=0.20)
    sscanf('FFA500','%2x%2x%2x',[1 3])/255,... % Mango Orange (γ=0.25)
    sscanf('ba81ee','%2x%2x%2x',[1 3])/255, ... % Violet (γ=0.30)
    sscanf('1F56A8','%2x%2x%2x',[1 3])/255,... % Blue (γ=1.00)
    sscanf('6B8E23','%2x%2x%2x',[1 3])/255 ... % Olive Green (FBA only)
};

% Labels for each condition
combined_labels = {
    'Both constraints, \gamma = 0.20',
    'Both constraints, \gamma = 0.25',
    'Both constraints, \gamma = 0.30',
    'No Membrane constraint, \gamma = 1.00',
    'No Enzymatic constraints, FBA only'
};

% Panel titles
panel_titles = {
    'Cumulative Distribution of Flux Variability',
    'Flux Variability of Membrane Reactions',
    'Flux Variability of Non-Membrane Reactions'
};

%% LOAD DATA
fprintf('\n=== LOADING FVA DATA ===\n');

% Get list of FVA result files
fva_files = dir(fullfile(datafolder, 'FVA_results_*.csv'));
if isempty(fva_files)
    error('No FVA result files found in %s', datafolder);
end

% Initialize combined data
all_data = table();

% Process each FVA file
for i = 1:length(fva_files)
    file_path = fullfile(datafolder, fva_files(i).name);
    fprintf('Processing file: %s\n', file_path);
    
    % Read the data
    opts = detectImportOptions(file_path);
    opts.VariableNamesLine = 1;
    opts.DataLines = 2;
    data = readtable(file_path, opts);
    
    % Extract model name from filename
    % Handles: FVA_results_ModelName.csv or FVA_results_v2_ModelName_timestamp.csv
    [~, filename] = fileparts(fva_files(i).name);
    
    % Try pattern 1: FVA_results_ModelName.csv
    model_name = regexp(filename, '^FVA_results_([^_\.]+)$', 'tokens');
    if isempty(model_name)
        % Try pattern 2: FVA_results_v2_ModelName_timestamp.csv
        model_name = regexp(filename, 'FVA_results_(?:v\d+_)?([^_]+)_', 'tokens');
    end
    
    if ~isempty(model_name)
        model_name = model_name{1}{1};
    else
        warning('Could not extract model name from %s', filename);
        continue;
    end
    
    fprintf('  Model: %s\n', model_name);
    
    % Check if this is one of our selected models
    if ~ismember(model_name, selected_models)
        fprintf('  Skipping (not in selected models)\n');
        continue;
    end
    
    % Check for required columns
    if ~ismember('Analysis_Type', data.Properties.VariableNames)
        warning('Analysis_Type column not found in %s', filename);
        continue;
    end
    
    % Filter based on model type
    switch model_name
        case 'MAFBA'
            % Use multiple gamma values for MAFBA
            valid_data = data(...
                strcmp(data.Model, model_name) & ...
                ismember(data.Gamma_fMSA, [1, 0.25, 0.3, 0.2]) & ...
                data.Uptake_Bound == uptake_bound & ...
                data.FVA_Percent == FVA_percentage & ...
                ismember(data.Analysis_Type, analysis_types), :);
        case 'irriML1515'
            % Use only gamma_fMSA = 1 for irreversible iML1515 (no enzymatic constraints)
            valid_data = data(...
                strcmp(data.Model, model_name) & ...
                data.Gamma_fMSA == 1 & ...
                data.Uptake_Bound == uptake_bound & ...
                data.FVA_Percent == FVA_percentage & ...
                ismember(data.Analysis_Type, analysis_types), :);
        otherwise
            warning('Unknown model: %s', model_name);
            continue;
    end
    
    if ~isempty(valid_data)
        fprintf('  Found %d valid rows\n', height(valid_data));
        all_data = [all_data; valid_data];
    else
        fprintf('  No valid data after filtering\n');
    end
end

% Check if we have data
if isempty(all_data)
    error('No data found after filtering. Check model names and parameters.');
end

fprintf('\n=== DATA SUMMARY ===\n');
fprintf('Total rows: %d\n', height(all_data));
fprintf('Models: %s\n', strjoin(unique(all_data.Model), ', '));
fprintf('Gamma values: %s\n', strjoin(string(unique(all_data.Gamma_fMSA)), ', '));
fprintf('Analysis types: %s\n', strjoin(unique(all_data.Analysis_Type), ', '));

%% GENERATE FIGURES
fprintf('\n=== GENERATING FVA CUMULATIVE DISTRIBUTION PLOTS ===\n');

for m = 1:length(analysis_types)
    analysis_type = analysis_types{m};
    panel_title = panel_titles{m};
    
    fprintf('\nProcessing Panel %d: %s\n', m, analysis_type);
    
    % Filter for this analysis type
    current_data = all_data(...
        strcmp(all_data.Analysis_Type, analysis_type), :);
    
    if isempty(current_data)
        warning('No data for analysis type: %s', analysis_type);
        continue;
    end
    
    fprintf('  Found %d rows\n', height(current_data));
    
    % Create figure
    fig = figure('Units', 'centimeters', 'Position', [1, 1, 20, 13]);
    hold on;
    set(gca, 'FontSize', 14);
    
    color_idx = 1;
    
    % Plot each model and its gamma values
    for i = 1:length(selected_models)
        model_name = selected_models{i};
        model_data = current_data(strcmp(current_data.Model, model_name), :);
        
        if isempty(model_data)
            continue;
        end
        
        model_gamma_values = unique(model_data.Gamma_fMSA);
        fprintf('  %s: gamma values = %s\n', model_name, ...
            strjoin(string(model_gamma_values), ', '));
        
        for j = 1:length(model_gamma_values)
            gamma_val = model_gamma_values(j);
            gamma_data = model_data(model_data.Gamma_fMSA == gamma_val, :);
            
            if isempty(gamma_data)
                continue;
            end
            
            % Extract and validate FVA differences
            fva_diff = gamma_data.FVA_Difference;
            fva_diff = fva_diff(~isnan(fva_diff) & ~isinf(fva_diff) & fva_diff > 0);
            
            if isempty(fva_diff)
                continue;
            end
            
            % Calculate cumulative distribution
            [f, x] = ecdf(fva_diff);
            cumulative_counts = f * length(fva_diff);
            
            % Get color and label
            color = base_colors{color_idx};
            if color_idx <= length(combined_labels)
                label_text = combined_labels{color_idx};
            else
                label_text = sprintf('%s (\\gamma=%.2f)', models_labels{i}, gamma_val);
            end
            
            % Plot
            plot(x, cumulative_counts, '-', 'Color', color, 'LineWidth', 2.5, ...
                'DisplayName', label_text);
            
            fprintf('    γ=%.2f: %d reactions, max=%d\n', gamma_val, ...
                length(fva_diff), max(cumulative_counts));
            
            color_idx = color_idx + 1;
        end
    end
    
    % Customize plot
    title(panel_title, 'FontSize', 18, 'FontWeight', 'bold');
    xlabel('Flux Variability (mmol gDW^{-1} h^{-1})', 'FontSize', 18);
    ylabel('Cumulative Count', 'FontSize', 18);
    legend('Location', 'northwest', 'FontSize', 16);
    
    % Set axes
    grid on;
    box on;
    set(gca, 'XScale', 'log');
    xlim([1e-3, 1e3]);
    ylim([0, inf]);
    
    % Save figure
    panel_tag = sprintf('panel_%c', char(64 + m));  % A, B, C
    filename_base = sprintf('%s_%s', figurePrefix, panel_tag);
    
    saveas(fig, fullfile(output_folder, sprintf('%s.svg', filename_base)), 'svg');
    saveas(fig, fullfile(output_folder, sprintf('%s.png', filename_base)), 'png');
    
    fprintf('  Saved: %s.svg/.png\n', filename_base);
    
    close(fig);
end

%% PANEL D - SPLIT-VIOLIN PLOT BY SUBSYSTEM
fprintf('\n=== GENERATING PANEL D: SUBSYSTEM SPLIT-VIOLINS ===\n');

% Load reaction→subsystem mapping
rxnAttrFile = fullfile(datafolder, 'iML1515_rxns_attrs.csv');
reactionSubsystemMap = containers.Map();

if exist(rxnAttrFile, 'file')
    fprintf('Loading reaction attributes from: %s\n', rxnAttrFile);
    Tmap = readtable(rxnAttrFile);
    if all(ismember({'Rxn_ID','Rxn_Subsystem'}, Tmap.Properties.VariableNames))
        reactionSubsystemMap = containers.Map(cellstr(Tmap.Rxn_ID), cellstr(Tmap.Rxn_Subsystem));
        fprintf('  Loaded %d reaction-subsystem mappings\n', length(reactionSubsystemMap));
    else
        warning('Required columns not found in reaction attributes file');
    end
else
    warning('Reaction attributes file not found: %s', rxnAttrFile);
end

if ~isempty(reactionSubsystemMap)
    % Filter for all_reactions and MAFBA model
    panel_d_data = all_data(...
        strcmp(all_data.Model, 'MAFBA') & ...
        strcmp(all_data.Analysis_Type, 'all_reactions'), :);
    
    if ~isempty(panel_d_data)
        % Compare γ=0.25 vs γ=1.00
        gamma_top = 0.25;
        gamma_bottom = 1.00;
        
        topData = panel_d_data(panel_d_data.Gamma_fMSA == gamma_top, {'Reaction','FVA_Difference'});
        botData = panel_d_data(panel_d_data.Gamma_fMSA == gamma_bottom, {'Reaction','FVA_Difference'});
        
        if ~isempty(topData) && ~isempty(botData)
            fprintf('  Comparing γ=%.2f (%d rxns) vs γ=%.2f (%d rxns)\n', ...
                gamma_top, height(topData), gamma_bottom, height(botData));
            
            % Join data
            J = outerjoin(topData, botData, 'Keys','Reaction', 'MergeKeys',true);
            
            % Handle column names
            if ismember('FVA_Difference_topData', J.Properties.VariableNames)
                FVA_top = J.FVA_Difference_topData;
                FVA_bot = J.FVA_Difference_botData;
            elseif ismember('FVA_Difference_left', J.Properties.VariableNames)
                FVA_top = J.FVA_Difference_left;
                FVA_bot = J.FVA_Difference_right;
            else
                FVA_top = J.FVA_Difference;
                FVA_bot = J.FVA_Difference;
            end
            
            rxNames = string(J.Reaction);
            
            % Compute absolute difference for ranking
            absDiff = abs(FVA_top - FVA_bot);
            
            % Map reactions to subsystems
            rxNorm = regexprep(rxNames, '^rx_', '');
            subs = strings(numel(rxNorm), 1);
            for k = 1:numel(rxNorm)
                rk = char(rxNorm(k));
                if isKey(reactionSubsystemMap, rk)
                    subs(k) = string(reactionSubsystemMap(rk));
                else
                    subs(k) = "Unknown";
                end
            end
            
            % Pick top 5 subsystems by absolute difference
            [~, ord] = sort(absDiff, 'descend', 'MissingPlacement','last');
            pick = strings(0,1);
            for k = 1:numel(ord)
                s = subs(ord(k));
                if ~ismember(s, pick) && s ~= "Unknown"
                    pick(end+1,1) = s; %#ok<AGROW>
                    if numel(pick) >= 5; break; end
                end
            end
            
            fprintf('  Top 5 subsystems: %s\n', strjoin(pick, ', '));
            
            % Create figure
            fig = figure('Units','centimeters','Position',[1,1,26,18]);
            hold on; set(gca,'FontSize',14);
            
            % Colors for split-violin
            col_mango = sscanf('FFA500','%2x%2x%2x',[1 3])/255; % Mango Orange (γ=0.25)
            col_blue  = sscanf('1F56A8','%2x%2x%2x',[1 3])/255; % Blue (γ=1.00)
            
            % Compute global x-limit
            allVar = [FVA_top(:); FVA_bot(:)];
            allVar = allVar(isfinite(allVar));
            xmax_full = max(allVar,[],'omitnan');
            if ~isfinite(xmax_full) || xmax_full<=0, xmax_full = 1; end
            
            % Plot each subsystem
            for si = 1:numel(pick)
                ss = pick(si);
                mask = strcmp(subs, ss);
                vTop_all = FVA_top(mask & isfinite(FVA_top));
                vBot_all = FVA_bot(mask & isfinite(FVA_bot));
                
                % Density for positive values only
                vt = vTop_all(vTop_all > 0);
                vb = vBot_all(vBot_all > 0);
                
                % Upper half (gamma_top = 0.25)
                if ~isempty(vt)
                    [f, xi] = ksdensity(vt, 'Support','positive', 'BoundaryCorrection','reflection');
                    xi = max(xi,0);
                    if max(f) > 0, f = f./max(f) * 0.30; end
                    yc = si*ones(size(xi));
                    X = [xi(:).' fliplr(xi(:).')];
                    Y = [ (yc(:).'+f(:).') fliplr(yc(:).') ];
                    patch(X, Y, col_mango, 'FaceAlpha', 0.5, 'EdgeColor', col_mango*0.8);
                    medv = median(vt,'omitnan');
                    plot([medv medv], [si si+0.18], '-', 'Color',[0.2 0.2 0.2], 'LineWidth',1.2);
                end
                jit = (rand(size(vTop_all))-0.5)*0.12;
                scatter(vTop_all, si + 0.18 + jit, 10, 'filled', ...
                    'MarkerFaceAlpha',0.25, 'MarkerFaceColor', col_mango, 'MarkerEdgeColor','none');
                
                % Lower half (gamma_bottom = 1.00)
                if ~isempty(vb)
                    [f2, xi2] = ksdensity(vb, 'Support','positive', 'BoundaryCorrection','reflection');
                    xi2 = max(xi2,0);
                    if max(f2) > 0, f2 = f2./max(f2) * 0.30; end
                    yc = si*ones(size(xi2));
                    X2 = [xi2(:).' fliplr(xi2(:).')];
                    Y2 = [ (yc(:).'-f2(:).') fliplr(yc(:).') ];
                    patch(X2, Y2, col_blue, 'FaceAlpha', 0.5, 'EdgeColor', col_blue*0.8);
                    medv2 = median(vb,'omitnan');
                    plot([medv2 medv2], [si-0.18 si], '-', 'Color',[0.2 0.2 0.2], 'LineWidth',1.2);
                end
                jit2 = (rand(size(vBot_all))-0.5)*0.12;
                scatter(vBot_all, si - 0.18 + jit2, 10, 'filled', ...
                    'MarkerFaceAlpha',0.25, 'MarkerFaceColor', col_blue, 'MarkerEdgeColor','none');
            end
            
            % Customize plot
            xlim([0, xmax_full*1.05]);
            xlabel('Flux Variability (mmol gDW^{-1} h^{-1})', 'FontSize', 18);
            ylabel('Subsystem', 'FontSize', 18);
            yticks(1:numel(pick));
            yticklabels(pick);
            ylim([0.5, numel(pick)+0.5]);
            grid on; box on;
            title('Variability Distributions by Subsystem', 'FontSize', 18, 'FontWeight', 'bold');
            ax = gca; ax.XAxis.FontSize = 16;
            
            % Legend
            h1 = plot(nan,nan,'s','MarkerFaceColor',col_mango,'MarkerEdgeColor','none');
            h2 = plot(nan,nan,'s','MarkerFaceColor',col_blue,'MarkerEdgeColor','none');
            lg = legend([h1 h2], ...
                {'Both constraints, \gamma = 0.25','No Membrane constraint, \gamma = 1.00'}, ...
                'Location','southoutside','Orientation','horizontal','Box','off');
            lg.FontSize = 16;
            
            % Save Panel D
            filename_base = sprintf('%s_panel_D', figurePrefix);
            saveas(fig, fullfile(output_folder, sprintf('%s.svg', filename_base)), 'svg');
            saveas(fig, fullfile(output_folder, sprintf('%s.png', filename_base)), 'png');
            
            fprintf('  Saved: %s.svg/.png\n', filename_base);
            close(fig);
        else
            warning('Missing data for Panel D gamma comparison');
        end
    else
        warning('No data found for Panel D');
    end
else
    warning('Skipping Panel D: reaction-subsystem mapping not available');
end

%% COMPLETION
fprintf('\n=== FIGURE 6 COMPLETE ===\n');
fprintf('Generated figures:\n');
fprintf('  Models compared:\n');
fprintf('    - %s (multiple γ values)\n', selected_models{1});
fprintf('    - %s (γ=1.0 only)\n', selected_models{2});
fprintf('  Uptake bound: %d mmol gDW^{-1} h^{-1}\n', uptake_bound);
fprintf('  Total panels: 4 (A, B, C, D)\n');
fprintf('  Total files: 4 panels × 2 formats = 8 files\n');
fprintf('  Output location: %s\n', output_folder);
