% In this script we increase the upper boundary of the reverse exchange 
% reaction of glucose to control the growth rate as the conventional
% method in Genome scale modeling. 
% This version focuses on sensitivity analysis of ALL reactions.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Fine grain approach - All Reactions                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% LOAD MODEL

% Set up of the MAFBA model
changeCobraSolver('gurobi'); %Using cplex !!Important to do this

%%
modelV = 0;
dateVersion = 'MAFBA_homogeneous';
if exist('ec_model', 'var') == 1
    ec_model = setUp_GlucoseGrowth(ec_model,dateVersion,modelV);
else
    ec_model = setUp_GlucoseGrowth(0,dateVersion,modelV);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   MODEL SETTING                                           %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Remove default bounds for glucose, oxygen and CO2
% Get the indexes for each reaction forward or reverse
ex_glc_r_f=find(strcmp(ec_model.rxns,'EX_glc__D_e_f'));
ex_glc_r_b=find(strcmp(ec_model.rxns,'EX_glc__D_e_b'));
ex_o2_r_f=find(strcmp(ec_model.rxns,'EX_o2_e_f'));
ex_o2_r_b=find(strcmp(ec_model.rxns,'EX_o2_e_b'));
ex_co2_r_f=find(strcmp(ec_model.rxns,'EX_co2_e_f'));
ex_co2_r_b=find(strcmp(ec_model.rxns,'EX_co2_e_b'));
% oxygen
ec_model.lb(ex_o2_r_b)=0; ec_model.ub(ex_o2_r_b)=20; %reverse
ec_model.ub(ex_o2_r_f)=0; ec_model.ub(ex_o2_r_f)=1000; %forward
% CO2
ec_model.lb(ex_co2_r_b)=0; ec_model.ub(ex_co2_r_b)=0; %reverse
ec_model.ub(ex_co2_r_f)=0; ec_model.ub(ex_co2_r_f)=1000; %forward
% Glucose
ec_model.lb(ex_glc_r_b)=0; ec_model.ub(ex_glc_r_b)=30; %reverse
ec_model.lb(ex_glc_r_f)=0; ec_model.ub(ex_glc_r_f)=1000; %forward

% Test the ecModel
solution = get_solution_ec_model(ec_model);
[phimax_values, phi_indices ]= get_phimax_values(ec_model);

%% Simulation Parameters
% Define analysis description and parameters
analysis_description = 'Sensitivity analysis of kcat values for all reactions through scaling of proteome constraint coefficients';
analysis_type = 'Coefficient scaling - All reactions';

% Add parameter for starting batch (default is 1 for new analysis)
start_batch = 1; % Change this value to continue from a specific batch

Npoints = 100;
% Allocate 30% of points below 10 and 70% above 10 (up to 35)
num_points_below = round(0.3 * Npoints);
num_points_above = Npoints - num_points_below;

% Generate points below 10 (from 0.2 to 9.8)
points_below_10 = linspace(0.1, 9.8, num_points_below);

% Generate points above 10 (from 10 to 35)
% Using square function to condense higher values
points_above_10 = 10 + 25 * (linspace(0, 1, num_points_above).^2);

% Combine the points
exbValues = [points_below_10, points_above_10];

sub_r=find(strcmp(ec_model.rxns,'EX_glc__D_e_b'));
bm_r = ec_model.protGroup(3).rxns;
ace_r=find(strcmp(ec_model.rxns,'EX_ac_e_f'));
resp_r=find(strcmp(ec_model.rxns,'ATPS4rpp_f'));

% indices reaction phi
phi_max_r = phi_indices(1);  % phi_max -- Overall Proteome Allocation
phi_Mmax_r = phi_indices(2); % phi_Mmax -- Membrane Allocation
phi_max_value = ec_model.ub(phi_max_r); % phi_max fixed value

% Define gamma values for sensitivity analysis (0.1 to 1)
% Total of 50 values: 35 below 0.3 and 15 above 0.3
total_gamma_values = 50;
num_values_below = 35;
num_values_above = 15;

% Generate values below 0.3 (from 0.1 to 0.3)
% Using more points near 0.1 with square root distribution
gamma_below = 0.1 + 0.2 * sqrt(linspace(0, 1, num_values_below));

% Generate values above 0.3 (from 0.3 to 1.0)
% Using more points near 0.3 with square distribution
gamma_above = 0.3 + 0.7 * (linspace(0, 1, num_values_above).^2);

% Combine the values
gammaValues = [gamma_below, gamma_above];
gammaValues = [1, 0.3, 0.25, 0.2, 0.15];
atpm_values = [6.89];

% Define coefficient change factors
change_factors = [1];
% change_factors = [1, 2, 5, 10, 20, 50, 1/2, 1/5, 1/10, 1/20, 1/50];

num_iterations = 1; % One iteration per change factor
batch_size = 5; % Process one change factor at a time

% Create timestamp-based folder structure
timestamp = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
mainFolder = 'simResults';
analysisFolder = sprintf('%s_all_reactions', dateVersion);
runFolder = fullfile(mainFolder, analysisFolder, timestamp);
datafolder = fullfile(runFolder, 'simData');

% Create folder structure
if ~exist(mainFolder, 'dir'), mkdir(mainFolder); end
if ~exist(fullfile(mainFolder, analysisFolder), 'dir'), mkdir(fullfile(mainFolder, analysisFolder)); end
if ~exist(runFolder, 'dir'), mkdir(runFolder); end
if ~exist(datafolder, 'dir'), mkdir(datafolder); end

% Save detailed analysis parameters and description
fid = fopen(fullfile(runFolder, 'analysis_metadata.txt'), 'w');
fprintf(fid, 'Analysis Description:\n%s\n\n', analysis_description);
fprintf(fid, 'Analysis Type: %s\n', analysis_type);
fprintf(fid, 'Timestamp: %s\n\n', timestamp);
fprintf(fid, 'Model Parameters:\n');
fprintf(fid, '- Model Version: %d\n', modelV);
fprintf(fid, '- Date Version: %s\n', dateVersion);
fprintf(fid, '- Number of Points: %d\n', Npoints);
fprintf(fid, '- Gamma Values: %s\n', mat2str(gammaValues));
fprintf(fid, '- ATPM Values: %s\n', mat2str(atpm_values));
fprintf(fid, '\nScaling Parameters:\n');
fprintf(fid, '- Change Factors: %s\n', mat2str(change_factors));
fprintf(fid, '- Total Iterations: %d\n', num_iterations);
fprintf(fid, '- Batch Size: %d\n', batch_size);

% Get original costs
original_costC = ec_model.S(end-1,:);
original_costM = ec_model.S(end,:);

atpm_index = find(strcmp(ec_model.rxns, 'ATPM'));

%% IDENTIFY ALL REACTIONS TO SCALE
% Exclude specific reactions from scaling
excluded_rxns = [ec_model.rxns(end-7:end); ec_model.rxns(ec_model.protGroup(3).rxns)];
all_rxns_idx = 1:length(ec_model.rxns);
rxnsToScale = all_rxns_idx(~ismember(all_rxns_idx, find(ismember(ec_model.rxns, excluded_rxns))));

% Log the number of reactions to scale
fprintf(fid, '\nReactions to Scale:\n');
fprintf(fid, '- Total model reactions: %d\n', length(ec_model.rxns));
fprintf(fid, '- After exclusions: %d\n', length(rxnsToScale));
fprintf(fid, '\nList of excluded reactions:\n');
for i = 1:length(excluded_rxns)
    fprintf(fid, '%s\n', excluded_rxns{i});
end
fclose(fid);

%% STARTING SIMULATIONS
% Initialize log file
logFile = fullfile(runFolder, 'simulation_log.txt');
logging = @(msg) write_to_log(logFile, msg);
logging('Starting simulation...');
logging(sprintf('Identified %d reactions for sensitivity analysis', length(rxnsToScale)));

% Initialize tables for storing results
% Add 'rx_' prefix to reaction IDs to avoid variable names starting with numbers
rxn_var_names = strcat('rx_', ec_model.rxns);
varTypes = repmat({'double'}, 1, numel(ec_model.rxns) + 4); % +1 for Change_Factor
dataTable = table('Size', [0, numel(ec_model.rxns)+5], ...
    'VariableTypes', [varTypes, {'double'}], ...
    'VariableNames', ["Gamma", "Iteration", "Change_Factor", "Uptake_Bound", rxn_var_names', "Random_Seed"]);

aceLTable = table('Size', [0, 6], ...
    'VariableTypes', {'double', 'double', 'double', 'double', 'double', 'double'}, ...
    'VariableNames', ["Gamma", "Iteration", "Change_Factor", "Acetate_threshold", "ATPM", "Random_Seed"]);

sectorsTable = table('Size', [0, 21], ...
    'VariableTypes', repmat({'double'}, 1, 21), ...
    'VariableNames', ["Gamma", "Iteration", "Change_Factor", "Uptake_Bound", "Growth", "ATPM", "Random_Seed", "phiCm", "phiCc", "phiR", "phiEc", "phiEr",...
        "phiEm", "phiP_m", "phiP", "phiM", "phiCm_m", "phiEr_m", "phiEm_m", "phiPmax", "phiMmax"]);

% Modify modelDataTable to include both scaling factors and actual coefficients
modelDataTable = table('Size', [0, 4 + 3*numel(ec_model.rxns)], ...
    'VariableTypes', [{'double', 'double', 'double', 'double'}, repmat({'double'}, 1, 3*numel(ec_model.rxns))], ...
    'VariableNames', ["Gamma", "Iteration", "Change_Factor", "Random_Seed", ...
                     strcat('rx_', ec_model.rxns, '_scale')', ...  % Scaling factors with rx_ prefix
                     strcat('rx_', ec_model.rxns, '_costC')', ...  % Cytosolic proteome coefficients with rx_ prefix
                     strcat('rx_', ec_model.rxns, '_costM')']);    % Membrane proteome coefficients with rx_ prefix

% Calculate total steps for progress tracking
total_steps = numel(atpm_values) * numel(gammaValues) * numel(change_factors) * Npoints;
step_count = 0;
batch_count = 1;

% Calculate which iteration to start from based on start_batch
start_iteration = ((start_batch - 1) * batch_size) + 1;

% Update the log file to indicate continuation if starting from later batch
if start_batch > 1
    logging(sprintf('Continuing analysis from batch %d (iteration %d)', start_batch, start_iteration));
end

% Initialize batch counter to start_batch
batch_count = start_batch;

try
    for atpm_idx = 1:numel(atpm_values)
        atpm_val = atpm_values(atpm_idx);
        logging(sprintf('Processing ATPM Value: %.2f (%d/%d)', atpm_val, atpm_idx, numel(atpm_values)));
        ec_model.ub(atpm_index) = atpm_val;
        
        for f = 1:numel(gammaValues)
            gamma = gammaValues(f);
            logging(sprintf('Processing Gamma: %.2f (%d/%d)', gamma, f, numel(gammaValues)));
            
            % Set phi_Mmax value
            phi_Mmax_value = phi_max_value*gamma;
            ec_model.ub(phi_Mmax_r) = phi_Mmax_value;

            % Process in batches, starting from the specified batch
            for batch_start = start_iteration:batch_size:numel(change_factors)
                batch_end = min(batch_start + batch_size - 1, numel(change_factors));
                logging(sprintf('Starting batch %d (change factors %d-%d)', batch_count, batch_start, batch_end));
                
                % Reset start_iteration after first gamma value to ensure we process all iterations
                % for subsequent gamma values
                if batch_start == start_iteration && (f > 1 || atpm_idx > 1)
                    start_iteration = 1;
                end

                for cf_idx = batch_start:batch_end
                    change_factor = change_factors(cf_idx);
                    logging(sprintf('Processing change factor %.4f (%d/%d)', change_factor, cf_idx, numel(change_factors)));
                    
                    % Store model coefficients for this iteration (initialize all to 1, then set reactions to scale)
                    modelCoeffs = ones(1, numel(ec_model.rxns));
                    modelCoeffs(rxnsToScale) = change_factor;
                    
                    % Get current proteome constraint coefficients
                    costC = ec_model.S(end-1,:);  % Cytosolic proteome coefficients
                    costM = ec_model.S(end,:);    % Membrane proteome coefficients
                    
                    % Combine all data for modelDataTable
                    modelRow = [gamma, 1, change_factor, 1, ...
                              modelCoeffs, ...    % Scaling factors
                              costC, ...          % Current cytosolic proteome coefficients
                              costM];            % Current membrane proteome coefficients
                    modelDataTable = vertcat(modelDataTable, array2table(modelRow, 'VariableNames', modelDataTable.Properties.VariableNames));
                    
                    % Reset to original costs first
                    ec_model.S(end-1,:) = original_costC;
                    ec_model.S(end,:) = original_costM;
                    
                    % Apply scaling to all reactions to scale
                    for rx_idx = 1:length(rxnsToScale)
                        rx = rxnsToScale(rx_idx);
                        ec_model.S(end-1,rx) = original_costC(rx) * change_factor;
                        ec_model.S(end,rx) = original_costM(rx) * change_factor;
                    end

                    % Simulate for each uptake bound
                    for i = 1:Npoints
                        step_count = step_count + 1;
                        if mod(step_count, 100) == 0
                            logging(sprintf('Progress: %.1f%% complete', (step_count/total_steps)*100));
                        end
                        
                        exlb = exbValues(i);
                        ec_model.ub(sub_r) = exlb;
                        solution = get_pFBAsolution_ec_model(ec_model);

                        % Store flux values
                        fluxValues = solution.x';
                        biomassValues(i) = solution.x(bm_r);
                        acetateValues(i) = solution.x(ace_r);
                        respirationValues(i) = solution.x(resp_r);
                        atpM_sol = solution.x(atpm_index);

                        % Store results in tables
                        rowValues = [gamma, 1, change_factor, exlb, fluxValues(:)', 1];
                        dataTable = vertcat(dataTable, array2table(rowValues, 'VariableNames', dataTable.Properties.VariableNames));

                        % Store sectors data
                        sectorValues = [gamma, 1, change_factor, exlb, solution.x(bm_r), atpM_sol, 1, ...
                            solution.phiCm, solution.phiCc, solution.phiR, solution.phiEc, solution.phiEr,...
                            solution.phiEm, solution.phiP_m, solution.phiP, solution.phiM, solution.phiCm_m,...
                            solution.phiEr_m, solution.phiEm_m, solution.phiPmax, solution.phiMmax];
                        sectorsTable = vertcat(sectorsTable, array2table(sectorValues, 'VariableNames', sectorsTable.Properties.VariableNames));
                    end

                    % Calculate acetate threshold for this iteration
                    T = table(biomassValues, acetateValues, respirationValues);
                    T = sortrows(T, 'biomassValues');

                    % Detect acetate threshold
                    dif_v = diff(T.respirationValues);
                    dif_g = diff(T.biomassValues);
                    slope_v = dif_v./dif_g;
                    m_s = max(slope_v);
                    memb_s = find(slope_v <= m_s*0.01);
                    t_i = find(T.acetateValues(memb_s) >= 0.01, 1, 'first');
                    ac_L = memb_s(t_i);

                    if isempty(ac_L)
                        memb_s = find(slope_v <= m_s*0.05);
                        t_i = find(T.acetateValues(memb_s) >= 0, 1, 'first');
                        ac_L = memb_s(t_i);
                    end
                    if isempty(ac_L)
                        memb_s = find(slope_v <= m_s*0.2);
                        t_i = find(T.acetateValues(memb_s) >= 0, 1, 'first');
                        ac_L = memb_s(t_i);
                    end
                    if isempty(ac_L)
                        ac_L = 1;
                    end

                    growth_ac = T.biomassValues(ac_L);

                    % Store acetate threshold
                    newRow = table(gamma, 1, change_factor, growth_ac, atpm_val, 1,...
                        'VariableNames', ["Gamma", "Iteration", "Change_Factor", "Acetate_threshold", "ATPM", "Random_Seed"]);
                    aceLTable = vertcat(aceLTable, newRow);
                end

                % Save batch results
                batch_suffix = sprintf('_batch%d', batch_count);
                save_batch_results(datafolder, modelV, dateVersion, batch_suffix, ...
                    dataTable, aceLTable, sectorsTable, modelDataTable);
                
                logging(sprintf('Batch %d completed and saved', batch_count));
                batch_count = batch_count + 1;
            end
        end
    end

    % Calculate and save final statistics
    logging('Calculating final statistics...');
    
    % Calculate averages and standard deviations
    avg_dataTable = grpstats(dataTable, {'Gamma', 'Change_Factor', 'Uptake_Bound'}, 'mean');
    avg_aceLTable = grpstats(aceLTable, {'Gamma', 'Change_Factor', 'ATPM'}, 'mean');
    avg_sectorsTable = grpstats(sectorsTable, {'Gamma', 'Change_Factor', 'Uptake_Bound'}, 'mean');
    avg_modelDataTable = grpstats(modelDataTable, {'Gamma', 'Change_Factor'}, 'mean');
    
    stdv_dataTable = grpstats(dataTable, {'Gamma', 'Change_Factor', 'Uptake_Bound'}, 'std');
    stdv_aceLTable = grpstats(aceLTable, {'Gamma', 'Change_Factor', 'ATPM'}, 'std');
    stdv_sectorsTable = grpstats(sectorsTable, {'Gamma', 'Change_Factor', 'Uptake_Bound'}, 'std');
    stdv_modelDataTable = grpstats(modelDataTable, {'Gamma', 'Change_Factor'}, 'std');

    % Save final results
    save_final_results(datafolder, modelV, dateVersion, ...
        dataTable, aceLTable, sectorsTable, modelDataTable, ...
        avg_dataTable, avg_aceLTable, avg_sectorsTable, avg_modelDataTable, ...
        stdv_dataTable, stdv_aceLTable, stdv_sectorsTable, stdv_modelDataTable);

    % Save workspace
    save(fullfile(runFolder, 'workspace.mat'));
    
    logging('Simulation completed successfully!');
    fprintf('Simulation completed successfully!\n');
    
catch ME
    logging(sprintf('Error occurred: %s', ME.message));
    rethrow(ME);
end

% Helper function to write to log file
function write_to_log(logFile, message)
    timestamp = datestr(now, 'yyyy-mm-dd HH:MM:SS');
    fprintf('%s: %s\n', timestamp, message);
    fid = fopen(logFile, 'a');
    fprintf(fid, '%s: %s\n', timestamp, message);
    fclose(fid);
end

% Helper function to save batch results
function save_batch_results(datafolder, modelV, dateVersion, batch_suffix, ...
    dataTable, aceLTable, sectorsTable, modelDataTable)
    
    baseFileName = sprintf('modelV%i_all_reactions_%s%s', modelV, dateVersion, batch_suffix);
    
    writetable(dataTable, fullfile(datafolder, ['data_' baseFileName '.csv']));
    writetable(aceLTable, fullfile(datafolder, ['aceL_' baseFileName '.csv']));
    writetable(sectorsTable, fullfile(datafolder, ['sectors_' baseFileName '.csv']));
    writetable(modelDataTable, fullfile(datafolder, ['model_' baseFileName '.csv']));
end

% Helper function to save final results
function save_final_results(datafolder, modelV, dateVersion, ...
    dataTable, aceLTable, sectorsTable, modelDataTable, ...
    avg_dataTable, avg_aceLTable, avg_sectorsTable, avg_modelDataTable, ...
    stdv_dataTable, stdv_aceLTable, stdv_sectorsTable, stdv_modelDataTable)
    
    baseFileName = sprintf('modelV%i_all_reactions_%s_final', modelV, dateVersion);
    
    % Save raw data
    writetable(dataTable, fullfile(datafolder, ['data_' baseFileName '.csv']));
    writetable(aceLTable, fullfile(datafolder, ['aceL_' baseFileName '.csv']));
    writetable(sectorsTable, fullfile(datafolder, ['sectors_' baseFileName '.csv']));
    writetable(modelDataTable, fullfile(datafolder, ['model_' baseFileName '.csv']));
    
    % Save averages
    writetable(avg_dataTable, fullfile(datafolder, ['avg_data_' baseFileName '.csv']));
    writetable(avg_aceLTable, fullfile(datafolder, ['avg_aceL_' baseFileName '.csv']));
    writetable(avg_sectorsTable, fullfile(datafolder, ['avg_sectors_' baseFileName '.csv']));
    writetable(avg_modelDataTable, fullfile(datafolder, ['avg_model_' baseFileName '.csv']));
    
    % Save standard deviations
    writetable(stdv_dataTable, fullfile(datafolder, ['stdv_data_' baseFileName '.csv']));
    writetable(stdv_aceLTable, fullfile(datafolder, ['stdv_aceL_' baseFileName '.csv']));
    writetable(stdv_sectorsTable, fullfile(datafolder, ['stdv_sectors_' baseFileName '.csv']));
    writetable(stdv_modelDataTable, fullfile(datafolder, ['stdv_model_' baseFileName '.csv']));
end

% Keep the original get_phimax_values function  
function [phimax_values, phi_indices] = get_phimax_values(ec_model)  
    phi_max_index=find(strcmp(ec_model.rxns,'Prot_Const'));  
    phi_Mmax_index=find(strcmp(ec_model.rxns,'Memb_Const'));  

    phimax_values = [ec_model.ub(phi_max_index),ec_model.ub(phi_Mmax_index)];  
    phi_indices = [phi_max_index,phi_Mmax_index];  

    fprintf('Upper boundary for phimax:  %1.4f \n', ec_model.ub(phi_max_index));  
    fprintf('Lower boundary for phimax:  %1.4f \n', ec_model.lb(phi_max_index));  
    fprintf('Upper boundary for phiMmax:  %1.4f \n', ec_model.ub(phi_Mmax_index));  
    fprintf('Lower boundary for phiMmax:  %1.4f \n', ec_model.lb(phi_Mmax_index));  
    fprintf('Gamma ratio:  %1.4f \n', ec_model.ub(phi_Mmax_index)/ec_model.ub(phi_max_index));  
end