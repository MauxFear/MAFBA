function [solutions, ec_models] = loadEcModel(gammas, varargin)
    % Set up the MAFBA model solver (important to use CPLEX)
    changeCobraSolver('ibm_cplex');
    if ~exist('gammas', 'var')
       gammas = [1,0.25]; 
    end
    
    % Default values
    defaults = struct(...
        'modelV', 0, ...
        'dateVersion', 'AlterOriginal', ...
        'saveMATLAB', 0, ...
        'saveSBML', 0, ...
        'opt', 1, ...
        'alearn', 0.1, ...
        'sigma_mode', 'parametrized_totalweight', ...
        'save_tag', '', ...
        'modelCode', '' ...
    );
    
    % Parse inputs
    p = inputParser;
    fields = fieldnames(defaults);
    for i = 1:length(fields)
        addOptional(p, fields{i}, defaults.(fields{i}));
    end
    parse(p, varargin{:});
    
    % Load Prepared Model
    modelID = 'iML1515';


    % Retrieve values
    params = p.Results;
    if isempty(params.modelCode)
        params.modelCode = deriveModelCode(modelID, params.dateVersion, params.sigma_mode);
    end
    

    
    % Construct paths
    baseDir = fileparts(mfilename('fullpath'));
    repoRoot = fullfile(baseDir, '..', '..');
    outputRoot = fullfile(repoRoot, '1_Model_Construction', '03_MAFBA_Models_Outputs', params.modelCode);
    modelsRoot = fullfile(outputRoot, 'models');
    attributesRoot = fullfile(outputRoot, 'attributes_spreadsheets');

    ensureDir(modelsRoot);
    ensureDir(attributesRoot);
    input_file = sprintf('prepared_%s_modelV%d_%s.mat', modelID, params.modelV, params.dateVersion);
    input_path = fullfile(repoRoot, '1_Model_Construction', '02_MATLAB_Build', 'output', input_file);
    
    % Load model with error handling
    if ~isfile(input_path)
        warning('Model file not found at: %s', input_path);
        input_path = fullfile(outputRoot, 'model_OutputData', ...
            sprintf('annotated_%s_irreversible_V%d_%s.xml', modelID, params.modelV, params.dateVersion));
    end
    
    % Load model based on file type
    if isMatlabDataFile(input_path)
        loadedData = load(input_path);
        modelFields = fieldnames(loadedData);
        modelField = modelFields{contains(modelFields, 'model')};
        if isempty(modelField)
            error('No model field found in loaded data');
        end
        model = loadedData.(modelField);
    else
        model = readCbModel(input_path);
    end

    % Check Model Annotations
    requiredFields = {'location', 'mw_c', 'mw_m', 'kcat', 'number_tMb'};
    if ~all(isfield(model, requiredFields))
        warning('Model missing required annotations. Preparing model...');
        attributesFileName = sprintf('%s_irreversible_reactionDataV%d_%s.csv', modelID, params.modelV, params.dateVersion);
        attributesPath = fullfile(outputRoot, 'models_InputData', attributesFileName);
        savePath = fullfile(repoRoot, '1_Model_Construction', '02_MATLAB_Build', 'output', ...
            sprintf('prepared_%s_modelV%d_%s.mat', modelID, params.modelV, params.dateVersion));
        
        model = prepareModelfromAnnotations(model, attributesPath, true, savePath);
        generateModelExcel(input_path, attributesRoot, params.dateVersion);
    end

    % MAFBA Model Configuration
    addZYprot = 1;
    sigma_param = 1;    
    gamma = gammas(1); 
    [solution, ec_model] = setUpModel_v1(model, params.sigma_mode, addZYprot, sigma_param, gamma);

    % Optimization of sigma
    if params.opt == 1
        disp('Optimizing w parameter');
        initial_sigma = getInitialSigma(params.modelV);
        opt_sigma_param = getOpt_sigma_param(initial_sigma, model, params.alearn, params.sigma_mode, gamma);
        fprintf('Optimal sigma parameter: %1.4f\n', opt_sigma_param);
    else
        opt_sigma_param = getDefaultSigma(params.modelV);
        fprintf('Using default w parameter: %1.4f\n', opt_sigma_param);
        [solution, ec_model] = setUpModel_v1(model, params.sigma_mode, addZYprot, opt_sigma_param, gamma);
    end

    % Process multiple gammas
    solutions = cell(1, length(gammas));
    ec_models = cell(1, length(gammas));
    gammaStrings = arrayfun(@formatGamma, gammas, 'UniformOutput', false);
    
    for i = 1:length(gammas)
        gamma = gammas(i);
        gammaString = gammaStrings{i};
        fprintf('Loading ecModel with gamma: %1.4f\n', gamma);
        
        [solutions{i}, ec_models{i}] = setUpModel_v1(model, params.sigma_mode, addZYprot, opt_sigma_param, gamma);
        
        if ~isempty(params.save_tag)
            model_id_suffix = ['_' params.save_tag];
        else
            model_id_suffix = '';
        end
        ec_models{i}.modelID = sprintf('ecModel_%s_%sV%d_g_%s%s', modelID, params.dateVersion, params.modelV, gammaString, model_id_suffix);

        gammaFolder = fullfile(modelsRoot, ['g' gammaString]);
        ensureDir(gammaFolder);
        
        if params.saveSBML
            fileName = fullfile(gammaFolder, sprintf('mafba_%s_ecV%d_g%s_%s%s.xml', modelID, params.modelV, gammaString, params.dateVersion, model_id_suffix));
            tagExistingFile(fileName);
            writeCbModel(ec_models{i}, 'format', 'sbml', 'fileName', fileName);
        end
        
        if params.saveMATLAB
            fileName = fullfile(gammaFolder, sprintf('mafba_%s_ecV%d_g%s_%s%s.mat', modelID, params.modelV, gammaString, params.dateVersion, model_id_suffix));
            tagExistingFile(fileName);
            ec_model = ec_models{i};
            save(fileName, 'ec_model');
        end
    end
end

function initial_sigma = getInitialSigma(modelV)
    switch modelV
        case 3
            initial_sigma = 8.4510;
        case 2
            initial_sigma = 3.23;
        case 1
            initial_sigma = 14.018;
        case 0
            initial_sigma = 0.000354;
        otherwise
            error('Invalid model version');
    end
end

function default_sigma = getDefaultSigma(modelV)
    switch modelV
        case 3
            default_sigma = 8.4701;
        case 2
            default_sigma = 3.2320;
        case 1
            default_sigma = 14.018;
        case 0
            default_sigma = 0.4048;
        otherwise
            error('Invalid model version');
    end
end

function gammaString = formatGamma(gamma)
    if abs(gamma - round(gamma)) < 1e-9
        gammaString = sprintf('%d', round(gamma));
    else
        gammaString = sprintf('%.12g', gamma);
        if startsWith(gammaString, '.')
            gammaString = ['0' gammaString];
        end
    end
end

function modelCode = deriveModelCode(modelID, dateVersion, sigmaMode)
    if strcmpi(sigmaMode, 'homogeneous')
        suffix = 'Homogeneous';
    elseif strcmpi(sigmaMode, 'parametrized_TMbw')
        suffix = 'AlterTMbw';
    else
        suffix = dateVersion;
    end
    modelCode = sprintf('%s_%s', modelID, suffix);
end

function isMatlabDataFile = isMatlabDataFile(input_path)
    [~, ~, extension] = fileparts(input_path);
    isMatlabDataFile = strcmpi(extension, '.mat');
end

function opt_sigma_param = getOpt_sigma_param(sigma_param, model,alearn,sigma_mode, gamma)
    %opt_objective_min minimize the difference between a target and the growth
    %rate obtained with a certain weight parameter
    %   Detailed explanation goes here

    target = 1; % max growth rate allowed
    
    
    % model_iMLirr = readCbModel('D:\Code\MATLAB\Projects\MAFBA\Utility\Models\modified_model_iML1515_may0523.mat');
    addZYprot = 1;
    [solution,~] = setUpModel_v1(model,sigma_mode,addZYprot,sigma_param, gamma);
%     % Printing the upper boundary of the glucose exchange
%     ex_glc_r_b = find(strcmp(model.rxns, 'EX_glc__D_e_b'));
%     fprintf('4 The upper boundary of the glucose exchange is: %f\n', model.ub(ex_glc_r_b));
    
    % [growthRate, ~, ~] = simpleOptECmodel(pre_ec_model,w_param, w_mode);
    growthRate = solution.f;
    objective = abs(target-growthRate);
    tolerance = 0.00005;
    i = 0;
    while objective > tolerance && i < 100
        i=i+1;
%         disp(i)
        dif_off = target-growthRate;
        % Adjust the weight parameter
        sigma_param = sigma_param + (-alearn*dif_off);

        fprintf('The w_param value is: %f  \n', sigma_param);
        % Re-optimize the model
        [solution,~] = setUpModel_v1(model,sigma_mode,addZYprot,sigma_param, gamma);
        growthRate = solution.f;

        objective = abs(target-growthRate);
    %     get_phimax_values(ec_model);
    %     disp(growthRate)

    end

    opt_sigma_param = sigma_param;
    disp(growthRate)
end

function ensureDir(pathToCreate)
    if ~exist(pathToCreate, 'dir')
        mkdir(pathToCreate);
    end
end

function tagExistingFile(filePath)
    if ~isfile(filePath)
        return;
    end

    [folder, baseName, ext] = fileparts(filePath);
    backupPath = fullfile(folder, sprintf('%s_backup%s', baseName, ext));
    if ~isfile(backupPath)
        movefile(filePath, backupPath);
        return;
    end

    idx = 1;
    while true
        candidate = fullfile(folder, sprintf('%s_backup_%02d%s', baseName, idx, ext));
        if ~isfile(candidate)
            movefile(filePath, candidate);
            return;
        end
        idx = idx + 1;
    end
end

function phimax_values = get_phimax_values(ec_model)
    %get_phimax_values get to 
    %   Detailed explanation goes here
    

    phi_max=find(strcmp(ec_model.rxns,'Prot_Const'));
    phi_Mmax=find(strcmp(ec_model.rxns,'Memb_Const'));

    phimax_values = [ec_model.ub(phi_max),ec_model.ub(phi_Mmax)];

    fprintf('Upper boundary for phimax:  %1.4f \n', ec_model.ub(phi_max));
    fprintf('Lower boundary for phimax:  %1.4f \n', ec_model.lb(phi_max));
    fprintf('Upper boundary for phiMmax:  %1.4f \n', ec_model.ub(phi_Mmax));
    fprintf('Lower boundary for phiMmax:  %1.4f \n', ec_model.lb(phi_Mmax));
    fprintf('Gamma ratio:  %1.4f \n', ec_model.ub(phi_Mmax)/ec_model.ub(phi_max));
end
