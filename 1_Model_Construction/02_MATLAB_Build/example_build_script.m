%% 
% Initialize variables
saveBool = true;
modelV = 1;
dateVersion = 'MAFBA';

% Define file names
modelFileName = sprintf('annotated_iML1515_irreversible_V%d_%s.xml', modelV, dateVersion);
attributesFileName = sprintf('iML1515_irreversible_reactionDataV%d_%s.csv', modelV, dateVersion);
outputModelName = sprintf('prepared_iML1515_modelV%d_%s.mat', modelV, dateVersion);

% Define relative paths
% Assuming script is run from MAFBA_Paper_Public_Materials/1_Model_Construction/02_MATLAB_Build/
baseDir = fileparts(mfilename('fullpath'));
modelCode = sprintf('iML1515_%s', dateVersion);
pythonOutputDir = fullfile(baseDir, '..', '03_MAFBA_Models_Outputs', modelCode);
utilityDir = fullfile(baseDir, '..', '..', 'utils');
inputDataDir = fullfile(baseDir, '..', '..', 'data', 'input');
outputDir = fullfile(baseDir, 'output');

% Add utility path
addpath(utilityDir);

% Construct full file paths
modelPath = fullfile(pythonOutputDir, 'model_OutputData', modelFileName);
attributesPath = fullfile(pythonOutputDir, 'models_InputData', attributesFileName);
savePath = fullfile(outputDir, outputModelName);
excelUpdatePath = '';

% Create output directory if it doesn't exist
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

%% Run Model Preparation
fprintf('Loading and preparing model from:\n %s\n', modelPath);
model_iMLirr = prepareModel(modelPath, attributesPath, saveBool, savePath, excelUpdatePath);

%% Post-processing
% Convert rules to human readable format
grRules = convertRulesWithGenes(model_iMLirr.rules, model_iMLirr.genes);
model_iMLirr.grRules = grRules;

% Generate Excel with Model's attributes
fprintf('Generating Excel export...\n');
generateModelExcel(model_iMLirr, outputDir, dateVersion, 0);

fprintf('Done.\n');

%% Helper Functions

function model = prepareModel(modelPath, attributesPath, saveBool, savePath, excelUpdatePath)
    % prepareModel - Prepare the iML1515 model with additional attributes.

    % Check the number of input arguments
    if nargin < 3
        error('Not enough input arguments. Please provide modelPath, attributesPath, and saveBool.');
    end

    % Load the model
    model = readCbModel(modelPath);

    % Check if the model structure is valid
    if ~isstruct(model) || ~isfield(model, 'rxns')
        error('Invalid model structure. Please provide a valid COBRA model.');
    end

    % Check if the attributes file exists
    if ~isfile(attributesPath)
        error('Attributes file not found. Please provide a valid path to the attributes CSV file: %s', attributesPath);
    end

    % Map attributes from SBML/XML annotations to model fields
    % Note: The Python step embeds these into the XML using specific IDs
    model.subSystems = model.rxnisstep__46__subsystemID;
    model.location = model.rxnislocationID;
    model.mw_c =  cellfun(@str2double, model.rxnisweight_mw_cID);
    model.mw_m =  cellfun(@str2double, model.rxnisweight_mw_mID);
    model.sa_m =  cellfun(@str2double, model.rxnissurface_area_tMbID);
    model.kcat =  cellfun(@str2double, model.rxniskcatID);
    model.number_tMb = cellfun(@str2double, model.rxnisnumber_tMbID);

    % Define the value to fill NaN with
    fillValue = 100000000;
    % Use isnan and logical indexing to fill NaN values
    model.kcat(isnan(model.kcat)) = fillValue;
    
    if isfield(model, 'rxnisprot_wCoef_cID')
        model.w_c =  cellfun(@str2double, model.rxnisprot_wCoef_cID);
    end
    if isfield(model, 'rxnisprot_wCoef_mID')
        model.w_m =  cellfun(@str2double, model.rxnisprot_wCoef_mID);
    end

    [nMets,nRxns]=size(model.S); 
    if strcmp(model.rxns{nRxns},'ER_Membpool_TG_')
            fprintf('model already with both protein constraints\n')
            model.location(end-1:end) = {'Protein_constraint'};
    elseif strcmp(model.rxns{nRxns},'ER_pool_TG_')
            fprintf('model already with one protein constraint\n')
            model.location{nRxns} = 'Protein_constraint';
    end

    % Adjust MSA parameter for given reactions using external Excel
    if exist('excelUpdatePath', 'var') && exist(excelUpdatePath, 'file')
        fprintf('Updating model params from %s\n', excelUpdatePath);
        model = updateModelFromExcel(model, excelUpdatePath);
    else
        warning('MSA update Excel file not found at %s', excelUpdatePath);
    end

    % Display model information
    disp(model.modelID);
    disp(length(model.rxns));

    % Optionally save the prepared model
    if saveBool && exist('savePath', 'var')
        save(savePath, 'model');
        disp(['Model saved at: ', savePath]);
    elseif saveBool
        warning('Save path not provided. Model will not be saved.');
    end
end


function ruleListWithGenes = convertRulesWithGenes(ruleList, geneList)
    % Initialize the output cell array
    ruleListWithGenes = cell(size(ruleList));
    % Loop through each rule in the rule list
    for i = 1:numel(ruleList)
        % Get the current rule
        currentRule = ruleList{i};
        geneIndexs = regexp(ruleList{i}, 'x\((\d+)\)', 'tokens');

        % If the rule is not empty
        if ~isempty(currentRule)
            % Replace 'x(#)' with the corresponding gene ID from the gene list
            % Extract gene IDs from the rule list
            geneIDs = cellfun(@(x) geneList{str2double(x{1})}, geneIndexs, 'UniformOutput', false);
            % Join the gene IDs into a single string
            geneIDs = strjoin(geneIDs, ', ');

            % Use regular expressions to find and replace patterns
            currentRule = regexprep(currentRule, 'x\(\d+\)', geneIDs);

            % Replace '&' with 'and'
            currentRule = strrep(currentRule, '&', ' and ');

            % Replace '|' with 'or'
            currentRule = strrep(currentRule, '|', ' or ');
        end

        % Store the modified rule in the output cell array
        ruleListWithGenes{i} = currentRule;
    end
end
