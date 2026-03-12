%%
% Step 2: Prepare the annotated model in MATLAB
saveBool = true;
modelV = 1;
dateVersion = 'MAFBA';
modelID = 'iML1515';
modelCode = sprintf('%s_%s', modelID, dateVersion);

baseDir = fileparts(mfilename('fullpath'));
repoRoot = fullfile(baseDir, '..', '..');
outputRoot = fullfile(repoRoot, '1_Model_Construction', '03_MAFBA_Models_Outputs', modelCode);

modelFileName = sprintf('annotated_%s_irreversible_V%d_%s.xml', modelID, modelV, dateVersion);
attributesFileName = sprintf('%s_irreversible_reactionDataV%d_%s.csv', modelID, modelV, dateVersion);

modelPath = fullfile(outputRoot, 'model_OutputData', modelFileName);
attributesPath = fullfile(outputRoot, 'models_InputData', attributesFileName);

outputDir = fullfile(baseDir, 'output');
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end
savePath = fullfile(outputDir, sprintf('prepared_%s_modelV%d_%s.mat', modelID, modelV, dateVersion));

attributesOut = fullfile(outputRoot, 'attributes_spreadsheets');
if ~exist(attributesOut, 'dir')
    mkdir(attributesOut);
end

updateExcel = ''; % optional; set to the update Excel path before running if needed

addpath(fullfile(repoRoot, 'utils'));

model_iMLirr = prepareModel(modelPath, attributesPath, saveBool, savePath, updateExcel);

grRules = convertRulesWithGenes(model_iMLirr.rules, model_iMLirr.genes);
model_iMLirr.grRules = grRules;

generateModelExcel(model_iMLirr, attributesOut, dateVersion, 0);

function model = prepareModel(modelPath, attributesPath, saveBool, savePath, excelUpdatePath)
% prepareModel - Prepare the iML1515 model with additional attributes.

if nargin < 4
    error('Not enough input arguments. Provide modelPath, attributesPath, saveBool, savePath, excelUpdatePath.');
end

model = readCbModel(modelPath);

if ~isstruct(model) || ~isfield(model, 'rxns')
    error('Invalid model structure. Please provide a valid COBRA model.');
end

if ~isfile(attributesPath)
    error('Attributes file not found: %s', attributesPath);
end

model.subSystems = model.rxnisstep__46__subsystemID;
model.location = model.rxnislocationID;
model.mw_c =  cellfun(@str2double, model.rxnisweight_mw_cID);
model.mw_m =  cellfun(@str2double, model.rxnisweight_mw_mID);
model.sa_m =  cellfun(@str2double, model.rxnissurface_area_tMbID);
model.kcat =  cellfun(@str2double, model.rxniskcatID);
model.number_tMb = cellfun(@str2double, model.rxnisnumber_tMbID);

fillValue = 100000000;
model.kcat(isnan(model.kcat)) = fillValue;

if isfield(model, 'rxnisprot_wCoef_cID')
    model.w_c =  cellfun(@str2double, model.rxnisprot_wCoef_cID);
end
if isfield(model, 'rxnisprot_wCoef_mID')
    model.w_m =  cellfun(@str2double, model.rxnisprot_wCoef_mID);
end

[nMets, nRxns] = size(model.S);
if strcmp(model.rxns{nRxns}, 'ER_Membpool_TG_')
    fprintf('model already with both protein constraints\n');
    model.location(end-1:end) = {'Protein_constraint'};
elseif strcmp(model.rxns{nRxns}, 'ER_pool_TG_')
    fprintf('model already with one protein constraint\n');
    model.location{nRxns} = 'Protein_constraint';
end

if ~isempty(excelUpdatePath) && isfile(excelUpdatePath)
    model = updateModelFromExcel(model, excelUpdatePath);
else
    warning('No update Excel provided. Skipping updateModelFromExcel.');
end

if saveBool
    save(savePath, 'model');
    disp(['Model saved at: ', savePath]);
end
end

function ruleListWithGenes = convertRulesWithGenes(ruleList, geneList)
    ruleListWithGenes = cell(size(ruleList));
    for i = 1:numel(ruleList)
        currentRule = ruleList{i};
        geneIndexs = regexp(ruleList{i}, 'x\((\d+)\)', 'tokens');

        if ~isempty(currentRule)
            geneIDs = cellfun(@(x) geneList{str2double(x{1})}, geneIndexs, 'UniformOutput', false);
            geneIDs = strjoin(geneIDs, ', ');
            currentRule = regexprep(currentRule, 'x\(\d+\)', geneIDs);
            currentRule = strrep(currentRule, '&', ' and ');
            currentRule = strrep(currentRule, '|', ' or ');
        end
        ruleListWithGenes{i} = currentRule;
    end
end
