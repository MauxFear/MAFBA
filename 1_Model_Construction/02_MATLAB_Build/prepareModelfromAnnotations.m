function model = prepareModelfromAnnotations(model, attributesPath, saveBool, savePath)
% prepareModelfromAnnotations - Prepare the model with additional attributes.

% Default values
if nargin < 2
    attributesPath = false;
end
if nargin < 3
    saveBool = true;
end
if nargin < 4
    savePath = 'prepared_Model.xml';
end

% Check the number of input arguments
if nargin < 1
    error('Not enough input arguments. Please provide a modelPath.');
end


% Check if the model structure is valid
if ~isstruct(model) || ~isfield(model, 'rxns')
    error('Invalid model structure. Please provide a valid COBRA model.');
end

% Check if the attributes file exists
if ~isempty(attributesPath) && ~isfile(attributesPath)
    error('Attributes file not found. Please provide a valid path to the attributes CSV file.');
else
    % Read the Excel file as a table
    dataTable = readtable(attributesPath);
end

% Assign fields if they exist
if isfield(model, 'rxnisstep__46__subsystemID')
    model.subSystems = model.rxnisstep__46__subsystemID;
elseif attributesPath ~= false
    model.subSystems = dataTable.('Rxn_Subsystem');
end
if isfield(model, 'rxnislocationID')
    model.location = model.rxnislocationID;
elseif attributesPath ~= false
    model.subSystems = dataTable.('Rxn_Loc');
end
if isfield(model, 'rxnisweight_mw_cID')
    model.mw_c = cellfun(@str2double, model.rxnisweight_mw_cID);
elseif attributesPath ~= false
    model.subSystems = dataTable.('Rxn_mW_c');
end
if isfield(model, 'rxnisweight_mw_mID')
    model.mw_m = cellfun(@str2double, model.rxnisweight_mw_mID);
elseif attributesPath ~= false
    model.subSystems = dataTable.('Rxn_mW_m');
end
if isfield(model, 'rxniskcatID')
    model.kcat = cellfun(@str2double, model.rxniskcatID);
elseif attributesPath ~= false
    model.subSystems = dataTable.('Rxn_kcat');
end
if isfield(model, 'rxnisnumber_tMbID')
    model.number_tMb = cellfun(@str2double, model.rxnisnumber_tMbID);
elseif attributesPath ~= false
    model.subSystems = dataTable.('Rxn_numTMb');
end

% Define the value to fill NaN with
fillValue = 100000000;
% Use isnan and logical indexing to fill NaN values
model.kcat(isnan(model.kcat)) = fillValue;

if isfield(model, 'rxnisprot_wCoef_cID')
    model.w_c = cellfun(@str2double, model.rxnisprot_wCoef_cID);
    model.w_m = cellfun(@str2double, model.rxnisprot_wCoef_mID);
    % Find indices where list2 is equal to a specific string
    stringToMatch = 'Inner Membrane';
    indices = strcmp(model.location, stringToMatch);
    % Set non-selected elements to zero
    model.w_m(~indices) = 0;
end

[~, nRxns] = size(model.S); 
if strcmp(model.rxns{nRxns}, 'ER_Membpool_TG_')
    fprintf('Model already with both protein constraints\n')
    model.location(end-1:end) = {'Protein_constraint'};
elseif strcmp(model.rxns{nRxns}, 'ER_pool_TG_')
    fprintf('Model already with one protein constraint\n')
    model.location{nRxns} = 'Protein_constraint';
end

% Display model information
disp(model.modelID);
disp(length(model.rxns));

% Optionally save the prepared model
if saveBool
    save(savePath, 'model');
    disp(['Model saved at: ', savePath]);
end
end



