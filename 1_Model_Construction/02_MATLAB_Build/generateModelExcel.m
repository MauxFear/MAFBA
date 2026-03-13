function [structuredData, fullFilePath] = generateModelExcel_v6(modelOrPath, outputFolder, version, costBool)
% generateModelExcel - Generate Excel spreadsheets from a COBRA model structure.
%
% Input:
%   - modelOrPath: Path to the COBRA model file (MAT or XML format) or a COBRA model structure.
%   - outputFolder: Path to the folder where the Excel file will be saved.
%   - version: String defining the version of the file generated. 
    if nargin < 4
        costBool = 0;
    end

% Load the model if a path is provided
    if ischar(modelOrPath)
        [~, ~, ext] = fileparts(modelOrPath);
        ext = lower(ext);
        switch ext
            case '.mat'
                model = load(modelOrPath);
                model = model.model; % assuming the variable in the MAT file is named 'model'
            case '.xml'
                model = readCbModel(modelOrPath);
            otherwise
                error('Unsupported file format. Please provide a .mat or .xml file.');
        end
    elseif isstruct(modelOrPath) && isfield(modelOrPath, 'rxns')
        % Assume modelOrPath is a valid model structure
        model = modelOrPath;
    else
        error('Invalid input. Please provide either a model structure or a file path.');
    end
    
    if iscell(model.kcat)
        model.kcat = cellfun(@(x) x(1), model.kcat);
    end

    if iscell(model.mw_c)
        model.mw_c = cellfun(@(x) x(1), model.mw_c);
    end

    if ~isfield(model, 'w_c') && isfield(model, 'w')
        model.w_c = model.w;
    end
    
    % Create a cell array to store protein group for each reaction
    reactionProtGroups = repmat({''}, numel(model.rxns), 1);
    if isfield(model, 'protGroup') && isstruct(model.protGroup)
        % Loop through the struct array
        for i = 1:numel(model.protGroup)
            % Check if the rxns field is not empty
            if ~isempty(model.protGroup(i).rxns)
                groupName = model.protGroup(i).name;
                rxnIndices = model.protGroup(i).rxns;
                
                % Assign group name to corresponding reactions
                for j = 1:numel(rxnIndices)
                    idx = rxnIndices(j);
                    if idx > 0 && idx <= numel(reactionProtGroups)
                        reactionProtGroups{idx} = groupName;
                    end
                end
            end
        end
    end
    model.reactionProtGroups = reactionProtGroups;

    % Convert the lists to strings and store them in the new fields
    newFieldLocs = cell(1, numel(model.rxns));
    newFieldMws = cell(1, numel(model.rxns));

    % Loop over reactions in model.rxns
    for r = 1:numel(model.rxns)
        if model.grRules{r}
            % Extract gene IDs from the grRules field
            geneIDs = regexp(model.grRules{r}, 'b\d+', 'match');

            % Preallocate arrays
            locList = cell(1, numel(geneIDs));
            mwList = cell(1, numel(geneIDs));

            % Iterate over gene IDs
            for j = 1:numel(geneIDs)
                % Find the index of the gene in model.genes using logical indexing
                geneIndex = find(strcmp(model.genes, geneIDs{j}), 1);

                % Check if the gene index is valid
                if ~isempty(geneIndex)
                    % Extract information from model.genes
                    loc = model.geneisstepdb__46__locationID{geneIndex};
                    mw = model.geneisstepdb__46__weightID{geneIndex};

                    % Append the information to the lists
                    locList{j} = loc;
                    % Split the string into a cell array of substrings using ';' as the delimiter
                    strParts = strsplit(mw, '; ');

                    % Concatenate the numbers into a single string without spaces
                    concatenatedString = strjoin(strParts, '');

                    mwList{j} = concatenatedString;
                end
            end
            % Convert the lists to strings and store them in the new fields
            newFieldLocs{r} = strjoin(locList, '; ');
            newFieldMws{r} = strjoin(mwList, ';');
        end
    end
    % Add new fields to the structure
    model.genesLocs = newFieldLocs';
    model.genesMws = newFieldMws';

% --- Verification Step ---
fprintf('\n--- Verifying first 20 reaction protein groups ---\n');
for i = 1:min(20, numel(model.rxns))
    fprintf('Reaction: %-20s, Group: %s\n', model.rxns{i}, model.reactionProtGroups{i});
end
fprintf('--- Verification complete ---\n\n');

%     Get reactionsCosts
    if costBool
        gamma_val = 1 ;
        % Reactions indices
        phi_max_r = find(strcmp(model.rxns,'Prot_Const'));  % phi_max -- Overall Proteome Allocation
        phi_Mmax_r = find(strcmp(model.rxns,'Memb_Const')); % phi_Mmax -- Membrane Allocation 
        phi_max_value = model.ub(phi_max_r); % phi_max fixed value 

        % Use gamma to calculate the phi_Mmax value
        phi_Mmax_value = phi_max_value*gamma_val;
        % Adjust the Membrane allocation reaction boundary
        model.ub(phi_Mmax_r) = phi_Mmax_value; 
        % Get the costs by optimizing the Cobra Model
        [costsP, costsM, solution] = getReactionCosts(model);


        model.costsP = costsP;
        model.costsM = costsM;
        model.fluxesFBA = solution.x;
        
        fieldsToKeep1 = {'rxns', 'rxnNames', 'kcat','subSystems', 'location', ...
        'mw_c', 'mw_m', 'sa_m', 'number_tMb', 'rules', 'grRules', 'reactionProtGroups', 'lb', 'ub', 'c', 'w_c','w_m', 'genesLocs',...
        'genesMws', 'costsP','costsM','fluxesFBA'};
        
    else
        fieldsToKeep1 = {'rxns', 'rxnNames', 'kcat', 'subSystems', 'location', ...
        'mw_c', 'mw_m', 'sa_m', 'number_tMb','rules', 'grRules', 'reactionProtGroups', 'lb', 'ub', 'c', 'w_c','w_m', 'genesLocs',...
        'genesMws'};
    end
    
% Fields for table 1
structM1 = struct();
for i = 1:numel(fieldsToKeep1)
    fieldName = fieldsToKeep1{i};
    if isfield(model, fieldName)
        structM1.(fieldName) = model.(fieldName);
    end
end

% Ensure that each field has the length of model.rxns
fieldskeeped = fieldnames(structM1);
for i = 1:numel(fieldskeeped)
    nField = numel(structM1.(fieldskeeped{i}));
    nRxns = numel(structM1.rxns);

    if isfield(model, fieldskeeped{i}) || nField ~= nRxns
        % If the field is missing or has a different length, pad with nan values
        difNum = nRxns - nField;
        if iscell (structM1.(fieldskeeped{i}))
            structM1.(fieldskeeped{i})(end+1:end+difNum) = {nan(1, difNum)};
        else
            structM1.(fieldskeeped{i})(end+1:end+difNum) = nan(1, difNum);
        end
    end
end
T1 = struct2table(structM1);

% Fields for table 2
fieldsToKeep2 = {'mets', 'metNames', 'metFormulas', 'metCharges', 'csense', 'b'};
structM2 = rmfield(model, setdiff(fieldnames(model), fieldsToKeep2));
T2 = struct2table(structM2);

% Fields for table 3
fieldsToKeep3 = {'genes', 'geneNames', 'geneEntrezID', 'geneEcoGeneID', 'proteins', 'proteinisuniprotID'};
structM3 = rmfield(model, setdiff(fieldnames(model), fieldsToKeep3));
T3 = struct2table(structM3);

% Create Excel file
filename = sprintf('attributes_model_%s_%s.xlsx', model.modelID, version);
fullFilePath = fullfile(outputFolder, filename);

% Ensure output folder exists
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

% Tag the file if it already exists
if isfile(fullFilePath)
    tagExistingFile(fullFilePath);
end

% Write each table to a separate sheet
writetable(T1, fullFilePath, 'Sheet', 'Reactions');
writetable(T2, fullFilePath, 'Sheet', 'Mets');
writetable(T3, fullFilePath, 'Sheet', 'Genes');

disp('Excel file generated successfully.');
disp(['File saved at: ', fullFilePath]);

structuredData= struct();
structuredData.rxnsTable = T1;
structuredData.metsTable = T2;
structuredData.genesTable = T3;

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
