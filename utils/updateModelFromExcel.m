function model = updateModelFromExcel(model, excelFile)
    % Read the Excel file
    excelData = readtable(excelFile);
    
    % Iterate over rows
    for i = 1:size(excelData, 1)
        reaction_id = excelData{i, 1}{1}; % the first column contains reaction IDs
        field_id = excelData{i, 2}{1}; % the second column contains field IDs
        field_value = excelData{i, 3}; % the third column contains field values
        
        % Find the index of the reaction_id in the model
        idx = find(strcmp(model.rxns, reaction_id));
        if isfield(model,field_id)
            % Update field value in the model
            if ~isempty(idx)
                model.(field_id)(idx) = field_value;
                fprintf('Reaction ID: %s was updated with a sa_m value of %d \n',reaction_id, field_value);
            else
                fprintf('Reaction ID %s was not found in the model. Skipping...\n', reaction_id);
            end
        else
            fprintf('Field ID %s was not found in the model. Skipping...\n', field_id);
        end
    end
end
