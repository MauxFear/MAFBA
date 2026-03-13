function newPath = constructPath(targetFolder, varargin)  
    % constructPath constructs a new path by finding the target folder within the parent directory  
    % and appending the provided subfolders. Also replaces user home directory with "~" when possible.  
    %  
    % Usage:  
    %   newPath = constructPath(targetFolder, subFolder1, subFolder2, ...)  
    %  
    % Inputs:  
    %   targetFolder - Target folder to locate within the parent directory  
    %   varargin - Additional subfolders to append to the new path  
    %  
    % Output:  
    %   newPath - Constructed new path with home directory replaced by "~" when applicable  

    % Get the current working directory  
    currentFolder = pwd;  
    
    % Initialize the new path variable  
    newPath = '';  

    % Split the current folder path into parts  
    if ispc  
        folderParts = regexp(currentFolder, '\\', 'split');  
    else  
        folderParts = regexp(currentFolder, '/', 'split');  
    end  
    
    % Remove empty elements  
    folderParts = folderParts(~cellfun('isempty', folderParts));  
    
    % Search for the target folder in the path  
    for i = 1:length(folderParts)  
        % Check if the target folder exists at this level  
        if strcmpi(targetFolder, folderParts{i})  
            % Construct base path up to the target folder  
            if ispc  
                basePath = strjoin(folderParts(1:i), '\');  
            else  
                basePath = strjoin(folderParts(1:i), '/');  
            end  
            
            % Add the target folder and any additional subfolders  
            if ~isempty(varargin)  
                parts = [{basePath} varargin];  
                if ispc  
                    newPath = strjoin(parts, '\');  
                else  
                    newPath = strjoin(parts, '/');  
                end  
            else  
                newPath = basePath;  
            end  
            
            % Debug output  
            fprintf('Found target folder at level %d\n', i);  
            break;  
        end  
    end  
    
    % Check if path was found  
    if isempty(newPath)  
        warning('Target folder "%s" not found in path: %s', targetFolder, currentFolder);  
        return;  
    end  
    
    % Ensure proper path format for OS  
    if ispc  
        newPath = strrep(newPath, '/', '\');  
        % Add drive letter if missing  
        if ~regexp(newPath, '^[A-Za-z]:\\')  
            newPath = [pwd(1:2) newPath];  
        end  
    else  
        newPath = strrep(newPath, '\', '/');  
        % Ensure absolute path  
        if ~startsWith(newPath, '/')  
            newPath = ['/' newPath];  
        end  
    end  
    
    % Get the user's home directory path (works on macOS/Linux/Windows)  
    userHome = char(java.lang.System.getProperty('user.home'));  
    
    % Convert path separators in userHome to match the OS  
    if ispc  
        userHome = strrep(userHome, '/', '\');  
    else  
        userHome = strrep(userHome, '\', '/');  
    end  
    
    % Replace user home directory with "~" if applicable  
    replacement = '~';  
    
    % Check if the path starts with the user's home directory before replacing  
    if startsWith(newPath, userHome)  
        newPath = replace(newPath, userHome, replacement);  
        fprintf('Replaced home directory with tilde: %s\n', newPath);  
    end  
    
    % Debug output  
    fprintf('Constructed path: %s\n', newPath);  
end  