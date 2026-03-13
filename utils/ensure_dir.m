function ensure_dir(targetDir)
    if ~exist(targetDir, 'dir')
        mkdir(targetDir);
    end
end
