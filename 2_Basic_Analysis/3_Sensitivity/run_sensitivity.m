% Sensitivity analysis using the bundled MAFBA model.
% Final staged outputs are written to:
%   data/outputs/iML1515_Sensitivity/
% Analysis-local raw intermediates are written to:
%   2_Basic_Analysis/3_Sensitivity/output/raw/sensitivity/<timestamp>/

baseDir = fileparts(mfilename('fullpath'));
repoRoot = fullfile(baseDir, '..', '..');
addpath(fullfile(repoRoot, 'utils'));

[model, info] = load_mafba_model();

gammaValues = [1.0, 0.3, 0.25, 0.2, 0.15];
changeFactors = [1, 2, 5, 10, 20, 50, 1/2, 1/5, 1/10, 1/20, 1/50];
uptakeBounds = get_default_uptake_grid('sweep');
iterationValue = 1;
randomSeedValue = 1;
atpmValue = 6.89;

timestamp = char(datetime('now', 'Format', 'yyyy-MM-dd_HH-mm-SS'));
rawFolder = fullfile(baseDir, 'output', 'raw', 'sensitivity', timestamp);
outputFolder = fullfile(repoRoot, 'data', 'outputs', 'iML1515_Sensitivity');
ensure_dir(rawFolder);
ensure_dir(outputFolder);

originalCostC = full(model.S(end-1, :));
originalCostM = full(model.S(end, :));

excludedRxns = unique([numel(model.rxns)-6:numel(model.rxns), model.protGroup(3).rxns']);
rxnsToScale = setdiff(1:numel(model.rxns), excludedRxns);

numRows = numel(gammaValues) * numel(changeFactors) * numel(uptakeBounds);
numModelRows = numel(gammaValues) * numel(changeFactors);

fluxMatrix = zeros(numRows, numel(model.rxns) + 5);
sectorMatrix = zeros(numRows, 21);
aceMatrix = zeros(numModelRows, 6);
modelMatrix = zeros(numModelRows, 4 + 3 * numel(model.rxns));

rowIdx = 1;
modelRowIdx = 1;

for g = 1:numel(gammaValues)
    gamma = gammaValues(g);
    gammaModel = model;
    gammaModel.ub(info.rxn.phiMmax) = info.basePhiPmax * gamma;
    gammaModel.ub(info.rxn.atpm) = atpmValue;

    for cf = 1:numel(changeFactors)
        changeFactor = changeFactors(cf);
        currentModel = gammaModel;

        currentModel.S(end-1, :) = originalCostC;
        currentModel.S(end, :) = originalCostM;
        currentModel.S(end-1, rxnsToScale) = originalCostC(rxnsToScale) * changeFactor;
        currentModel.S(end, rxnsToScale) = originalCostM(rxnsToScale) * changeFactor;

        scalingVector = ones(1, numel(model.rxns));
        scalingVector(rxnsToScale) = changeFactor;
        modelMatrix(modelRowIdx, :) = [ ...
            gamma, iterationValue, changeFactor, randomSeedValue, ...
            scalingVector, full(currentModel.S(end-1, :)), full(currentModel.S(end, :))];

        biomassValues = zeros(numel(uptakeBounds), 1);
        acetateValues = zeros(numel(uptakeBounds), 1);
        respirationValues = zeros(numel(uptakeBounds), 1);

        for i = 1:numel(uptakeBounds)
            uptakeBound = uptakeBounds(i);
            currentModel.ub(info.rxn.glc_b) = uptakeBound;

            solution = get_pFBAsolution_ec_model(currentModel);

            fluxMatrix(rowIdx, :) = [ ...
                gamma, iterationValue, changeFactor, uptakeBound, ...
                solution.x(:)', randomSeedValue];
            sectorMatrix(rowIdx, :) = [ ...
                gamma, iterationValue, changeFactor, uptakeBound, solution.f, ...
                solution.x(info.rxn.atpm), randomSeedValue, solution.phiCm, ...
                solution.phiCc, solution.phiR, solution.phiEc, solution.phiEr, ...
                solution.phiEm, solution.phiP_m, solution.phiP, solution.phiM, ...
                solution.phiCm_m, solution.phiEr_m, solution.phiEm_m, ...
                solution.phiPmax, solution.phiMmax];

            biomassValues(i) = solution.x(info.rxn.biomass);
            acetateValues(i) = solution.x(info.rxn.ac);
            respirationValues(i) = solution.x(info.rxn.resp);
            rowIdx = rowIdx + 1;
        end

        aceMatrix(modelRowIdx, :) = [ ...
            gamma, iterationValue, changeFactor, ...
            detect_acetate_threshold(biomassValues, acetateValues, respirationValues), ...
            atpmValue, randomSeedValue];
        modelRowIdx = modelRowIdx + 1;
    end
end

fluxTable = array2table(fluxMatrix, ...
    'VariableNames', ["Gamma", "Iteration", "Change_Factor", "Uptake_Bound", ...
    strcat("rx_", string(model.rxns')), "Random_Seed"]);
sectorsTable = array2table(sectorMatrix, ...
    'VariableNames', ["Gamma", "Iteration", "Change_Factor", "Uptake_Bound", ...
    "Growth", "ATPM", "Random_Seed", "phiCm", "phiCc", "phiR", "phiEc", ...
    "phiEr", "phiEm", "phiP_m", "phiP", "phiM", "phiCm_m", "phiEr_m", ...
    "phiEm_m", "phiPmax", "phiMmax"]);
aceTable = array2table(aceMatrix, ...
    'VariableNames', ["Gamma", "Iteration", "Change_Factor", ...
    "Acetate_threshold", "ATPM", "Random_Seed"]);
modelTable = array2table(modelMatrix, ...
    'VariableNames', ["Gamma", "Iteration", "Change_Factor", "Random_Seed", ...
    strcat("rx_", string(model.rxns'), "_scale"), ...
    strcat("rx_", string(model.rxns'), "_costC"), ...
    strcat("rx_", string(model.rxns'), "_costM")]);

writetable(fluxTable, fullfile(rawFolder, 'dataFlux_Sensitivity.csv'));
writetable(sectorsTable, fullfile(rawFolder, 'sectors_Sensitivity.csv'));
writetable(aceTable, fullfile(rawFolder, 'aceT_Sensitivity.csv'));
writetable(modelTable, fullfile(rawFolder, 'model_Sensitivity.csv'));

writetable(fluxTable, fullfile(outputFolder, 'dataFlux_Sensitivity.csv'));
writetable(sectorsTable, fullfile(outputFolder, 'sectors_Sensitivity.csv'));
writetable(aceTable, fullfile(outputFolder, 'aceT_Sensitivity.csv'));

fprintf('Sensitivity analysis complete.\n');
fprintf('Raw outputs: %s\n', rawFolder);
fprintf('Staged outputs: %s\n', outputFolder);
