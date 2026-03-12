% FVA workflow for the bundled MAFBA model.
% This script generates the staged MAFBA FVA dataset used by Figure 6.
% Outputs:
%   2_Basic_Analysis/4_FVA/output/raw/fva/<timestamp>/
%   data/outputs/FVA_Results/FVA_results_MAFBA.csv

baseDir = fileparts(mfilename('fullpath'));
repoRoot = fullfile(baseDir, '..', '..');
addpath(fullfile(repoRoot, 'utils'));

[model, info] = load_mafba_model();

gammaValues = [1.0, 0.3, 0.25, 0.2];
uptakeBounds = get_default_uptake_grid('fva');
atpmValue = 6.89;
analysisTypes = {'all_reactions', 'membrane', 'nonmembrane'};
storedFvaPercent = 0;

timestamp = char(datetime('now', 'Format', 'yyyy-MM-dd_HH-mm-SS'));
rawFolder = fullfile(baseDir, 'output', 'raw', 'fva', timestamp);
outputFolder = fullfile(repoRoot, 'data', 'outputs', 'FVA_Results');
ensure_dir(rawFolder);
ensure_dir(outputFolder);

excludedRxns = unique([numel(model.rxns)-6:numel(model.rxns), model.protGroup(3).rxns']);
allIdx = setdiff(1:numel(model.rxns), excludedRxns);
membraneIdx = allIdx(info.membraneMask(allIdx));
nonMembraneIdx = setdiff(allIdx, membraneIdx);

rowsPerCondition = numel(allIdx) + numel(membraneIdx) + numel(nonMembraneIdx);
totalRows = numel(gammaValues) * numel(uptakeBounds) * rowsPerCondition;

modelCol = strings(totalRows, 1);
analysisCol = strings(totalRows, 1);
reactionCol = strings(totalRows, 1);
fvaPercentCol = zeros(totalRows, 1);
gammaCol = zeros(totalRows, 1);
atpmCol = zeros(totalRows, 1);
uptakeCol = zeros(totalRows, 1);
originalCol = zeros(totalRows, 1);
minCol = zeros(totalRows, 1);
maxCol = zeros(totalRows, 1);
diffCol = zeros(totalRows, 1);

rowIdx = 1;
for g = 1:numel(gammaValues)
    gamma = gammaValues(g);
    currentModel = model;
    currentModel.ub(info.rxn.phiMmax) = info.basePhiPmax * gamma;
    currentModel.ub(info.rxn.atpm) = atpmValue;

    for u = 1:numel(uptakeBounds)
        uptakeBound = uptakeBounds(u);
        currentModel.ub(info.rxn.glc_b) = uptakeBound;

        solution = optimizeCbModel(currentModel);
        if solution.stat ~= 1
            warning('Skipping gamma %.3f uptake %.2f: infeasible optimization.', gamma, uptakeBound);
            continue;
        end

        exactModel = currentModel;
        exactModel.lb(info.rxn.biomass) = solution.f;

        analysisIdx = {allIdx, membraneIdx, nonMembraneIdx};
        for a = 1:numel(analysisTypes)
            currentIdx = analysisIdx{a};
            [minFlux, maxFlux] = fluxVariability(exactModel, 100, 'max', exactModel.rxns(currentIdx), 0, true, '2-norm');
            currentRows = rowIdx:(rowIdx + numel(currentIdx) - 1);

            modelCol(currentRows) = "MAFBA";
            analysisCol(currentRows) = string(analysisTypes{a});
            reactionCol(currentRows) = string(exactModel.rxns(currentIdx));
            fvaPercentCol(currentRows) = storedFvaPercent;
            gammaCol(currentRows) = gamma;
            atpmCol(currentRows) = atpmValue;
            uptakeCol(currentRows) = uptakeBound;
            originalCol(currentRows) = solution.x(currentIdx);
            minCol(currentRows) = minFlux;
            maxCol(currentRows) = maxFlux;
            diffCol(currentRows) = maxFlux - minFlux;
            rowIdx = rowIdx + numel(currentIdx);
        end
    end
end

validRows = modelCol ~= "";
fvaResults = table( ...
    cellstr(modelCol(validRows)), fvaPercentCol(validRows), cellstr(analysisCol(validRows)), ...
    gammaCol(validRows), atpmCol(validRows), uptakeCol(validRows), ...
    cellstr(reactionCol(validRows)), originalCol(validRows), minCol(validRows), ...
    maxCol(validRows), diffCol(validRows), ...
    'VariableNames', {'Model', 'FVA_Percent', 'Analysis_Type', 'Gamma_fMSA', ...
    'ATPM', 'Uptake_Bound', 'Reaction', 'Original_Value', 'Min_Value', ...
    'Max_Value', 'FVA_Difference'});

writetable(fvaResults, fullfile(rawFolder, 'FVA_results_MAFBA.csv'));
writetable(fvaResults, fullfile(outputFolder, 'FVA_results_MAFBA.csv'));

attrsInput = fullfile(repoRoot, 'data', 'input', 'attributes_model_ecModel_iML1515_MAFBA_g_1_v3.xlsx');
if isfile(attrsInput)
    attrTable = readtable(attrsInput, 'Sheet', 'Reactions');
    rxnAttrTable = table( ...
        (0:height(attrTable)-1)', attrTable.rxns, attrTable.rxnNames, attrTable.location, ...
        attrTable.subSystems, attrTable.mw_c, attrTable.mw_m, attrTable.number_tMb, ...
        attrTable.kcat, attrTable.grRules, strings(height(attrTable), 1), ...
        'VariableNames', {'Index', 'Rxn_ID', 'Rxn_name', 'Rxn_Loc', 'Rxn_Subsystem', ...
        'Rxn_Mw_c', 'Rxn_Mw_t', 'Rxn_Mw_p', 'Rxn_kcat', 'Rxn_genes_rules', 'Id_matchs'});
    writetable(rxnAttrTable, fullfile(outputFolder, 'iML1515_rxns_attrs.csv'));
end

fprintf('FVA analysis complete.\n');
fprintf('Raw outputs: %s\n', rawFolder);
fprintf('Staged outputs: %s\n', outputFolder);
