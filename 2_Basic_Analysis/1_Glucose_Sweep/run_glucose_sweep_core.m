function run_glucose_sweep_core(varargin)
% Shared glucose sweep runner for homogeneous and parametrized analyses.

    opts = parse_inputs(varargin{:});

    baseDir = fileparts(mfilename('fullpath'));
    repoRoot = fullfile(baseDir, '..', '..');
    addpath(fullfile(repoRoot, 'utils'));

    [model, info] = load_mafba_model( ...
        'sigmaMode', opts.sigmaMode, ...
        'alearn', opts.alearn);

    gammaValues = [1.0, 0.3, 0.25, 0.2];
    acetateThreshold = 0.01;
    uptakeBounds = generate_canonical_glucose_grid();

    timestamp = char(datetime('now', 'Format', 'yyyy-MM-dd_HH-mm-SS'));
    rawFolder = fullfile(baseDir, 'output', 'raw', opts.rawFolderName, timestamp);
    logFolder = fullfile(rawFolder, 'logs');
    outputFolder = fullfile(repoRoot, 'data', 'outputs', opts.outputFolderName);
    ensure_dir(rawFolder);
    ensure_dir(logFolder);
    ensure_dir(outputFolder);

    logFile = fullfile(logFolder, 'uptake_grid.log');
    metadataFile = fullfile(rawFolder, 'run_metadata.json');
    fluxRawFile = fullfile(rawFolder, 'flux_data.csv');
    sectorsRawFile = fullfile(rawFolder, 'sectors_data.csv');
    aceRawFile = fullfile(rawFolder, 'aceL_data.csv');
    kcatsRawFile = fullfile(rawFolder, 'kcats_data.csv');
    coeffsRawFile = fullfile(rawFolder, 'coeffs_data.csv');

    initialize_csv(fluxRawFile, create_flux_raw_template(model));
    initialize_csv(sectorsRawFile, create_sectors_raw_template());
    initialize_csv(aceRawFile, create_acetate_raw_template());
    initialize_csv(kcatsRawFile, create_summary_template());
    initialize_csv(coeffsRawFile, create_coeffs_template(model));

    numRows = numel(gammaValues) * numel(uptakeBounds);
    stagedFluxMatrix = zeros(numRows, numel(model.rxns) + 2);
    stagedSectorMatrix = zeros(numRows, 17);
    stagedAceMatrix = NaN(numel(gammaValues), 2);
    stagedRowIdx = 1;

    membraneCoeffs = full(model.S(end, :));

    write_log(logFile, sprintf('Starting %s glucose sweep with %d gamma values and %d uptake points.', ...
        opts.modeLabel, numel(gammaValues), numel(uptakeBounds)));
    write_log(logFile, sprintf('Loader: %s', info.loader));
    write_log(logFile, sprintf('Model path: %s', info.modelPath));
    write_log(logFile, sprintf('Sigma mode: %s', opts.sigmaMode));

    for g = 1:numel(gammaValues)
        gamma = gammaValues(g);
        scenarioModel = model;
        scenarioModel.ub(info.rxn.phiMmax) = info.basePhiPmax * gamma;

        write_log(logFile, sprintf('Gamma %.3f: phiMmax bound = %.6f', ...
            gamma, scenarioModel.ub(info.rxn.phiMmax)));

        lambdaRecorded = false;
        lambdaMu = NaN;
        lambdaGlucose = NaN;
        lambdaAcetate = NaN;
        phiMAtLambda = NaN;
        muMax = -Inf;

        for i = 1:numel(uptakeBounds)
            glucoseUptake = uptakeBounds(i);
            scenarioModel.lb(info.rxn.glc_b) = 0;
            scenarioModel.ub(info.rxn.glc_b) = glucoseUptake;

            solution = get_pFBAsolution_ec_model(scenarioModel);
            if solution.stat == 1
                objectiveStatus = "optimal";
                mu = solution.x(info.rxn.biomass);
                acetateFlux = solution.x(info.rxn.ac);
                muMax = max(muMax, mu);
            else
                objectiveStatus = "failed";
                mu = NaN;
                acetateFlux = NaN;
            end

            isLambdaAc = false;
            if ~lambdaRecorded && solution.stat == 1 && acetateFlux >= acetateThreshold
                lambdaRecorded = true;
                isLambdaAc = true;
                lambdaMu = mu;
                lambdaGlucose = glucoseUptake;
                lambdaAcetate = acetateFlux;
                phiMAtLambda = solution.phiM;

                aceRow = table(gamma, lambdaMu, lambdaGlucose, lambdaAcetate, phiMAtLambda, ...
                    "Acetate onset detected", ...
                    'VariableNames', {'Gamma', 'lambda_ac_mu', 'lambda_ac_glucose', ...
                    'lambda_ac_acetate', 'phi_M_at_lambda', 'notes'});
                append_table(aceRawFile, aceRow);
                write_log(logFile, sprintf('Lambda_ac detected at gamma %.3f, mu %.6f, glucose %.6f, acetate %.6f.', ...
                    gamma, lambdaMu, lambdaGlucose, lambdaAcetate));
            end

            fluxRow = build_flux_raw_row(model, gamma, glucoseUptake, mu, objectiveStatus, isLambdaAc, solution);
            sectorsRow = build_sectors_raw_row(gamma, glucoseUptake, mu, objectiveStatus, isLambdaAc, solution);
            append_table(fluxRawFile, fluxRow);
            append_table(sectorsRawFile, sectorsRow);

            if solution.stat == 1
                stagedFluxMatrix(stagedRowIdx, :) = [gamma, glucoseUptake, solution.x(:)'];
                stagedSectorMatrix(stagedRowIdx, :) = [ ...
                    gamma, glucoseUptake, mu, solution.phiCm, solution.phiCc, ...
                    solution.phiR, solution.phiEc, solution.phiEr, solution.phiEm, ...
                    solution.phiP_m, solution.phiP, solution.phiM, solution.phiCm_m, ...
                    solution.phiEr_m, solution.phiEm_m, solution.phiPmax, solution.phiMmax];
            else
                stagedFluxMatrix(stagedRowIdx, :) = [gamma, glucoseUptake, NaN(1, numel(model.rxns))];
                stagedSectorMatrix(stagedRowIdx, :) = [gamma, glucoseUptake, NaN(1, 15)];
            end
            stagedRowIdx = stagedRowIdx + 1;
        end

        if lambdaRecorded
            stagedAceMatrix(g, :) = [gamma, lambdaMu];
        else
            stagedAceMatrix(g, :) = [gamma, NaN];
            write_log(logFile, sprintf('No acetate onset detected for gamma %.3f.', gamma));
        end

        summaryRow = table(gamma, muMax, scenarioModel.ub(info.rxn.phiMmax), info.basePhiPmax, ...
            lambdaMu, lambdaAcetate, string(opts.summaryNote), ...
            'VariableNames', {'Gamma', 'mu_max', 'phi_Mmax_bound', 'phi_Pmax_bound', ...
            'lambda_ac_mu', 'lambda_ac_acetate', 'notes'});
        append_table(kcatsRawFile, summaryRow);

        coeffsRow = array2table([gamma, membraneCoeffs], ...
            'VariableNames', ["Gamma", string(model.rxns')]);
        append_table(coeffsRawFile, coeffsRow);
    end

    stagedFluxTable = array2table(stagedFluxMatrix, ...
        'VariableNames', ["Gamma", "Uptake_Bound", string(model.rxns')]);
    stagedSectorTable = array2table(stagedSectorMatrix, ...
        'VariableNames', ["Gamma", "Uptake_Bound", "Growth", "phiCm", "phiCc", "phiR", ...
        "phiEc", "phiEr", "phiEm", "phiP_m", "phiP", "phiM", "phiCm_m", ...
        "phiEr_m", "phiEm_m", "phiPmax", "phiMmax"]);
    stagedAceTable = array2table(stagedAceMatrix, ...
        'VariableNames', ["Gamma", "Acetate_threshold"]);

    writetable(stagedFluxTable, fullfile(outputFolder, opts.fluxFileName));
    writetable(stagedSectorTable, fullfile(outputFolder, opts.sectorsFileName));
    writetable(stagedAceTable, fullfile(outputFolder, opts.aceFileName));

    metadata = struct();
    metadata.timestamp = timestamp;
    metadata.model_path = info.modelPath;
    metadata.loader = info.loader;
    metadata.sigma_mode = opts.sigmaMode;
    metadata.output_folder = outputFolder;
    metadata.gamma_values = gammaValues;
    metadata.acetate_threshold = acetateThreshold;
    metadata.glucose_grid = uptakeBounds;
    metadata.n_flux_rows = height(stagedFluxTable);
    metadata.n_gamma_values = numel(gammaValues);
    jsonText = jsonencode(metadata, 'PrettyPrint', true);
    fid = fopen(metadataFile, 'w');
    fprintf(fid, '%s', jsonText);
    fclose(fid);

    write_log(logFile, sprintf('Completed %s glucose sweep. Raw outputs: %s', opts.modeLabel, rawFolder));
    write_log(logFile, sprintf('Staged outputs: %s', outputFolder));

    fprintf('%s glucose sweep complete.\n', opts.modeLabel);
    fprintf('Loader: %s\n', info.loader);
    fprintf('Model path: %s\n', info.modelPath);
    fprintf('Raw outputs: %s\n', rawFolder);
    fprintf('Staged outputs: %s\n', outputFolder);
end

function opts = parse_inputs(varargin)
    parser = inputParser;
    addParameter(parser, 'sigmaMode', 'parametrized_totalweight');
    addParameter(parser, 'modeLabel', 'Parametrized');
    addParameter(parser, 'rawFolderName', 'glucose_sweep_parametrized');
    addParameter(parser, 'outputFolderName', 'Parametrized');
    addParameter(parser, 'fluxFileName', 'data_modelV0_kcatscomb_Parametrized.csv');
    addParameter(parser, 'sectorsFileName', 'sectors_modelV0_kcatscomb_Parametrized.csv');
    addParameter(parser, 'aceFileName', 'aceTData_modelV0_kcatscomb_Parametrized.csv');
    addParameter(parser, 'summaryNote', 'Parametrized glucose sweep');
    addParameter(parser, 'alearn', 0.001);
    parse(parser, varargin{:});
    opts = parser.Results;
end

function uptakeBounds = generate_canonical_glucose_grid()
    coarse = linspace(0.2, 9.8, 15);
    dense = 10 + 25 * (linspace(0, 1, 25) .^ 2);
    uptakeBounds = [coarse, dense];
end

function template = create_flux_raw_template(model)
    template = table('Size', [0, numel(model.rxns) + 5], ...
        'VariableTypes', [{'double', 'double', 'double', 'string', 'logical'}, ...
        repmat({'double'}, 1, numel(model.rxns))], ...
        'VariableNames', ["Gamma", "glc_uptake", "mu", "objective_status", ...
        "is_lambda_ac", string(model.rxns')]);
end

function template = create_sectors_raw_template()
    template = table('Size', [0, 13], ...
        'VariableTypes', {'double', 'double', 'double', 'string', 'logical', ...
        'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double'}, ...
        'VariableNames', {'Gamma', 'glc_uptake', 'mu', 'objective_status', 'is_lambda_ac', ...
        'phiEc', 'phiEm', 'phiEr', 'phiCm', 'phiR', 'phiCc', 'unassigned', 'phiM'});
end

function template = create_acetate_raw_template()
    template = table('Size', [0, 6], ...
        'VariableTypes', {'double', 'double', 'double', 'double', 'double', 'string'}, ...
        'VariableNames', {'Gamma', 'lambda_ac_mu', 'lambda_ac_glucose', ...
        'lambda_ac_acetate', 'phi_M_at_lambda', 'notes'});
end

function template = create_summary_template()
    template = table('Size', [0, 7], ...
        'VariableTypes', {'double', 'double', 'double', 'double', 'double', 'double', 'string'}, ...
        'VariableNames', {'Gamma', 'mu_max', 'phi_Mmax_bound', 'phi_Pmax_bound', ...
        'lambda_ac_mu', 'lambda_ac_acetate', 'notes'});
end

function template = create_coeffs_template(model)
    template = array2table(zeros(0, numel(model.rxns) + 1), ...
        'VariableNames', ["Gamma", string(model.rxns')]);
end

function row = build_flux_raw_row(model, gamma, glucoseUptake, mu, objectiveStatus, isLambdaAc, solution)
    if solution.stat == 1
        fluxValues = solution.x(:)';
    else
        fluxValues = NaN(1, numel(model.rxns));
    end
    row = array2table([gamma, glucoseUptake, mu, double(isLambdaAc), fluxValues], ...
        'VariableNames', ["Gamma", "glc_uptake", "mu", "is_lambda_ac", string(model.rxns')]);
    row = addvars(row, string(objectiveStatus), 'After', 'mu', 'NewVariableNames', 'objective_status');
    row.is_lambda_ac = logical(row.is_lambda_ac);
end

function row = build_sectors_raw_row(gamma, glucoseUptake, mu, objectiveStatus, isLambdaAc, solution)
    if solution.stat == 1
        unassigned = solution.phiPmax - solution.phiP;
        row = table(gamma, glucoseUptake, mu, string(objectiveStatus), logical(isLambdaAc), ...
            solution.phiEc, solution.phiEm, solution.phiEr, solution.phiCm, ...
            solution.phiR, solution.phiCc, unassigned, solution.phiM, ...
            'VariableNames', {'Gamma', 'glc_uptake', 'mu', 'objective_status', ...
            'is_lambda_ac', 'phiEc', 'phiEm', 'phiEr', 'phiCm', 'phiR', ...
            'phiCc', 'unassigned', 'phiM'});
    else
        row = table(gamma, glucoseUptake, mu, string(objectiveStatus), logical(isLambdaAc), ...
            NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, ...
            'VariableNames', {'Gamma', 'glc_uptake', 'mu', 'objective_status', ...
            'is_lambda_ac', 'phiEc', 'phiEm', 'phiEr', 'phiCm', 'phiR', ...
            'phiCc', 'unassigned', 'phiM'});
    end
end

function initialize_csv(filePath, template)
    writetable(template, filePath);
end

function append_table(filePath, row)
    writetable(row, filePath, 'WriteMode', 'append');
end

function write_log(logFile, message)
    fid = fopen(logFile, 'a');
    fprintf(fid, '[%s] %s\n', char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')), message);
    fclose(fid);
end
