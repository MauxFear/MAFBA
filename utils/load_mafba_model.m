function [model, info] = load_mafba_model(varargin)
    baseDir = fileparts(mfilename('fullpath'));
    repoRoot = fullfile(baseDir, '..');
    workspaceRoot = fileparts(repoRoot);
    utilityDir = fullfile(workspaceRoot, 'Utility');
    buildDir = fullfile(repoRoot, '1_Model_Construction', '02_MATLAB_Build');
    opts = parse_inputs(varargin{:});

    addpath(baseDir);
    if isfolder(utilityDir)
        addpath(utilityDir);
    end
    if isfolder(buildDir)
        addpath(buildDir);
    end

    configure_solver();
    [model, modelSource] = build_mafba_ec_model(repoRoot, workspaceRoot, opts);

    info = struct();
    info.repoRoot = repoRoot;
    info.modelPath = modelSource.modelPath;
    info.loader = modelSource.loader;
    info.rxn.glc_f = find(strcmp(model.rxns, 'EX_glc__D_e_f'), 1);
    info.rxn.glc_b = find(strcmp(model.rxns, 'EX_glc__D_e_b'), 1);
    info.rxn.o2_f = find(strcmp(model.rxns, 'EX_o2_e_f'), 1);
    info.rxn.o2_b = find(strcmp(model.rxns, 'EX_o2_e_b'), 1);
    info.rxn.co2_f = find(strcmp(model.rxns, 'EX_co2_e_f'), 1);
    info.rxn.co2_b = find(strcmp(model.rxns, 'EX_co2_e_b'), 1);
    info.rxn.ac = find(strcmp(model.rxns, 'EX_ac_e_f'), 1);
    info.rxn.resp = find(strcmp(model.rxns, 'ATPS4rpp_f'), 1);
    info.rxn.atpm = find(strcmp(model.rxns, 'ATPM'), 1);
    info.rxn.biomass = find(model.c, 1);
    info.rxn.phiPmax = find(strcmp(model.rxns, 'Prot_Const'), 1);
    info.rxn.phiMmax = find(strcmp(model.rxns, 'Memb_Const'), 1);

    requiredRxns = [info.rxn.glc_f, info.rxn.glc_b, info.rxn.o2_f, info.rxn.o2_b, ...
        info.rxn.co2_f, info.rxn.co2_b, info.rxn.ac, info.rxn.resp, ...
        info.rxn.atpm, info.rxn.biomass, info.rxn.phiPmax, info.rxn.phiMmax];
    if any(requiredRxns == 0)
        error('One or more required reactions are missing from the bundled model.');
    end

    model = apply_default_bounds(model, info.rxn);

    info.basePhiPmax = model.ub(info.rxn.phiPmax);
    info.basePhiMmax = model.ub(info.rxn.phiMmax);
    info.membraneMask = abs(full(model.S(end, :))) > 1e-12;
end

function [model, source] = build_mafba_ec_model(repoRoot, workspaceRoot, opts)
    source = struct();
    source.loader = '';
    source.modelPath = '';

    % The prepared MAT file retains its existing build tag in the filename.
    preparedPath = fullfile(repoRoot, '1_Model_Construction', '02_MATLAB_Build', ...
        'output', 'prepared_iML1515_modelV0_AlterDfv2.mat');

    if exist('loadEcModel', 'file') == 2
        try
            gammaValues = 1;
            modelV = 0;
            dateVersion = 'AlterDfv2';  % matches the build tag of the prepared MAT file
            saveMATLAB = 0;
            saveSBML = 0;
            optimizeSigma = 1;
            learningRate = opts.alearn;
            sigmaMode = opts.sigmaMode;
            saveTag = '';

            [~, ecModels] = loadEcModel(gammaValues, modelV, dateVersion, ...
                saveMATLAB, saveSBML, optimizeSigma, learningRate, ...
                sigmaMode, saveTag);
            model = ecModels{1};
            source.loader = sprintf('loadEcModel (%s)', sigmaMode);
            source.modelPath = preparedPath;
            return;
        catch loadErr
            warning('load_mafba_model:LoadEcModelFailed', ...
                ['Falling back to direct setup after loadEcModel failed: %s'], ...
                loadErr.message);
        end
    end

    if ~isfile(preparedPath)
        error('Prepared model not found: %s', preparedPath);
    end

    preparedModel = load_prepared_model(preparedPath);

    if exist('setUpModel_v1', 'file') == 2
        sigmaParam = get_default_sigma(0, opts.sigmaMode);
        addZYprot = 1;
        gamma = 1;
        testPrint = 0;
        [~, model] = setUpModel_v1(preparedModel, opts.sigmaMode, addZYprot, ...
            sigmaParam, gamma, testPrint);
        source.loader = sprintf('prepared MAT + Utility/setUpModel_v1 (%s)', opts.sigmaMode);
        source.modelPath = preparedPath;
        return;
    end

    error(['Unable to construct the ec_model. Expected one of these loaders to be ' ...
        'available on the MATLAB path: loadEcModel or setUpModel_v1.']);
end

function model = load_prepared_model(preparedPath)
    loadedData = load(preparedPath);
    fieldNames = fieldnames(loadedData);

    model = [];
    for idx = 1:numel(fieldNames)
        candidate = loadedData.(fieldNames{idx});
        if isstruct(candidate) && isfield(candidate, 'rxns') && isfield(candidate, 'S')
            model = candidate;
            return;
        end
    end

    error('No COBRA model structure found in %s', preparedPath);
end

function sigmaParam = get_default_sigma(modelV, sigmaMode)
    switch sigmaMode
        case 'homogeneous'
            switch modelV
                case 0
                    sigmaParam = 0.000354;
                otherwise
                    error('Unsupported homogeneous sigma default for model version: %d', modelV);
            end
        case {'parametrized_totalweight', 'parametrized_TMbw', 'parametrized_MSA', 'parametrized_MSA_only'}
            switch modelV
                case 3
                    sigmaParam = 8.4701;
                case 2
                    sigmaParam = 3.2320;
                case 1
                    sigmaParam = 14.018;
                case 0
                    sigmaParam = 0.4048;
                otherwise
                    error('Unsupported model version: %d', modelV);
            end
        otherwise
            error('Unsupported sigma mode: %s', sigmaMode);
    end
end

function opts = parse_inputs(varargin)
    parser = inputParser;
    addParameter(parser, 'sigmaMode', 'parametrized_totalweight');
    addParameter(parser, 'alearn', 0.001);
    parse(parser, varargin{:});
    opts = parser.Results;
end

function model = apply_default_bounds(model, rxn)
    model.lb(rxn.o2_b) = 0;
    model.ub(rxn.o2_b) = 20;
    model.lb(rxn.o2_f) = 0;
    model.ub(rxn.o2_f) = 1000;

    model.lb(rxn.co2_b) = 0;
    model.ub(rxn.co2_b) = 0;
    model.lb(rxn.co2_f) = 0;
    model.ub(rxn.co2_f) = 1000;

    model.lb(rxn.glc_b) = 0;
    model.ub(rxn.glc_b) = 30;
    model.lb(rxn.glc_f) = 0;
    model.ub(rxn.glc_f) = 1000;

    glcdIdx = find(strcmp(model.rxns, 'GLCDpp'), 1);
    if ~isempty(glcdIdx)
        model.lb(glcdIdx) = 0;
        model.ub(glcdIdx) = 0;
    end

    ldhIdx = find(strcmp(model.rxns, 'LDH_D2'), 1);
    if ~isempty(ldhIdx)
        model.lb(ldhIdx) = 0;
        model.ub(ldhIdx) = 0;
    end
end

function configure_solver()
    solverCandidates = {'gurobi', 'ibm_cplex'};
    configured = false;

    for idx = 1:numel(solverCandidates)
        solverName = solverCandidates{idx};
        try
            changeCobraSolver(solverName, 'LP');
            configured = true;
            if strcmp(solverName, 'gurobi') && exist('changeCobraSolverParams', 'file') == 2
                changeCobraSolverParams('LP', 'feasTol', 1e-5);
                changeCobraSolverParams('LP', 'optTol', 1e-5);
            end
            break;
        catch
        end
    end

    if ~configured
        warning('No preferred LP solver was configured. Using the current COBRA solver selection.');
    end
end
