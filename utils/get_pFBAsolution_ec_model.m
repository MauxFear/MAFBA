function solution = get_pFBAsolution_ec_model(model)
    % Solve parsimonious FBA for MAFBA ec-models
    % First maximizes growth, then minimizes protein cost while maintaining growth
    
    FBAsoln = optimizeCbModel(model);
    if FBAsoln.stat == 1
        model.lb(model.c==1) = FBAsoln.f;
        model_min = changeObjective(model, 'Prot_Const');
        % minimize the flux measuring demand (netFlux)
        MinimizedFlux = optimizeCbModel(model_min,'min');
    %     [GeneClasses, RxnClasses, modelIrrevFM, MinimizedFlux] = pFBA(model, 'skipclass', 1);
        solution = MinimizedFlux;
    else
       solution =  FBAsoln;
    end
    solution.phiPmax = model.ub(end-1);
    solution.phiMmax = model.ub(end);
    
    if solution.stat == 1
        biomassRxn = find(model.c);
        solution.f = solution.x(biomassRxn);
%         fprintf('The maximum growth rate is %1.4f \n', solution.f);
        solution.x = solution.x(1:length(model.rxns));
        Mrxns = [model.protGroup(1).rxns; model.protGroup(4).rxns; model.protGroup(5).rxns];
        coefficients_Pc = full(model.S(end-1,:))';
        coefficients_Pm = full(model.S(end,:))';
        %Adding the sum of each sector in the overall proteome
        costP_m = coefficients_Pc(Mrxns).* abs(solution.x(Mrxns));
        solution.phiP_m = sum(costP_m); % phiP_m
        costP_c = coefficients_Pc(model.protGroup(1).rxns).* abs(solution.x(model.protGroup(1).rxns));
        solution.phiCm = sum(costP_c); %phiCm
        costP_r = coefficients_Pc(model.protGroup(3).rxns).* abs(solution.x(model.protGroup(3).rxns));
        solution.phiR = sum(costP_r); %phiR
        costP_ec = coefficients_Pc(model.protGroup(2).rxns).* abs(solution.x(model.protGroup(2).rxns));
        solution.phiEc = sum(costP_ec); %phiEc
        costP_er = coefficients_Pc(model.protGroup(4).rxns).* abs(solution.x(model.protGroup(4).rxns));
        solution.phiEr = sum(costP_er); %phiEr
        costP_em = coefficients_Pc(model.protGroup(5).rxns).* abs(solution.x(model.protGroup(5).rxns));
        solution.phiEm = sum(costP_em); %phiEm
        subset = solution.x(1:length(solution.x)-2);     % Access elements from startIndex to endIndex
        costP_p = coefficients_Pc(1:length(solution.x)-2).* abs(subset);
        solution.phiP = sum(costP_p); %phiP used
        solution.phiCc = solution.phiPmax - solution.phiP ; %phiCc

        %Adding the sum of each sector in the membrane proteome
        costM_m = coefficients_Pm(1:length(solution.x)-2).* abs(subset);
        solution.phiM = sum(costM_m); % phiM 
        costM_c = coefficients_Pm(model.protGroup(1).rxns).* abs(solution.x(model.protGroup(1).rxns));
        solution.phiCm_m = sum(costM_c); %phiCm
        costM_em = coefficients_Pm(model.protGroup(5).rxns).* abs(solution.x(model.protGroup(5).rxns));
        solution.phiEm_m = sum(costM_em); %phiEm
        costM_er = coefficients_Pm(model.protGroup(4).rxns).* abs(solution.x(model.protGroup(4).rxns));
        solution.phiEr_m = sum(costM_er); %phiEr

    else
        solution.x = zeros( length(model.rxns), 1);
        solution.f = 0;
        
        solution.phiP_m = 0; solution.phiCm = 0; solution.phiCc = 0;solution.phiR = 0;...
        solution.phiEc = 0; solution.phiEr = 0; solution.phiEm = 0; solution.phiP = 0;
        solution.phiM =0; solution.phiCm_m =0; solution.phiEm_m =0; solution.phiEr_m =0;
        
    end
    
end
