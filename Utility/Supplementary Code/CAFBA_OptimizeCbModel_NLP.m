function [sol,convergenceData]=CAFBA_OptimizeCbModel_NLP(model,varargin)
% [sol,convergenceData]=CAFBA_OptimizeCbModel_NLP(model,varargin)
%
% This program iteratively optimizes the fluxes in the iJR904 model,
% with a growth rate dependent biomass function.
% The functions are described in the SI text of the paper. They can be
% modified in the 'updateBiomassComposition' here below.
%
% This function can also handle standard FBA optimization. In this case
% the 'protGroup' and 'w' fields are not necessary.
%
% INPUT
%   model : a valid iJR904 model with the CAFBA-specific fields
%           (see CAFBA_OptimizeCbModel.m)
%
% OUTPUT
%   sol  : solution to the NLP CAFBA problem
%          (see again CAFBA_OptimizeCbModel.m)
%
%   convergenceData : some informations about the convergence of the
%                     iterative procedure

% default settings
  % Maximum number of iterations
  maxIter=100;  
  
  % Solver. 'FBA', 'CAFBA_glpk' or 'CAFBA_gurobi'
  solver='CAFBA_glpk';
  
  % Accuracy of the iterative solver
  accuracy=1e-4;
  
  % Initial growth rate, in 1/h
  referenceBiomassLambda=0.6;  
  
  % verbosity. 0=silent
  verbose=0;  
% end of default settings

if nargin>1
    % each additional argument is a string identifying the variable, and
    % its value.
    
    % Make sure the number of optional arguments is even!
    nArgs=length(varargin);
    if round(nArgs/2)~=nArgs/2
        error('The number of optional arguments is not even!')
    end
    
    for i=1:2:(nargin-1)
        switch varargin{i}
            case 'maxIter'
                maxIter=varargin{i+1};
            case 'solver'
                solver=varargin{i+1};
            case 'accuracy'
                accuracy=varargin{i+1};
            case 'referenceBiomassLambda'
                referenceBiomassLambda=varargin{i+1};
            case 'verbose'
                verbose=varargin{i+1};
            otherwise
                error('Unrecognized optional parameter ''%s''',varargin{i});
        end
    end
end


% Initialize
exitFlag=0;

if strcmp(solver,'FBA')
    optimizingFun = @(model)optimizeCbModel(model);
elseif strcmp(solver,'CAFBA_glpk')
    optimizingFun = @(model)CAFBA_OptimizeCbModel_glpk(model);
elseif strcmp(solver,'CAFBA_gurobi')
    optimizingFun = @(model)CAFBA_OptimizeCbModel_gurobi(model);
end

% Initial condition
convergenceData.iter(1)=0;
convergenceData.lambda(1)=referenceBiomassLambda;
convergenceData.diff(1)=NaN;

% First step of the iteration procedure
model=updateBiomassComposition(model,referenceBiomassLambda);
sol=optimizingFun(model); 
lambdaOld=sol.f;
model=updateBiomassComposition(model,sol.f);

iter=1;
convergenceData.iter(iter+1)=1;
convergenceData.lambda(iter+1)=sol.f;
convergenceData.diff(iter+1)=abs(sol.f - referenceBiomassLambda);

while exitFlag==0    
    lambdaOld=sol.f;
    sol=optimizingFun(model);
    iter=iter+1; 
    convergenceData.iter(iter+1)=iter;
    convergenceData.lambda(iter+1)=sol.f;
    convergenceData.diff(iter+1)=abs(sol.f - lambdaOld);
    if abs(sol.f-lambdaOld)<accuracy || iter==maxIter
        exitFlag=1;
    else
        if verbose
            fprintf('%4i | Df = %f | newf = %f\n',iter,sol.f-lambdaOld,sol.f)
        end
        model=updateBiomassComposition(model,sol.f); 
    end
end

if iter==maxIter
    warning(['CAFBA_NLP_MaxIter',...
             'Warning: maximum number of iterations reached.\n',...
             'Last difference: Delta f = %f\t(accuracy=%f)\n'],...
             sol.f-lambdaOld,accuracy);
end
end

function model=updateBiomassComposition(model,lambda)
% here we set the growth-dependent biomass coefficients.
% See: J. Pramanik, J. D. Keasling, "Stoichiometric Model of Escherichia coli
% Metabolism: Incorporation of Growth-Rate Dependent Biomass Composition
% and MechanisticEnergy Requirements" (1997)

if ~isempty(strfind(model.description,'904'))
    r_bm=150;
    
    % Basic stoichiometric coefficients
    iJR904_betaAA=5.081; % sum of all st. coeffs.
    iJR904_betaDNA=0.1001; % sum of all st. coeffs.
    iJR904_betaRNA=0.636; % sum of all st. coeffs.
    iJR904_betaLIP=[0.000129;0.001935;0.000464;0.00052;0.154;0.0084;0.0276;0.035;0.007];
    
    % Average molecular masses
    mu_AA=0.110;
    mu_NT_DNA=0.303;
    mu_NT_RNA=0.320;
    mu_LIP=[69.708;35.656;37.155;37.805;0.162;3.877;0.9904;0.090;0.148]; % in g/mmol
    mu_LIP_avg=(iJR904_betaLIP' *mu_LIP)/sum(iJR904_betaLIP);
    
    % These are the functions specified in the SI
    % Definition of the Psi's
    psiDNA= 0.06 *(1+lambda.*lambda )./(1 + 3*lambda.*lambda);
    psiLIP= 0.07 + 0.2./(1+2*lambda);
    R=0.087 + lambda/4.5; % Scott et al. 2010
    psiOTHER=0.025;
    psiAA=(1-psiOTHER-psiDNA-psiLIP)./(1+R);
    psiRNA=(1-psiOTHER-psiDNA-psiLIP).*R./(1+R);
    
       
    % Inversion for the Beta's
    betaAA=psiAA/mu_AA;
    betaDNA=psiDNA/mu_NT_DNA;
    betaRNA=psiRNA/mu_NT_RNA;
    betaLIP=psiLIP/mu_LIP_avg;
    
    % ATP hydrolysis flux
    betaEnergy=45;

    % % iJR904 defaults
    % betaAA=5.081;
    % betaDNA=0.1001;
    % betaRNA=0.636;
    % betaGlyc=0.154;
    % betaEnergy=45.5608; % they are not exactly balanced... this is the
    %                     % coefficient for adp and h20
    
    %%%%%%%%%% Here we vary the individual stoichiometric coefficients
    
    % Energy: atp+h2o->adp+h+pi flux=45.5608*lambda mmol/g_{DW}h (see adp coefficient)
    r_atp=200; r_adp=151; r_h2o=424; r_h=427; r_pi=589;
    
    model.S(r_atp,r_bm)= - betaEnergy;
    model.S(r_h2o,r_bm)= - betaEnergy;
    model.S(r_adp,r_bm)= + betaEnergy;
    model.S(r_h  ,r_bm)= + betaEnergy;
    model.S(r_pi ,r_bm)= + betaEnergy;
    
    % DNA: datp, dctp, dgtp, dttp
    r_DNA=[254 ; 259 ; 267 ; 302];
    f_DNA=[0.0254; 0.0257; 0.0257; 0.0254];
    f_DNA=f_DNA/sum(f_DNA);    
    for i=1:4
        model.S(r_DNA(i),r_bm)=-f_DNA(i)*betaDNA;
    end
    
    % RNA: atp+ctp+gtp+utp flux = 0.636*lambda mmol/g_{DW}h
    % The slope in extracted from Pramanik paper
    r_RNA=[200; 240; 417; 749];
    f_RNA=[0.171; 0.126; 0.203; 0.136];
    f_RNA=f_RNA/sum(f_RNA);
    
    model.S(r_RNA(1),r_bm)= - f_RNA(1)*betaRNA - betaEnergy; % atp
    for i=2:4
        model.S(r_RNA(i),r_bm)= - f_RNA(i)*betaRNA;
    end
    
    % "Lipids"
    lip_r=[227;576;580;614;406;489;579;616;662];
    lip_f=iJR904_betaLIP./sum(iJR904_betaLIP);
    for i=1:length(lip_r)
        model.S(lip_r(i),r_bm)= - betaLIP * lip_f(i);
    end
    
    %lip_m=[69708;35656;37155;162;3877;990.4;90;148]; % molecular weight
    
    % Glycogen flux 0.154000*lambda mmol/g_{DW}h
    %r_glycogen=406;
    %model.S(r_glycogen,r_bm)= - betaGlyc;
    
    % List of AAs and biomass coefficients in iJR904
    % 169	  ala_L[c]	-0.488000
    % 192	  arg_L[c]	-0.281000
    % 195	  asn_L[c]	-0.229000
    % 197	  asp_L[c]	-0.229000
    % 244	  cys_L[c]	-0.087000
    % 382	  gln_L[c]	-0.250000
    % 388	  glu_L[c]	-0.250000
    % 393	    gly[c]	-0.582000
    % 436	  his_L[c]	-0.090000
    % 455	  ile_L[c]	-0.276000
    % 481	  leu_L[c]	-0.428000
    % 490	  lys_L[c]	-0.326000
    % 517	  met_L[c]	-0.146000
    % 583	  phe_L[c]	-0.176000
    % 611	  pro_L[c]	-0.210000
    % 651	  ser_L[c]	-0.205000
    % 692	  thr_L[c]	-0.241000
    % 707	  trp_L[c]	-0.054000
    % 714	  tyr_L[c]	-0.131000
    % 750	  val_L[c]	-0.402000
    aa_r=[169;192;195;197;244;...
        382;388;393;436;455;...
        481;490;517;583;611;...
        651;692;707;714;750];
    aa_f=[0.488000;0.281000;0.229000;0.229000;0.087000;...
        0.250000;0.250000;0.582000;0.090000;0.276000;...
        0.428000;0.326000;0.146000;0.176000;0.210000;...
        0.205000;0.241000;0.054000;0.131000;0.402000];  %stoich. coeff. from iJR904   
    aa_f=aa_f/sum(aa_f);
    for i=1:length(aa_r)
        model.S(aa_r(i),r_bm)= - betaAA * aa_f(i);
    end



else
    error(['This version of CAFBA_OptimizeCbModel_NLP only works with the',...
        'iJR904 E. coli model. Make sure the model is correct and there is a'...
        '''904'' string in the description field of the structure.'] );
end
end

% Useful things
% myAAs={'ala_L[c]','arg_L[c]','asn_L[c]','asp_L[c]','cys_L[c]',...
%        'gln_L[c]' ,'glu_L[c]' ,'gly[c]' ,'his_L[c]','ile_L[c]',...
%        'leu_L[c]','lys_L[c]' ,'met_L[c]','phe_L[c]','pro_L[c]',...
%        'ser_L[c]','thr_L[c]','trp_L[c]' ,'tyr_L[c]','val_L[c]'};
% myAAs={'ala-L[c]','arg-L[c]','asn-L[c]','asp-L[c]','cys-L[c]',...
%        'gln-L[c]' ,'glu-L[c]' ,'gly[c]' ,'his-L[c]','ile-L[c]',...
%        'leu-L[c]','lys-L[c]' ,'met-L[c]','phe-L[c]','pro-L[c]',...
%        'ser-L[c]','thr-L[c]','trp-L[c]' ,'tyr-L[c]','val-L[c]'};
