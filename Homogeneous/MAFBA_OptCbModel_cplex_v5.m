function solution = MAFBA_OptCbModel_cplex_v4(model,varargin)
% solution = MAFBA_OptimizeCbModel_cplex(model,varargin)
%
% This program computes a CAFBA optimization using the cplex solver.
% CPLEX 12.8 version
%
% The only necessary input is the model. Optional arguments can be
% specified to modify default values.
%
% INPUT:
%
% model: this must be a model returned by the Cobra Toolbox function
%        'readCbModel()' with additional 'w' and 'protGroup' fields.
%
% Optional arguments can be used to vary default settings, as:
%          varargin{i}=varargin{i+1}  with i=1,3,5,...
% The odd args are the variabile names, the even ones are the values.
%
% OUTPUT:
%
% solution: structure with the following fields
%  x:         vector of optimized fluxes
%  f:         value of the objective function (e.g. biomass production rate)
%  status:    status of the optimization. Its value depends on the LP solver.
%  extras:    other outputs from the LP solver (e.g. shadow prices)
%
% The output is very similar to the output of the glpk routine.
% See in the glpk help page for further information.
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% Defaults

    % Direction of the optimization. It can be either 'min' or 'max'
    optSense='max';

    % optMode is a string which can attain the following values:
    % 'flux'
    %     the standard objective function (the scalar product
    %     between the vector model.c and the fluxes) is maximized or
    %     minimized, depending on the optimization sense.
    % 'phi'
    %     A protein group is minimized or maximized, depending on the
    %     optimization sense. CAFBA constraint is not used here.
    %     The variable optProtGroup is the index of the proteome fraction
    %     to be optimized.
    optMode='flux';
    optProtGroup=2;
    
    % Only for the 'phi' mode: if the following flag is 1, contributions to
    % phi arising from loops between forward and backward fluxes are removed 
    removePhiFromFluxSplitting=1;
    
    % Allocation constraint: 'S' for equality, 'U' and 'L' for inequality
    % (<= or >=, respectively). These are the same labels used by the GLPK
    % solver.
    cSense_p='S';
    cSense_m = 'U';
    % Verbosity
    verbose=1; % 1 = only GLPK error messages are displayed.
    
    %PhiM_max
    phiM_max = 0.20;
    
    %Gamma ratio
    gamma = 0;
    
% End of defaults

% Optional arguments
optargin=size(varargin,2);
if floor(optargin/2)~=optargin/2
    error('Error, optargin is not an even number.');
elseif optargin>0
    for i=1:2:(optargin-1)
        switch varargin{i}
            case 'optSense'
                optSense=varargin{i+1};
            case 'optMode'
                optMode=varargin{i+1};
            case 'optProtGroup'
                optProtGroup=varargin{i+1};
            case 'cSense_p'
                cSense_p=varargin{i+1};
            case 'cSense_m'
                cSense_m=varargin{i+1};
            case 'removePhiFromFluxSplitting'
                removePhiFromFluxSplitting=varargin{i+1};
            case 'verbose'
                verbose=varargin{i+1};
            case 'phiMmax'
                phiM_max = varargin{i+1};
            case 'gamma'
                gamma = varargin{i+1};
            otherwise
                error('Optional argument not recognized.');
        end
    end
end

% TECHNICAL STUFF
%
% ADDITIONAL RXNS
% if x1,x2,...,xnRxns are the original fluxes, now we have
%    x1,x2,...,z1,z2,...,d1,d2,...dnRxns
% so we have nRxnsExt=3*nRxns.
%
% ADDITIONAL METS:
% if lb<0 && ub>0
%    put xp+xm in the allocation constraint
%    split the fluxes in the stoichiometric matrix
% if lb >= 0
%   put x in the allocation constraint
%
% if ub<=0
%   put -x in the allocation constraint
%
% for every reaction in the CAFBA constraint. We have thus a number of fake
% mets in the interval [nMets+2*nConstr,nMets+4*nConstr], depending on the
% bounds. The extra constraint is: a (z_1 + z_2 + ... +z_n )/inKc + b
% x_growth = c

[nMets,nRxns]=size(model.S);
[solution.x,solution.f,solution.status]=deal(0);

% Checking input model
if ~isfield(model,'w')
    error('Invalid MAFBA model: ''w'' field missing.');
end
if ~isfield(model,'protGroup')
    error('Invalid MAFBA model: ''protGroup'' field missing.');
end

% Start the cplex Object
cplex = Cplex('MAFBA');
cplex.DisplayFunc = [];


% SETTING VARS AND S MATRIX

% b: fill the first elements of the RHS vector if not provided
if (~isfield(model,'b'))
    bVec = zeros(nMets,1);
else
    bVec = model.b;
end

% c: first constraints are equalities %
lhsVec= bVec;

% Direction of the optimization
if strcmp(optSense,'max')
    cplex.Model.sense='maximize'; % Maximization
elseif strcmp(optSense,'min')
    cplex.Model.sense='minimize'; % Minimization
else
    error('optSense not recognized.');
end

% Finding all non-zero w entries and adding the corresponding fluxes 
% to the three lists (depending on the bounds of the reactions)
posVec= zeros(nRxns); nPos=0;
negVec= zeros(nRxns); nNeg=0;
mixVec= zeros(nRxns); nMix=0;
Mrxns = [model.protGroup(4).rxns;model.protGroup(1).rxns;model.protGroup(5).rxns]; % Membrane Reactions
posVec_m= zeros(nRxns); nPos_m=0;
negVec_m= zeros(nRxns); nNeg_m=0;
mixVec_m= zeros(nRxns); nMix_m=0;

for r=1:nRxns
    if model.w(r)~=0
        if model.ub(r)>=0 && model.lb(r)>=0 % irreversible positive reactions
            nPos=nPos+1;
            posVec(nPos)=r;
            if ismember(r, Mrxns, 'rows')
                nPos_m = nPos_m+1;
                posVec_m(nPos_m) = r;
            end
        elseif model.ub(r)<=0 && model.lb(r)<0 % irreversible negative reactions
            nNeg=nNeg+1;
            negVec(nNeg)=r;
            if ismember(r, Mrxns, 'rows')
                nNeg_m = nNeg_m+1;
                negVec_m(nNeg_m) = r;
            end
        else                                % reversible reactions
            nMix=nMix+1;
            mixVec(nMix)=r;
            if ismember(r, Mrxns, 'rows')
                nMix_m = nMix_m+1;
                mixVec_m(nMix_m) = r;
            end
        end
    end
end
posVec = posVec(1:nPos);
negVec = negVec(1:nNeg);
mixVec = mixVec(1:nMix);

posVec_m = posVec_m(1:nPos_m);
negVec_m = negVec_m(1:nNeg_m);
mixVec_m = mixVec_m(1:nMix_m);


Smixed = sparse(-model.S(:,mixVec)); % stoich. mat. for extra Reactions

% Adding nConstr variables, corresponding to the splitted vars
%
% Example. If one starts with (x1,x2,x3,x4), but only the second and
% the fourth have a positive weight w and the fourth is irreversibile, then
% nPositive=1, nNegative=0, nMixed=3 and the new set of variables is:
%   (x1,xp2,x3,x4,xm2)
% where the second flux has been splitted. In this case:
%   sum_X phi_X = w2*xp2 + w2*xm2 + w4*x4
% where X runs over all protein groups.

switch optMode
    case 'flux'
        nMetsExt=nMets+1;
    case 'phi'
        nMetsExt=nMets;
end
nRxnsExt=nRxns+nMix; % extend the number of reactions

% Splitting the fluxes
vartype = cell(nRxns,1);
vartype(:)={'C'};
lb=model.lb;
ub=model.ub;
c=model.c;

for n=1:nMix
    vartype(nRxns+n)={'C'};
    lb(mixVec(n))=0; % the original one is made irreversible, set value of v(+)
    lb(nRxns+n)=0; % a new, positive, variable is introduced
    ub(nRxns+n)=-model.lb(mixVec(n)); % set value of v(-)
    c(nRxns+n)=-model.c(mixVec(n)); % obj. vec.
end

% SETTING UP THE ALLOCATION CONSTRAINT
if strcmp(optMode,'flux')
    
    % here we go with the additional constraint.
    
    % Coefficients of the fluxes.
    % Note: to insert the absolute value of the fluxes, we put some
    % minus signs somewhere here
    rrxn = (nRxns+1):nRxnsExt;
    
    wv = model.w(posVec);
    wv = [wv; (-model.w(negVec))]; % minus sign!
    wv = [wv; (model.w(mixVec))];
    wv = [wv; model.w(mixVec)]; % here a + is needed, since both x+ and x- are >0
    
    Sac = sparse(1,[posVec,negVec,mixVec, rrxn],wv);
    
    wv_m = [model.w(posVec_m); (-model.w(negVec_m));... 
        (model.w(mixVec_m)); (model.w(mixVec_m))];
    mxr = find(ismember(mixVec,mixVec_m));
    mxr = nRxns + mxr;
    
    Sac_m = sparse(1,[posVec_m,negVec_m,mixVec_m, mxr],wv_m,1,nRxnsExt);
    
    
    
    % right hand side
    mySum=0;
    for n=1:length(model.protGroup)
        mySum=mySum+model.protGroup(n).phi0;
    end
    %adjusting phiM_max according to the gamma ratio : phiM_max / phi_max
    if gamma > 0
        phiM_max = gamma* (1-mySum);
    elseif gamma < 0
        error('gamma value cannot be negative.');
    end
    bVec(nMets+1)=1-mySum; % add another constraint value which is phi_max
    bVec(nMets+2)=phiM_max; % add another constraint value which is phiM_max
    
    % type of constraint ( 'S', 'U', 'L', as in GLPK ) https://octave.sourceforge.io/octave/function/glpk.html
    if cSense_p == 'S'
        lhsVec(nMets+1)=bVec(nMets+1);
    elseif cSense_p == 'U'
        lhsVec(nMets+1)=0;
    else
        error('Proteome constraint type not recognized, use S for equalities or U for inequalities');
    end
    if cSense_m == 'S'
        lhsVec(nMets+2)= bVec(nMets+2);
    elseif cSense_m == 'U'
        lhsVec(nMets+2)=0;
    else
        error('Membrane constraint type not recognized, use S for equalities or U for inequalities');
    end
else
    Sac=sparse(0,nRxnsExt);
    Sac_m = sparse(0,nRxnsExt);
end

% OTHER SETTINGS AND CALL TO SOLVER
% building the extended S matrix, as well as the constraints parameters
A=[[model.S , Smixed] ; Sac;Sac_m];

% if we are not maximizing the objective function but some phi group
% we have to change the objective function

if strcmp(optMode,'phi')
    % Maximizes/minimizes the proteome fraction of the sector labeled by
    % optProtGroup variable
    c=zeros(nRxnsExt,1);
    for n=1:length(posVec)
        r=posVec(n);
        if ~isempty(find(model.protGroup(optProtGroup).rxns==r, 1))
            c(r)  =  model.w(r);
        end
    end
    for n=1:length(negVec)
        r=negVec(n);
        if ~isempty(find(model.protGroup(optProtGroup).rxns==r, 1))
            c(r)  = -model.w(r);
        end
    end
    for n=1:length(mixVec)
        rp=mixVec(n);
        if ~isempty(find(model.protGroup(optProtGroup).rxns==rp, 1))
            rm=nRxns+n;
            c(rp) =  model.w(rp);
            c(rm) =  model.w(rp);
        end
    end
elseif ~strcmp(optMode,'flux')
    error('optMode variable not recognized. Exiting.');
end

% Other settings
% LPproblem.param.msglev=verbose;

%Add the data to the CPLEX model
cplex.addCols(c, [], lb, ub);
cplex.addRows(lhsVec(1:nMets), A(1:nMets,:), bVec(1:nMets));
cplex.Param.lpmethod.Cur = 3;
cplex.solve();

% Add the complicated constraints
cplex.addRows(lhsVec(nMets+2), A(nMets+2,:), bVec(nMets+2));

% Re-solve the model with the dual simplex method
cplex.Param.lpmethod.Cur = 2;
cplex.solve();

cplex.addRows(lhsVec(nMets+1), A(nMets+1,:), bVec(nMets+1));
cplex.Param.lpmethod.Cur = 2;
cplex.solve();

%% 

%Construct the solution structure
solution.f = cplex.Solution.objval;
solution.x = cplex.Solution.x;
solution.status = cplex.Solution.status;
solution.statusstring = cplex.Solution.statusstring;
solution.method = cplex.Solution.method;
solution.extras.time = cplex.Solution.time;
solution.extras.dettime = cplex.Solution.dettime;
solution.extras.dual = cplex.Solution.dual;

solution.extras.basis = cplex.Solution.basis;
solution.extras.reducedcost = cplex.Solution.reducedcost;
solution.extras.ax = cplex.Solution.ax;

% Remove contributions to phi due to flux splitting
if strcmp(optMode,'phi') && removePhiFromFluxSplitting
    for i=1:nMix
        if ~isempty(find(model.protGroup(optProtGroup).rxns==mixVec(i), 1))            
            solution.f=solution.f- ...
                2*min([solution.x(mixVec(i)),solution.x(nRxns+i)])...
                *model.w(mixVec(i));            
        end
    end
end

% Put together fluxes after splitting
for i=1:nMix
    solution.x(mixVec(i))=solution.x(mixVec(i))-solution.x(nRxns+i);
end

solution.x=solution.x(1:nRxns,1); 
%Adding the sum of each sector in the solution
wx_m = model.w(Mrxns).* abs(solution.x(Mrxns));
solution.phiM = sum(wx_m); % phiM
wx_c = model.w(model.protGroup(1).rxns).* abs(solution.x(model.protGroup(1).rxns));
solution.phiCm = sum(wx_c); %phiCm
wx_cc = model.w(model.protGroup(7).rxns).* abs(solution.x(model.protGroup(7).rxns));
solution.phiCc = sum(wx_cc); %phiCc
wx_r = model.w(model.protGroup(3).rxns).* abs(solution.x(model.protGroup(3).rxns));
solution.phiR = sum(wx_r); %phiR
wx_ec = model.w(model.protGroup(2).rxns).* abs(solution.x(model.protGroup(2).rxns));
solution.phiEc = sum(wx_ec); %phiEc
wx_er = model.w(model.protGroup(4).rxns).* abs(solution.x(model.protGroup(4).rxns));
solution.phiEr = sum(wx_er); %phiEr
wx_em = model.w(model.protGroup(5).rxns).* abs(solution.x(model.protGroup(5).rxns));
solution.phiEm = sum(wx_em); %phiEm
wx_p = model.w.* abs(solution.x);
solution.phi_max = sum(wx_p); %phi max
end