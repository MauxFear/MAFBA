function solution = CAFBA_OptimizeCbModel_glpk(model,varargin)
% solution = CAFBA_OptimizeCbModel_glpk(model,varargin)
%
% This program computes a CAFBA optimization using the glpk solver.
% Make sure that the function glpk() actually works!
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
    constraintSense='S';
    
    % Verbosity
    verbose=1; % 1 = only GLPK error messages are displayed.
    
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
            case 'constraintSense'
                constraintSense=varargin{i+1};
            case 'removePhiFromFluxSplitting'
                removePhiFromFluxSplitting=varargin{i+1};
            case 'verbose'
                verbose=varargin{i+1};
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
% for every reaction in the CAFBA constraint. We have thus a number of fake mets in the
% interval [nMets+2*nConstr,nMets+4*nConstr], depending on the bounds. The extra
% constraint is:
% a (z_1 + z_2 + ... +z_n )/inKc + b x_growth = c

[nMets,nRxns]=size(model.S);
[solution.x,solution.f,solution.status,solution.extra]=deal(0);

% Checking input model
if ~isfield(model,'w')
    error('Invalid CAFBA model: ''w'' field missing.');
end
if ~isfield(model,'protGroup')
    error('Invalid CAFBA model: ''protGroup'' field missing.');
end

% SETTING VARS AND S MATRIX

% b: fill the first elements of the RHS vector if not provided
if (~isfield(model,'b'))
    bVec = zeros(nMets,1);
else
    bVec = model.b;
end

% c: first constraints are equalities %
ctypeVec(1:nMets,1)= 'S';

% Direction of the optimization
if strcmp(optSense,'max')
    LPproblem.sense=-1; % Maximization
elseif strcmp(optSense,'min')
    LPproblem.sense=1; % Minimization
else
    error('optSense not recognized.');
end

% Finding all non-zero w entries and adding the corresponding fluxes 
% to the three lists (depending on the bounds of the reactions)
positiveVec=[]; nPositive=0;
negativeVec=[]; nNegative=0;
mixedVec=[]; nMixed=0;
for r=1:nRxns
    if model.w(r)~=0
        if model.ub(r)>=0 && model.lb(r)>=0
            nPositive=nPositive+1;
            positiveVec(nPositive)=r;
        elseif model.ub(r)<=0 && model.lb(r)<0
            nNegative=nNegative+1;
            negativeVec(nNegative)=r;
        else
            nMixed=nMixed+1;
            mixedVec(nMixed)=r;
        end
    end;
end;

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
nRxnsExt=nRxns+nMixed;

% Splitting the fluxes
for n=1:nRxns
    LPproblem.vartype(n)='C';
    LPproblem.lb(n)=model.lb(n);
    LPproblem.ub(n)=model.ub(n);
    LPproblem.c(n)=model.c(n);
end
Smixed = sparse(nMets,nMixed); 
for n=1:nMixed
    LPproblem.vartype(nRxns+n)='C';
    LPproblem.lb(mixedVec(n))=0; % the original one is made irreversible
    LPproblem.lb(nRxns+n)=0; % a new, positive, variable is introduced
    LPproblem.ub(nRxns+n)=-model.lb(mixedVec(n));
    LPproblem.c(nRxns+n)=-model.c(mixedVec(n)); % obj. vec.
    
    Smixed(:,n)=-model.S(:,mixedVec(n)); % stoich. mat.
end

% SETTING UP THE ALLOCATION CONSTRAINT
if strcmp(optMode,'flux')
    
    % here we go with the additional constraint.
    
    % Coefficients of the fluxes.
    % Note: to insert the absolute value of the fluxes, we put some
    % minus signs somewhere here
    Sac=sparse(1,nRxnsExt);
    for n=1:nPositive
        Sac(1,positiveVec(n))=  model.w(positiveVec(n));
    end
    for n=1:nNegative
        Sac(1,negativeVec(n))= -model.w(negativeVec(n)); % minus sign!
    end
    for n=1:nMixed
        Sac(1,mixedVec(n))   =  model.w(mixedVec(n));
        Sac(1,nRxns+n)       =  model.w(mixedVec(n)); % here a + is needed, since both x+ and x- are >0
    end
    
    % right hand side
    mySum=0;
    for n=1:length(model.protGroup)
        mySum=mySum+model.protGroup(n).phi0;
    end
    bVec(nMets+1)=1-mySum;
    
    % type of constraint ( 'S', 'U', 'L', as in GLPK )
    ctypeVec(nMets+1)=constraintSense;
else
    Sac=sparse(0,nRxnsExt);
end

% OTHER SETTINGS AND CALL TO SOLVER
% building the extended S matrix, as well as the constraints parameters
LPproblem.A=[[model.S , Smixed] ; Sac];

LPproblem.b=bVec; % constraints rhs
LPproblem.ctype=ctypeVec; % constraint type

% if we are not maximizing the objective function but some phi group
% we have to change the objective function

if strcmp(optMode,'phi')
    % Maximizes/minimizes the proteome fraction of the sector labeled by
    % optProtGroup variable
    LPproblem.c=zeros(nRxnsExt,1);
    for n=1:length(positiveVec)
        r=positiveVec(n);
        if ~isempty(find(model.protGroup(optProtGroup).rxns==r, 1))
            LPproblem.c(r)  =  model.w(r);
        end
    end
    for n=1:length(negativeVec)
        r=negativeVec(n);
        if ~isempty(find(model.protGroup(optProtGroup).rxns==r, 1))
            LPproblem.c(r)  = -model.w(r);
        end
    end
    for n=1:length(mixedVec)
        rp=mixedVec(n);
        if ~isempty(find(model.protGroup(optProtGroup).rxns==rp, 1))
            rm=nRxns+n;
            LPproblem.c(rp) =  model.w(rp);
            LPproblem.c(rm) =  model.w(rp);
        end
    end
elseif ~strcmp(optMode,'flux')
    error('optMode variable not recognized. Exiting.');
end

% Other settings
LPproblem.param.msglev=verbose;

% Solving the problem with a call to GLPK
[solution.x,solution.f,solution.status,solution.extra] = ...
    glpk(LPproblem.c, LPproblem.A, LPproblem.b, LPproblem.lb, ...
    LPproblem.ub, LPproblem.ctype, LPproblem.vartype, ...
    LPproblem.sense, LPproblem.param );

% Remove contributions to phi due to flux splitting
if strcmp(optMode,'phi') && removePhiFromFluxSplitting
    for i=1:nMixed
        if ~isempty(find(model.protGroup(optProtGroup).rxns==mixedVec(i), 1))            
            solution.f=solution.f- ...
                2*min([solution.x(mixedVec(i)),solution.x(nRxns+i)])...
                *model.w(mixedVec(i));            
        end
    end
end

% Put together fluxes after splitting
for i=1:nMixed
    solution.x(mixedVec(i))=solution.x(mixedVec(i))-solution.x(nRxns+i);
end
solution.x=solution.x(1:nRxns,1);

end