% Read an E. coli SBML model (e.g. iJR904 model) using COBRA.
model = readCbModel('Ec_iJR904_flux1.xml');

% Remove default bounds for glucose and oxygen
ex_glc_r=find(strcmp(model.rxns,'EX_glc_e_'));
ex_o2_r=find(strcmp(model.rxns,'EX_o2_e_'));
model.lb(ex_glc_r)=-1000; model.ub(ex_glc_r)=1000;
model.lb(ex_o2_r)=-1000; model.ub(ex_o2_r)=1000;

% Knocking out the GLCDe flux
glcd_r=find(strcmp(model.rxns,'GLCDe'));
model.lb(glcd_r)=0;
model.ub(glcd_r)=0;


% Add protein groups to the model. If called without the second argument,
% this function includes the glucose exchange reaction 'EX_glc' into the
% first proteome group. By default, the function also removes the standard
% lower bound to the glucose exchange reaction ( -10 mmol/gDWh, which is
% brought to -1000). The groups are numbered in this way:
% 1 : C
% 2 : E
% 3 : R
% 4 : Q
model=addProteinGroupsToModel(model);

% For example, this is the code to model lactose minimal medium:
%   model=addProteinGroupsToModel(model,'C-lim lac');
% Look inside the function's source code for more informations.

% We now set the offsets for the different groups.
model.protGroup(1).phi0 = 0;
model.protGroup(2).phi0 = 0;
model.protGroup(3).phi0 = 0.066;
model.protGroup(4).phi0 = 0.45;

% Now we have to set the weights for the different groups.
% (note that the function creates the field model.weights if needed)
model=setWeights(model,1,0);
model=setWeights(model,2,0.00083);
model=setWeights(model,3,0.169);
model=setWeights(model,4,0);

% Clear useless variables
clear('ex_glc_r','ex_o2_r');