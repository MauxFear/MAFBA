% Read an E. coli SBML model (e.g. iJR904 model) using COBRA.
model_iJR_M3 = readCbModel('D:\Code\MATLAB\Projects\CAFBA\Models\Ec_iJR904_flux1.xml');

% Remove default bounds for glucose and oxygen
ex_glc_r=find(strcmp(model_iJR_M3.rxns,'EX_glc_e_'));
ex_o2_r=find(strcmp(model_iJR_M3.rxns,'EX_o2_e_'));
ex_co2_r=find(strcmp(model_iJR_M3.rxns,'EX_co2_e_'));

model_iJR_M3.lb(ex_o2_r)=0; model_iJR_M3.ub(ex_o2_r)=1000;
model_iJR_M3.lb(ex_co2_r)=0; model_iJR_M3.ub(ex_co2_r)=1000;


% Knocking out the GLCDe flux
glcd_r=find(strcmp(model_iJR_M3.rxns,'GLCDe'));
model_iJR_M3.lb(glcd_r)=0;
model_iJR_M3.ub(glcd_r)=0;

model_iJR_M3 = addReaction(model_iJR_M3, 'EX_pre_z','reactionFormula', 'pre_z <=> ');
model_iJR_M3 = addReaction(model_iJR_M3, 'RxnZ','reactionFormula', 'pre_z -> prot_Z', 'subsystem', 'Unassigned', 'upperBound', 1);
model_iJR_M3 = addReaction(model_iJR_M3, 'EX_Z','reactionFormula', 'prot_Z <=> ');

model_iJR_M3 = addReaction(model_iJR_M3, 'EX_pre_y','reactionFormula', 'pre_y <=> ');
model_iJR_M3 = addReaction(model_iJR_M3, 'RxnY','reactionFormula', 'pre_y -> prot_Y', 'subsystem', 'Putative Transporters', 'upperBound', 1);
model_iJR_M3 = addReaction(model_iJR_M3, 'EX_Y','reactionFormula', 'prot_Y <=> ');

%% 


% Add protein groups to the model. 
model_iJR_M3.lb(ex_glc_r)=-1000; model_iJR_M3.ub(ex_glc_r)=1000;

model_iJR_M3=addProteinGroups_MAFBA_v3(model_iJR_M3);
% % 
% model_iJR_M3.lb(ex_glc_r)=0; model_iJR_M3.ub(ex_glc_r)=1000;
% 
% model_iJR_M3=addProteinGroups_MAFBA_v3(model_iJR_M3, 'g6p');


% For example, this is the code to model lactose minimal medium:
%   model=addProteinGroupsToModel(model,'C-lim lac');
% Look inside the function's source code for more informations.
%% 


% We now set the offsets for the different groups.
model_iJR_M3.protGroup(1).phi0 = 0; % phiCm_0
model_iJR_M3.protGroup(2).phi0 = 0; % phiEc_0
model_iJR_M3.protGroup(3).phi0 = 0.066; % phiR_0
model_iJR_M3.protGroup(4).phi0 = 0; % phiEr_0
model_iJR_M3.protGroup(5).phi0 = 0; % phiEm_0
model_iJR_M3.protGroup(6).phi0 = 0.45; % phiQ_0
model_iJR_M3.protGroup(7).phi0 = 0; % phiCc_0

weight = 0.0007596; % 0.0007596
% Now we have to set the weights for the different groups.
% (note that the function creates the field model.weights if needed)
model_iJR_M3=setWeights(model_iJR_M3,1,weight);
model_iJR_M3=setWeights(model_iJR_M3,2,weight ); % 0.00083
model_iJR_M3=setWeights(model_iJR_M3,4,weight );
model_iJR_M3=setWeights(model_iJR_M3,5,weight ); % 0.00062
model_iJR_M3=setWeights(model_iJR_M3,3,0.169);
model_iJR_M3=setWeights(model_iJR_M3,6,0);
model_iJR_M3=setWeights(model_iJR_M3,7,0);


%Run the MAFBA
sol_iJR_m3 = MAFBA_OptCbModel_cplex_v4(model_iJR_M3, 'gamma', 1, 'cSense_m', 'U');
fprintf('The maximum growth rate is %1.4f \n', sol_iJR_m3.f);

% Clear useless variables
clear('ex_glc_r','ex_o2_r', 'ex_co2_r', 'glcd_r', 'weight');