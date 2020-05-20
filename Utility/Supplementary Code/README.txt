%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This package is given as Supporting Material software to the article:

              CONSTRAINED ALLOCATION FLUX BALANCE ANALYSIS

Copyright (C) 2015  M. Mori, T. Hwa, O. C. Martin, A. De Martino & E. Marinari

Last edit: June, 3rd, 2015

Bug reports, comments and suggestions are welcome, and can be sent to
[ matteo.mori87_at_gmail.com ]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Contents:
1) DESCRIPTION
2) PREREQUISITES
3) CAFBA WORKFLOW

=== 1: DESCRIPTION ===

This package provides COBRA-compatible Matlab code in order to perform
Constraint Allocation Flux Balance Analysis (CAFBA) on constraint-based
metabolic models. It contains a set of COBRA-compatible Matlab functions
(.m files):

CAFBA_OptimizeCbModel_glpk.m
   Solves a CAFBA problem using the GLPK solver.

CAFBA_OptimizeCbModel_gurobi.m
   Solves a CAFBA problem using the GUROBI solver.
      
CAFBA_OptimizeCbModel_NLP.m
   Iteratively solves an FBA or CAFBA problem with growth-dependent biomass
   functions. It uses either COBRA's optimizeCbModel() function or one of
   the two CAFBA functions.
   
addProteinGroupsToModel.m
   This function adds a field 'protGroup' to the COBRA model with the
   informations about proteome groups (phi_C, phi_E, phi_R, phi_Q)
   
setWeights.m
   Adds or modifies the 'w' field to the COBRA model. This allows to set
   the weights for each different proteome group.
   
setWeights_rand.m
   Adds or modifies the 'w' field in the COBRA model. This allows to
   randomly generate the weights within a proteome group using some
   predefined probability distributions.
   
=== 2: PREREQUISITES ===

CAFBA optimizes models which are produced by the COBRA readCbModel(),
which converts SMBL models to a Matlab structure.

CAFBA is an LP problem, therefore requiring an appropriate solver.
Two versions of the main function CAFBA_optimizeCbModel are provided.
The first one (CAFBA_optimizeCbModel_glpk) uses the GLPK solver, with its
Matlab interface GLPKMEX. Both are downloadable at the following websites:

    https://www.gnu.org/software/glpk/
    http://glpkmex.sourceforge.net/

The second version (CAFBA_optimizeCbModel_gurobi) uses the non-free
Gurobi solver (it is available to academics with a free license, though).

CAFBA has been tested with various E. coli models. Other models can be used
if the function addProteinGroupsToModel() is modified accordingly, in
particular in the definitions of the E-sector.

=== 3: CAFBA WORKFLOW ===

We show in this section some examples of CAFBA workflow. A detailed "recipe" for
running CAFBA in general substrates can be found in Supplementary Note EV3 in
the paper.
The following codes can be directly copied and pasted into the Matlab command
window, or in Matlab scripts. Alternatively, we also provide the codes as
standalone scripts (example_setUpModel.m and example_carbonLimitation.m)

=== 3.1: SET UP THE MODEL FOR CAFBA (example_setUpModel.m) ===

% Read an E. coli SBML model (e.g. iJR904 model) using COBRA.
model = readCbModel();

% Remove bounds for glucose and oxygen
ex_glc_r=find(strcmp(model.rxns,'EX_glc(e)'));
ex_o2_r=find(strcmp(model.rxns,'EX_o2(e)'));
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
% (note that the function creates the field model.weights the first time it
% is called)
model=setWeights(model,1,0);
model=setWeights(model,2,0.00083);
model=setWeights(model,3,0.169);
model=setWeights(model,4,0);

% Clear useless variables
clear('ex_glc_r','ex_o2_r','glcd_r');

=== 3.2: YOUR FIRST CAFBA OPTIMIZATION ===

% After we have set up the model, we can start optimizing fluxes.
% In this case we use the glpk solver.

sol = CAFBA_OptimizeCbModel_glpk(model);

% With the gurobi solver it would be:
% sol = CAFBA_OptimizeCbModel_gurobi(model);

% As in usual FBA solutions obtained with COBRA function optimizeCbModel(),
% the optimal growth rate is sol.f (it should be 0.9964/h),
% while the optimal flux configuration is given by sol.x.
% Complete description of the output is contained into the source code of
% the CAFBA_optimizeCbModel() function .

=== 3.3: CARBON LIMITATION (example_carbonLimitation.m) ===

% Here we progressively increase the weight of the first proteome group
% (therefore reducing the glucose uptake rate). Read and set up the
% iJR904 model as done in section 3.1.

% This is a good wC mesh
Npoints=100;
wvec= (0.01*linspace(0,1,Npoints)  + 0.99*linspace(0,1,Npoints) .^3 );

% Computation will take some seconds
for i=1:Npoints
  model=setWeights(model,wvec(i),1);
  sol(i) = CAFBA_OptimizeCbModel_glpk(model);   
end

% We can now plot fluxes against growth rate. 
% In the iJR904 model reaction indexes are as follows:
%   Biomass --> 150 
%   Glucose exchange --> 344 
%   Acetate exchange --> 294 
%   AKGDH flux --> 78

[gr,glc,ac,akgdh]=deal(zeros(Npoints,1));
for i=1:Npoints
	gr(i)=sol(i).x(150);
	glc(i)=-sol(i).x(344); % note the minus sign
	ac(i) = sol(i).x(294);
	akgdh(i)=sol(i).x(78);
end 
figure(); hold on;
plot(gr,glc,'-+','Color',[0 .9 0],'LineWidth',2);
plot(gr,ac,'-+','Color',[.9 0 0],'LineWidth',2);
plot(gr,akgdh,'-+','Color',[0 0 .9],'LineWidth',2);
hold off; box on;
legend({'Glucose uptake','Acetate excretion','AKGDH (TCA) flux'},...
        'Location','NorthWest');
xlabel('Growth rate (1/h)');
ylabel('Flux (mmol/g_{DW}h)');




   
   

