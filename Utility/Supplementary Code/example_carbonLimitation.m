% Here we progressively increase the weight of the first proteome group
% (therefore reducing the glucose uptake rate). Read and set up the
% iJR904 model as done in section 3.1.

% This is a good wC mesh
Npoints=100;
wvec= (0.01*linspace(0,1,Npoints)  + 0.99*linspace(0,1,Npoints) .^3 );

% Computation will take some seconds
for i=1:Npoints
  model=setWeights(model,1,wvec(i));
  sol(i) = CAFBA_OptimizeCbModel_glpk(model);   
end

% We can now plot fluxes against growth rate. 
% In the iJR904 model reaction indexes are as follows:
%   Biomass --> 150 
%   Glucose exchange --> 344 
%   Acetate exchange --> 294 
%   AKGDH flux --> 78
glc_r=model.protGroup(1).rxns;
bm_r=model.protGroup(3).rxns;
akgdh_r=find(strcmp(model.rxns,'AKGDH'));
ace_r=find(strcmp(model.rxns,'EX_ac_e_'));

[gr,glc,ac,akgdh]=deal(zeros(Npoints,1));
for i=1:Npoints
	gr(i)=sol(i).x(bm_r);
	glc(i)=-sol(i).x(glc_r); % note the minus sign
	ac(i) = sol(i).x(ace_r);
	akgdh(i)=sol(i).x(akgdh_r);
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
