% Here we progressively increase the weight of the first proteome group
% (therefore reducing the glucose uptake rate). Read and set up the
% iJR904 model as done in section 3.1.

% We now set the offsets for the different groups.
model_iJR.protGroup(1).phi0 = 0;
model_iJR.protGroup(2).phi0 = 0;
model_iJR.protGroup(3).phi0 = 0.066;
model_iJR.protGroup(4).phi0 = 0.45;

% Now we have to set the weights for the different groups.
% (note that the function creates the field model.weights if needed)
model_iJR=setWeights(model_iJR,1,0.00083);
model_iJR=setWeights(model_iJR,2,0.00083); % 0.00083
model_iJR=setWeights(model_iJR,3,0.169);
model_iJR=setWeights(model_iJR,4,0);



%% 

% This is a good wC mesh
Npoints=100;
wvec= (0.01*linspace(0,1,Npoints)  + 0.99*linspace(0,1,Npoints) .^3 );

% We can now plot fluxes against growth rate. 
% In the iJR904 model reaction indexes are as follows:
%   Biomass --> 150 
%   Glucose exchange --> 344 
%   Acetate exchange --> 294 
%   AKGDH flux --> 78
glc_r=find(strcmp(model_iJR.rxns,'EX_glc_e_'));
bm_r=model_iJR.protGroup(3).rxns;
akgdh_r=find(strcmp(model_iJR.rxns,'AKGDH'));
ace_r=find(strcmp(model_iJR.rxns,'EX_ac_e_'));
mals_r=find(strcmp(model_iJR.rxns,'MALS'));
edd_r=find(strcmp(model_iJR.rxns,'EDD'));
memb_r=find(strcmp(model_iJR.rxns,'ATPS4r'));
co2_r=find(strcmp(model_iJR.rxns,'EX_co2_e_'));
% Glycolitic Carbon sources
EX_g6p=find(strcmp(model_iJR.rxns,'EX_g6p_e_'));
EX_man=find(strcmp(model_iJR.rxns,'EX_man_e_')); % mannose
EX_lcts=find(strcmp(model_iJR.rxns,'EX_lcts_e_') ); % lactose transport
EX_glyc=find(strcmp(model_iJR.rxns,'EX_glyc_e_')); % glycerol
EX_gal=find(strcmp(model_iJR.rxns,'EX_gal_e_'));
% Non-glycolitic carbon sources
EX_fum=find(strcmp(model_iJR.rxns,'EX_fum_e_'));
EX_succ=find(strcmp(model_iJR.rxns,'EX_succ_e_')); % succinate
EX_pyr=find(strcmp(model_iJR.rxns,'EX_pyr_e_')); % pyruvate
EX_lac=find(strcmp(model_iJR.rxns,'EX_lac_D_e_')); % lactate transport
EX_glcn=find(strcmp(model_iJR.rxns,'EX_glcn_e_'));



cs_v = {'glc','g6p','lcts','man','glyc'};
cs_ng = {'glc','fum','succ','pyr','lac','glcn'};
color = {sscanf('65AF65','%2x%2x%2x',[1 3])/255, sscanf('C74848','%2x%2x%2x',[1 3])/255, sscanf('9065D1','%2x%2x%2x',[1 3])/255,...
    sscanf('266DD7','%2x%2x%2x',[1 3])/255, sscanf('F2A1E4','%2x%2x%2x',[1 3])/255, sscanf('D2C725','%2x%2x%2x',[1 3])/255,...
    sscanf('B2B2B2','%2x%2x%2x',[1 3])/255,sscanf('EC8327','%2x%2x%2x',[1 3])/255};

    
    fig=figure(); 
    
    set(fig, 'Position',[10 10 800 600]);
    title({'CAFBA: Flux profile vs growth rate','varying w_{C} with different carbon sources'})
    
    for cs = 1:length(cs_ng)
        if cs == 1
            hold on
        end
        model_iJR.lb(EX_g6p)=0; model_iJR.lb(EX_gal)=0;
        model_iJR.lb(EX_lcts)=0; model_iJR.lb(EX_glyc)=0;
        model_iJR.lb(EX_fum)=0; model_iJR.lb(EX_succ)=0;
        model_iJR.lb(EX_pyr)=0; model_iJR.lb(EX_lac)=0;
        model_iJR.lb(EX_glcn)=0; model_iJR.lb(glc_r)=-1000;
        
        model_iJR=addProteinGroupsToModel(model_iJR, cs_ng(cs));

        % Computation will take some seconds
        [gr,glc,ac,akgdh,mals,edd,memb,co2]=deal(zeros(Npoints,1));
        for i=1:Npoints
            model_iJR=setWeights(model_iJR,1,wvec(i));
            sol_M(i) = CAFBA_OptCbModel_cplex(model_iJR); %0.2123966  
            gr(i)=sol_M(i).x(bm_r);
            glc(i)=-sol_M(i).x(glc_r); % note the minus sign
            ac(i) = sol_M(i).x(ace_r);
            akgdh(i)=sol_M(i).x(akgdh_r);
            mals(i)=sol_M(i).x(mals_r);
            edd(i)=sol_M(i).x(edd_r);
            memb(i)=0.25* sol_M(i).x(memb_r);
            co2(i)=sol_M(i).x(co2_r);
        end 
        T = table(wvec', gr, ac, glc, memb, akgdh, mals, edd , co2);
        T = sortrows(T,'gr');
        
        %Add detect acetate threshold >=0.1
        dif_v = diff(T.memb);
        memb_s = find(dif_v <= 0.1);
        t_i = find(T.ac(memb_s) >= 0.01, 1, 'first');
        ac_t = memb_s(t_i);
        
        if isempty(ac_t)
            memb_s = find(dif_v <= 0.2);
            t_i = find(T.ac(memb_s) >= 0, 1, 'first');
            ac_t = memb_s(t_i);
        end
        if isempty(ac_t)
            ac_t = 1;
        end
        
        gr_ac=T.gr(ac_t);
        ctr = char(cs_ng(cs));
        str_cs = ['cs: ', ctr]; 
        plot(T.gr,T.ac,'-^','Color',color{cs},'LineWidth',1, 'markerfacecolor', color{cs}, 'MarkerSize', 5, 'DisplayName',['acetate flux (', str_cs, ')' ] );
        
%         plot(T.gr,T.memb,'-s','Color',color{cs},'LineWidth',1, 'markerfacecolor', color{cs}, 'MarkerSize', 5 , 'DisplayName',['ATPS4r flux 0.25x (', str_cs, ')' ]);
        max_ac(cs) = max(ac);
        plot([gr_ac,gr_ac], [0, 100],'--', 'Color',color{cs}, 'LineWidth',1, 'DisplayName',['\lambda_{ac} (', str_cs, ')' ]);
    
    end
    set(gca,'FontSize',12);  
    xlabel('Growth rate (h^{-1})', 'Fontsize', 14);
    ylabel('Flux (mmol/g_{DW}h)', 'Fontsize', 14);
    lgd = legend;
    set(lgd, 'Box', 'off', 'Location', 'west','Fontsize', 9,...
        'Position',[0.125214521452145 0.519773720831774 0.303217827211512 0.270695370396242]);
    max_y = max(max_ac) + 0.5;
    axis([0 inf 0 max_y]);
    hold off;
%     filename = ['MAFBA_C_ace_g_', str,'.png'];
%     saveas(fig,filename,'png');



