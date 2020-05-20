% Here we progressively increase the weight of the first proteome group
% (therefore reducing the glucose uptake rate). Read and set up the
% iJR904 model as done in section 3.1.

% We now set the offsets for the different groups.
model_iJR_M3.protGroup(1).phi0 = 0; % phiC_0
model_iJR_M3.protGroup(2).phi0 = 0; % phiEc_0
model_iJR_M3.protGroup(3).phi0 = 0.066; % phiR_0
model_iJR_M3.protGroup(4).phi0 = 0; % phiEr_0
model_iJR_M3.protGroup(5).phi0 = 0; % phiEm_0
model_iJR_M3.protGroup(6).phi0 = 0.45; % phiQ_0
model_iJR_M3.protGroup(7).phi0 = 0; % phiCc_0

weight = 0.0007596;

% Now we have to set the weights for the different groups.
% (note that the function creates the field model.weights if needed)
model_iJR_M3=setWeights(model_iJR_M3,1,weight);
model_iJR_M3=setWeights(model_iJR_M3,2,weight ); % 0.00083
model_iJR_M3=setWeights(model_iJR_M3,4,weight );
model_iJR_M3=setWeights(model_iJR_M3,5,weight ); % 0.00062
model_iJR_M3=setWeights(model_iJR_M3,3,0.169);
model_iJR_M3=setWeights(model_iJR_M3,6,0);
model_iJR_M3=setWeights(model_iJR_M3,7,0);
%% 

Npoints=100;
wvec= (0.4+ 5* linspace(0,1,Npoints)  + 20.5*linspace(0,1,Npoints) .^3 );
ex_glc_r=find(strcmp(model_iJR_M3.rxns,'EX_glc_e_'));

% We can now plot fluxes against growth rate. 
% In the iJR904 model reaction indexes are as follows:
%   Biomass --> 150 
%   Glucose exchange --> 344 
%   Acetate exchange --> 294 
%   AKGDH flux --> 78
glc_r=find(strcmp(model_iJR_M3.rxns,'EX_glc_e_'));
bm_r=model_iJR_M3.protGroup(3).rxns;
akgdh_r=find(strcmp(model_iJR_M3.rxns,'AKGDH'));
ace_r=find(strcmp(model_iJR_M3.rxns,'EX_ac_e_'));
mals_r=find(strcmp(model_iJR_M3.rxns,'MALS'));
edd_r=find(strcmp(model_iJR_M3.rxns,'EDD'));
memb_r=find(strcmp(model_iJR_M3.rxns,'ATPS4r'));
co2_r=find(strcmp(model_iJR_M3.rxns,'EX_co2_e_'));

B(1) = 0.4741 ; B(2) = -0.0435;
gamma_v = [ 0.20,0.15];

for f=1:length(gamma_v)
    str = num2str(gamma_v(f));
    gamma_string = ['\gamma = ', str ];
     
    % Computation will take some seconds
    [gr,glc,ac,akgdh,mals,edd,memb,co2]=deal(zeros(Npoints,1));
    for i=1:Npoints
        model_iJR_M3.lb(ex_glc_r)=-wvec(i);
%         y = B(1)/wvec(i) + B(2);
%         if y < 0
%            y = 0; 
%         end
%         model_iJR_M3=setWeights(model_iJR_M3,7,y);

        sol_M3(i,f) = MAFBA_OptCbModel_cplex_v4(model_iJR_M3, 'gamma', gamma_v(f), 'cSense_m', 'U');
        
        gr(i)=sol_M3(i,f).x(bm_r);
        glc(i)=-sol_M3(i,f).x(glc_r); % note the minus sign
        ac(i) = sol_M3(i,f).x(ace_r);
        akgdh(i)=sol_M3(i,f).x(akgdh_r);
        mals(i)=sol_M3(i,f).x(mals_r);
        edd(i)=sol_M3(i,f).x(edd_r);
        memb(i)=0.5*sol_M3(i,f).x(memb_r);
        co2(i)=sol_M3(i,f).x(co2_r);

    end 
    T = table(wvec', gr, ac, glc, memb, akgdh, mals, edd , co2);
        T = sortrows(T,'gr');
        

        %Add detect acetate threshold >=0.1
        dif_v = diff(T.memb);
        dif_g = diff(T.gr);
        slope_v = dif_v./dif_g;
        m_s = max(slope_v);
        
        
        memb_s = find(slope_v <= m_s*0.01);
        t_i = find(T.ac(memb_s) >= 0.01, 1, 'first');
        ac_t = memb_s(t_i);

        if isempty(ac_t)
            fprintf('if 1 \n');
            memb_s = find(slope_v <= m_s*0.05);
            t_i = find(T.ac(memb_s) >= 0, 1, 'first');
            ac_t = memb_s(t_i);
        end
        if isempty(ac_t)
            fprintf('if 2 \n');
            memb_s = find(slope_v <= m_s*0.2);
            t_i = find(T.ac(memb_s) >= 0, 1, 'first');
            ac_t = memb_s(t_i);
        end
        
        if isempty(ac_t)
            fprintf('if 3 \n');
            ac_t = 1;
        end

        gr_ac=T.gr(ac_t);
    
    color = {sscanf('65AF65','%2x%2x%2x',[1 3])/255, sscanf('C74848','%2x%2x%2x',[1 3])/255, sscanf('9065D1','%2x%2x%2x',[1 3])/255,...
    sscanf('266DD7','%2x%2x%2x',[1 3])/255, sscanf('F2A1E4','%2x%2x%2x',[1 3])/255, sscanf('D2C725','%2x%2x%2x',[1 3])/255,...
    sscanf('B2B2B2','%2x%2x%2x',[1 3])/255,sscanf('EC8327','%2x%2x%2x',[1 3])/255};

    fig=figure(); hold on; 
    set(fig,'defaultAxesColorOrder',[color{4}; color{8}], 'Position',[10 10 800 600]);   
    title({'MAFBA: Flux profile vs growth rate','varying EX_{glc} reaction boundaries', gamma_string})
    
    yyaxis left
    set(gca,'FontSize',16)
    l1= plot(T.gr,T.glc,'-o','Color',color{1},'LineWidth',1, 'markerfacecolor', color{1}, 'MarkerSize', 8 );
    
    l3= plot(T.gr,T.akgdh,'-^','Color',color{3},'LineWidth',1, 'markerfacecolor', color{3}, 'MarkerSize', 8);
    l4= plot(T.gr,T.memb,'-d','Color',color{4}, 'LineWidth',1, 'markerfacecolor', color{4}, 'MarkerSize', 8);

    l5= plot(T.gr,T.mals,'-v','Color',color{5}, 'LineWidth',1, 'markerfacecolor', color{5}, 'MarkerSize', 8);
    l6= plot(T.gr,T.edd,'->','Color',color{6}, 'LineWidth',1, 'markerfacecolor', color{6}, 'MarkerSize', 8);
    l7= plot(T.gr,T.co2,'-<','Color',color{7}, 'LineWidth',1, 'markerfacecolor', color{7}, 'MarkerSize', 8);
    l2= plot(T.gr,T.ac,'-s','Color',color{2},'LineWidth',1, 'markerfacecolor', color{2}, 'MarkerSize', 8);
    y = ylim; 

    lac = plot([gr_ac,gr_ac], [0, y(2)],'--','Color',color{2},'LineWidth',1, 'markerfacecolor', color{2});
    xlabel('Growth rate (h^{-1})', 'Fontsize', 18);
    ylabel('Flux (mmol/g_{DW}h)', 'Color', color{4}, 'Fontsize', 18);

    yyaxis right
    l8= plot(T.gr,-T.Var1,'-+','Color',color{8},'LineWidth',1, 'markerfacecolor', color{8}, 'MarkerSize', 5 );
    ylabel('low boundary of EX_{glc} (mmol/g_{DW}h)', 'Color', color{8}, 'Fontsize', 18);
    lgd = legend([l1;l2;l3;l4;l5;l6;l7;l8;lac] ,'Glucose uptake','Acetate excretion','AKGDH (TCA)','ATPS4r flux 0.5x (Membrane)',...
        'Glyoxylate shunt ', 'ED pathway','CO_{2} excretion','lb of EX_{glc}', '\lambda_{ac}');
    set(lgd, 'Box', 'off', 'Location', 'northwest','Fontsize', 16);
    hold off;
end

