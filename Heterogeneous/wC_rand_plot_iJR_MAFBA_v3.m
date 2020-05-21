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
model_iJR_M3.protGroup(7).phi0 = 0;
weight = 0.0007596;
% weight = 0.00083;
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

% This is a good wC mesh
Npoints=100;
Wpoints=1000;
wvec= (0.01*linspace(0,1,Npoints)  + 0.99*linspace(0,1,Npoints) .^3 );
x_matrix = zeros(length(model_iJR_M3.rxns),Wpoints);
sol_arr = zeros(length(model_iJR_M3.rxns),Npoints);

ex_glc_r=find(strcmp(model_iJR_M3.rxns,'EX_glc_e_'));
bm_r=model_iJR_M3.protGroup(3).rxns;
akgdh_r=find(strcmp(model_iJR_M3.rxns,'AKGDH'));
ace_r=find(strcmp(model_iJR_M3.rxns,'EX_ac_e_'));
mals_r=find(strcmp(model_iJR_M3.rxns,'MALS'));
edd_r=find(strcmp(model_iJR_M3.rxns,'EDD'));
memb_r=find(strcmp(model_iJR_M3.rxns,'ATPS4r'));
co2_r=find(strcmp(model_iJR_M3.rxns,'EX_co2_e_'));

color = {sscanf('65AF65','%2x%2x%2x',[1 3])/255, sscanf('C74848','%2x%2x%2x',[1 3])/255, sscanf('9065D1','%2x%2x%2x',[1 3])/255,...
sscanf('266DD7','%2x%2x%2x',[1 3])/255, sscanf('F2A1E4','%2x%2x%2x',[1 3])/255, sscanf('D2C725','%2x%2x%2x',[1 3])/255,...
sscanf('B2B2B2','%2x%2x%2x',[1 3])/255,sscanf('EC8327','%2x%2x%2x',[1 3])/255};

gamma_v = [0.25];

%% 

for f=1:length(gamma_v)
    str = num2str(gamma_v(f));
    gamma_string = ['\gamma = ', str ];
        % Computation will take aprox one hour
    for i=1:Npoints
        model_iJR_M3=setWeights(model_iJR_M3,7,wvec(i));
        avg_w = 0.000822;
        for j=1:Wpoints
            model_iJR_M3=setWeights_rand(model_iJR_M3,1,'expbox',avg_w,1);
            model_iJR_M3=setWeights_rand(model_iJR_M3,2,'expbox',avg_w,1);
            model_iJR_M3=setWeights_rand(model_iJR_M3,4,'expbox',avg_w,1);
            model_iJR_M3=setWeights_rand(model_iJR_M3,5,'expbox',avg_w,1);
            sol = MAFBA_OptCbModel_cplex_v3(model_iJR_M3, 'gamma', gamma_v(f), 'cSense_m', 'U');
            x_matrix(:,j) = sol.x;
        end
        x_avg = mean(x_matrix,2);
        sol_arr(:,i) = x_avg;
        fprintf('%i--%i\n', f,i)
        str = strrep(str,'.','');
        save(['rand_sol_arr_g_', str,'.mat'],'sol_arr');

    end    
    gr=sol_arr(bm_r,:);
    glc=-sol_arr(ex_glc_r,:); % note the minus sign
    ac = sol_arr(ace_r,:);
    akgdh =sol_arr(akgdh_r,:);
    memb =sol_arr(memb_r,:);
    mals =sol_arr(mals_r,:);
    edd =sol_arr(edd_r,:);
    co2 =sol_arr(co2_r,:);
    
    T = table(wvec, gr, ac, glc, memb, akgdh, mals, edd , co2);
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
    
    
%     fig=figure(); hold on; 
%     set(fig,'defaultAxesColorOrder',[color{4}; color{8}], 'Position',[10 10 780 560]);
%     
%     title({'MAFBA: Flux profile vs growth rate','Heterogeneous weights',['varying w_{Cc}', '   ' ,gamma_string]})
%     yyaxis left
%     set(gca,'FontSize',18)
%        
% 
%     l1= plot(T.gr,T.glc,'-o','Color',color{1},'LineWidth',1, 'markerfacecolor', color{1}, 'MarkerSize', 5 );
%     
%     l3= plot(T.gr,T.akgdh,'-^','Color',color{3},'LineWidth',1, 'markerfacecolor', color{3}, 'MarkerSize', 5);
%     l4= plot(T.gr,0.5*T.memb,'-d','Color',color{4}, 'LineWidth',1, 'markerfacecolor', color{4}, 'MarkerSize', 5);
% 
%     l5= plot(T.gr,T.mals,'-v','Color',color{5}, 'LineWidth',1, 'markerfacecolor', color{5}, 'MarkerSize', 5);
%     l6= plot(T.gr,T.edd,'->','Color',color{6}, 'LineWidth',1, 'markerfacecolor', color{6}, 'MarkerSize', 5);
%     l7= plot(T.gr,T.co2,'-<','Color',color{7}, 'LineWidth',1, 'markerfacecolor', color{7}, 'MarkerSize', 5);
%     l2= plot(T.gr,T.ac,'-s','Color',color{2},'LineWidth',1, 'markerfacecolor', color{2}, 'MarkerSize', 7);
%     y = ylim;
%     lac = plot([gr_ac,gr_ac], [0, y(2)],'--','Color',color{2},'LineWidth',1); 
%     
%     xlabel('Growth rate (h^{-1})', 'Fontsize', 18);
%     ylabel('Flux (mmol/g_{DW}h)', 'Color', color{4}, 'Fontsize', 18);
%     axis([0 inf 0 inf])
% 
%     yyaxis right
%     l8= plot(T.gr,T.wvec,'-+','Color',color{8},'LineWidth',1, 'markerfacecolor', color{8}, 'MarkerSize', 5 );
%     ylabel('Weight (gh/mmol)', 'Color', color{8}, 'Fontsize', 18);
%     lgd = legend([l1;l2;l3;l4;l5;l6;l7;l8;lac] ,'Glucose uptake','Acetate excretion','AKGDH (TCA)',...
%     'ATPS4r flux 0.5x (Membrane)','Glyoxylate shunt ', 'ED pathway','CO_{2} excretion','w_{C}',' \lambda_{ac}');
%     set(lgd, 'Box', 'off', 'Location', 'northwest','Fontsize', 16,...
%         'Position',[0.164529914529915 0.532440503084436 0.292307685544858 0.277678563552244]);
%     
%     axis([0 inf 0 inf])
%     hold off;
%     
%     filename = ['MAFBA_C_rand_g_', str,'.png'];
%     
%     saveas(fig,filename,'png');

end
    


