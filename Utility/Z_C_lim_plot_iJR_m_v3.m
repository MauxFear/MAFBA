
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
wvec= (linspace(0.09124203971,1,Npoints) .^3 );
ex_glc_r=find(strcmp(model_iJR_M3.rxns,'EX_glc_e_'));

ex_c = model_iJR_M3.protGroup(3).rxns;
glc_r=find(strcmp(model_iJR_M3.rxns,'EX_glc_e_'));
bm_r=model_iJR_M3.protGroup(3).rxns;
akgdh_r=find(strcmp(model_iJR_M3.rxns,'AKGDH'));
ace_r=find(strcmp(model_iJR_M3.rxns,'EX_ac_e_'));
mals_r=find(strcmp(model_iJR_M3.rxns,'MALS'));
edd_r=find(strcmp(model_iJR_M3.rxns,'EDD'));
memb_r=find(strcmp(model_iJR_M3.rxns,'ATPS4r'));
co2_r=find(strcmp(model_iJR_M3.rxns,'EX_co2_e_'));

z_r = find(strcmp(model_iJR_M3.rxns,'RxnZ'));

% B(1) = 0.4741 ; B(2) = -0.0435;

%Reactions to simulate the expression of a cytosolic protein Z and a
%membrane protein Y, changing the percentage of fraction Z or Y we are
%simulating protein overexpression. 

model_iJR_M3.lb(z_r) = 1;
model_iJR_M3.ub(z_r) = 1;
model_iJR_M3.w(z_r)=0.08; % percentage of fraction Z

y_r = find(strcmp(model_iJR_M3.rxns,'RxnY'));
model_iJR_M3.lb(y_r) = 1;
model_iJR_M3.ub(y_r) = 1;
model_iJR_M3.w(y_r)=0; % percentage of fraction Y

% model_iJR_M3.protGroup(6).phi0 = 0.45 + 0.16; % Increasing Q fraction

gamma_v = [0.20];
% gamma_v = [0.15,0.10,0.05];


for f=1:length(gamma_v)
    str = num2str(gamma_v(f));
    gamma_string = ['\gamma = ', str , '    \phi Z = 8%'];
     
    % Computation will take some seconds
    [gr,glc,ac,akgdh,mals,edd,memb,co2]=deal(zeros(Npoints,1));
    for i=1:Npoints
        % Here we progressively increase the weight of the Cc group
        % (therefore reducing the glucose uptake rate).
        model_iJR_M3=setWeights(model_iJR_M3,7,wvec(i));
        sol_M(i,f) = MAFBA_OptCbModel_cplex_v4(model_iJR_M3, 'gamma', gamma_v(f), 'cSense_m', 'U');
        
        gr(i)=sol_M(i,f).x(bm_r);
        glc(i)=-sol_M(i,f).x(glc_r); % note the minus sign
        ac(i) = sol_M(i,f).x(ace_r);
        akgdh(i)=sol_M(i,f).x(akgdh_r);
        mals(i)=sol_M(i,f).x(mals_r);
        edd(i)=sol_M(i,f).x(edd_r);
        memb(i)=sol_M(i,f).x(memb_r);
        co2(i)=sol_M(i,f).x(co2_r);

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
    
    
    color = {sscanf('65AF65','%2x%2x%2x',[1 3])/255, sscanf('C74848','%2x%2x%2x',[1 3])/255, sscanf('9065D1','%2x%2x%2x',[1 3])/255,...
    sscanf('266DD7','%2x%2x%2x',[1 3])/255, sscanf('F2A1E4','%2x%2x%2x',[1 3])/255, sscanf('D2C725','%2x%2x%2x',[1 3])/255,...
    sscanf('B2B2B2','%2x%2x%2x',[1 3])/255,sscanf('EC8327','%2x%2x%2x',[1 3])/255};

    fig=figure(); hold on; 
    set(fig,'defaultAxesColorOrder',[color{4}; color{8}], 'Position',[10 10 800 600]);   
    title({'MAFBA: Flux profile vs growth rate','varying EX_{glc} reaction boundaries', gamma_string})
    
    yyaxis left
    set(gca,'FontSize',12);
    l1= plot(T.gr,T.glc,'-o','Color',color{1},'LineWidth',1, 'markerfacecolor', color{1}, 'MarkerSize', 5 );
    
    l3= plot(T.gr,T.akgdh,'-^','Color',color{3},'LineWidth',1, 'markerfacecolor', color{3}, 'MarkerSize', 5);
    l4= plot(T.gr,T.memb,'-d','Color',color{4}, 'LineWidth',1, 'markerfacecolor', color{4}, 'MarkerSize', 5);

    l5= plot(T.gr,T.mals,'-v','Color',color{5}, 'LineWidth',1, 'markerfacecolor', color{5}, 'MarkerSize', 5);
    l6= plot(T.gr,T.edd,'->','Color',color{6}, 'LineWidth',1, 'markerfacecolor', color{6}, 'MarkerSize', 5);
    l7= plot(T.gr,T.co2,'-<','Color',color{7}, 'LineWidth',1, 'markerfacecolor', color{7}, 'MarkerSize', 5);
    l2= plot(T.gr,T.ac,'-s','Color',color{2},'LineWidth',1, 'markerfacecolor', color{2}, 'MarkerSize', 7);
    y = ylim; 

    lac = plot([gr_ac,gr_ac], [0, y(2)],'--','Color',color{2},'LineWidth',1, 'markerfacecolor', color{2});
    xlabel('Growth rate (h^{-1})', 'Fontsize', 14);
    ylabel('Flux (mmol/g_{DW}h)', 'Color', color{4}, 'Fontsize', 14);
    axis([0 inf 0 inf]);
%     axis([0.8 inf 0 10])

    yyaxis right
    l8= plot(T.gr,T.Var1,'-+','Color',color{8},'LineWidth',1, 'markerfacecolor', color{8}, 'MarkerSize', 5 );
    ylabel('low boundary of EX_{glc} (mmol/g_{DW}h)', 'Color', color{8}, 'Fontsize', 14);
    lgd = legend([l1;l2;l3;l4;l5;l6;l7;l8;lac] ,'Glucose uptake','Acetate excretion','AKGDH (TCA)','ATPS4r flux (Membrane)',...
        'Glyoxylate shunt ', 'ED pathway','CO_{2} excretion','lb of EX_{glc}', '\lambda_{ac}');
    set(lgd, 'Box', 'off', 'Location', 'west','Fontsize', 9,...
        'Position',[0.125214521452145 0.519773720831774 0.303217827211512 0.270695370396242]);
    axis([0 inf 0 inf]);
%     axis([0.8 inf -inf 0])

     
    hold off;
    filename = ['MAFBA_Z_08_g_', str,'.png'];
    saveas(fig,filename,'png');
    
end

