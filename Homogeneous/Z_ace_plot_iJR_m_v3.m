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
color = {sscanf('65AF65','%2x%2x%2x',[1 3])/255, sscanf('C74848','%2x%2x%2x',[1 3])/255, sscanf('9065D1','%2x%2x%2x',[1 3])/255,...
    sscanf('266DD7','%2x%2x%2x',[1 3])/255, sscanf('F2A1E4','%2x%2x%2x',[1 3])/255, sscanf('D2C725','%2x%2x%2x',[1 3])/255,...
    sscanf('B2B2B2','%2x%2x%2x',[1 3])/255,sscanf('EC8327','%2x%2x%2x',[1 3])/255};

Npoints=100;
wvec= (0.4+ 5* linspace(0,1,Npoints)  + 20.5*linspace(0,1,Npoints) .^3 );
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
y_r = find(strcmp(model_iJR_M3.rxns,'RxnY'));

B(1) = 0.4741 ; B(2) = -0.0435;

gamma_v = [1];

Rxn_U = 'z';
u_v = [0, 0.02, 0.04, 0.10, 0.16]; % Z
color = {color{1},color{2},color{4},color{3},color{7},color{8},color{5}}; % Z
markers = {'-o', '-^', '-d', '-s', '-v', '-<'}; % Z

% Rxn_U = 'y';
% u_v = [0, 0.01, 0.02, 0.04];  % Y
% color = {color{1},color{8},color{2},color{4},color{3},color{7}}; %Y
% markers = {'-o', '->', '-^', '-d', '-s', '-v', '-<'}; % Y



for f=1:length(gamma_v)
    str = num2str(gamma_v(f));
    fig=figure(); 
    
    set(fig, 'Position',[10 10 800 600]);
    gamma_string = ['\gamma = ', str , '    \phi Z ' ];
    title({'MAFBA: Flux profile vs growth rate','varying EX_{glc} reaction boundaries', gamma_string})
        
        
    for u=1:length(u_v)
        if u == 1
            hold on
        end
        num_u = u_v(u)*100;
        utr = num2str(num_u);
        
        if Rxn_U == 'z'
            model_iJR_M3.lb(z_r) = 1;
            model_iJR_M3.ub(z_r) = 1;
            model_iJR_M3.w(z_r)= u_v(u); % percentage of fraction Z

            model_iJR_M3.lb(y_r) = 1;
            model_iJR_M3.ub(y_r) = 1;
            model_iJR_M3.w(y_r)= 0; % percentage of fraction Y
        
        elseif Rxn_U == 'y'
            model_iJR_M3.lb(z_r) = 1;
            model_iJR_M3.ub(z_r) = 1;
            model_iJR_M3.w(z_r)= 0; % percentage of fraction Z

            y_r = find(strcmp(model_iJR_M3.rxns,'RxnY'));
            model_iJR_M3.lb(y_r) = 1;
            model_iJR_M3.ub(y_r) = 1;
            model_iJR_M3.w(y_r)= u_v(u); % percentage of fraction Y
        end

        % Computation will take some seconds
        [gr,glc,ac,akgdh,mals,edd,memb,co2]=deal(zeros(Npoints,1));
        for i=1:Npoints
%             for r=1:length(ex_c)
%                 model_iJR_M3.lb(ex_c(r))=-wvec(i);
%             end
              model_iJR_M3.lb(ex_glc_r)=-wvec(i); 
%             y = B(1)/wvec(i) + B(2);
%             if y < 0
%                y = 0; 
%             end
%             model_iJR_M3=setWeights(model_iJR_M3,7,y);
            sol_M(i,f) = MAFBA_OptCbModel_cplex_v4(model_iJR_M3, 'gamma', gamma_v(f), 'cSense_m', 'U');

            gr(i)=sol_M(i,f).x(bm_r);
            glc(i)=-sol_M(i,f).x(glc_r); % note the minus sign
            ac(i) = sol_M(i,f).x(ace_r);
            akgdh(i)=sol_M(i,f).x(akgdh_r);
            mals(i)=sol_M(i,f).x(mals_r);
            edd(i)=sol_M(i,f).x(edd_r);
            memb(i)=0.5*sol_M(i,f).x(memb_r);
            co2(i)=sol_M(i,f).x(co2_r);
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
        
        str_u = ['\phi Z ', utr, '%']; 
        plot(T.gr,T.ac,markers{u},'Color',color{u},'LineWidth',1, 'markerfacecolor', color{u}, 'MarkerSize', 8, 'DisplayName',str_u );
        plot([gr_ac,gr_ac], [0, 100],'--','Color',color{u},'LineWidth',1, 'HandleVisibility','off');
%         plot(gr,memb,markers{u},'Color',color{u},'LineWidth',1, 'markerfacecolor', color{u}, 'MarkerSize', 5 , 'DisplayName',[str_u,' membrane (0.5x)']);
        max_ac(u) = max(ac);
%         max_memb(u) = max(memb);
        
        
    end   
    
     set(gca,'FontSize',18);
        
    xlabel('Growth rate (h^{-1})', 'Fontsize', 20);
    ylabel('Flux (mmol/g_{DW}h)', 'Fontsize', 20);
    plot([0,0], [100, 100],'--','Color','k','LineWidth',1, 'DisplayName', '\lambda_{ac}' );
    
    lgd = legend();
    set(lgd, 'Box', 'off', 'Location', 'west','Fontsize', 18,...
        'Position',[0.127714521452145 0.273038080643847 0.303217827211512 0.547499984105428]);
    acmax = max(max_ac);
    membmax = 0;
%     membmax = max(max_memb);
    
    if membmax > acmax
        max_y = membmax + 0.5;
    else 
        max_y = max(max_ac) + 0.5;
    end
    
    axis([0 inf 0 max_y]);

    str = strrep(str,'.',''); 
    hold off;
    filename = ['MAFBA_Z_lb_ace_g_', str,'.png'];
    saveas(fig,filename,'png');
    
end

