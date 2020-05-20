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
model_iJR_M3.protGroup(7).phi0 = 0; % phiQ_0

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
wvec= (0.1*linspace(0,1,Npoints)  + 0.99*linspace(0,1,Npoints) .^3 );
% wvec= (70*linspace(0,1,Npoints) .^3 );

bm_r=model_iJR_M3.protGroup(3).rxns;
gamma_v = [1];
ex_glc_r=find(strcmp(model_iJR_M3.rxns,'EX_glc_e_'));


for f=1:length(gamma_v)
    str = num2str(gamma_v(f));
    gamma_string = ['\gamma = ', str ];
     
    % Computation will take some seconds
    [gr,phiCm,phiR,phiEc,phiEr,phiM,phi_max,phiEm,phiCc]=deal(zeros(Npoints,1));
    
    for i=1:Npoints
        model_iJR_M3=setWeights(model_iJR_M3,7,wvec(i));
        sol_M3(i,f) = MAFBA_OptCbModel_cplex_v4(model_iJR_M3, 'gamma', gamma_v(f), 'cSense_m', 'U');
        gr(i)=sol_M3(i,f).x(bm_r);
        phiCm(i)=sol_M3(i,f).phiCm;
        phiCc(i)=sol_M3(i,f).phiCc;
        phiR(i)=sol_M3(i,f).phiR; 
        phiEc(i) = sol_M3(i,f).phiEc;
        phiEr(i)=sol_M3(i,f).phiEr;
        phiEm(i)=sol_M3(i,f).phiEm;
        phiM(i)=sol_M3(i,f).phiM;
        phi_max(i)=sol_M3(i,f).phi_max;

    end 
    
        
    color = {sscanf('65AF65','%2x%2x%2x',[1 3])/255, sscanf('C74848','%2x%2x%2x',[1 3])/255, sscanf('9065D1','%2x%2x%2x',[1 3])/255,...
    sscanf('266DD7','%2x%2x%2x',[1 3])/255, sscanf('F2A1E4','%2x%2x%2x',[1 3])/255, sscanf('D2C725','%2x%2x%2x',[1 3])/255,...
    sscanf('B2B2B2','%2x%2x%2x',[1 3])/255,sscanf('EC8327','%2x%2x%2x',[1 3])/255, sscanf('35B2B2','%2x%2x%2x',[1 3])/255};

    fig=figure(); hold on; 
    set(fig,'Position',[10 10 800 600]); 

    title({'MAFBA: Sectors profile vs growth rate','varying w_{C}', gamma_string})
%     yyaxis right
    
    set(gca,'FontSize',16)
    

    l2= plot(gr,phiR,'-^','Color',color{3},'LineWidth',1, 'markerfacecolor', color{3}, 'MarkerSize', 5);

    l4= plot(gr,phiEr,'-d','Color',color{4}, 'LineWidth',1, 'markerfacecolor', color{4}, 'MarkerSize', 5);
    l6= plot(gr,phiM,'->','Color',color{9}, 'LineWidth',1, 'markerfacecolor', color{9}, 'MarkerSize', 5);
    
    l5= plot(gr,phiEm,'-v','Color',color{6}, 'LineWidth',1, 'markerfacecolor', color{6}, 'MarkerSize', 5);
    l1= plot(gr,phiCm,'-o','Color',color{1},'LineWidth',1, 'markerfacecolor', color{1}, 'MarkerSize', 5 );
    
    xlabel('Growth rate (h^{-1})', 'Fontsize', 18);
%     ylabel('Proteome fraction (\phi R, \phi E_{r}, \phi E_{m} and \phi^{m})', 'Color', color{8}, 'Fontsize', 14);
%     axis([0 inf 0 0.2])
    
%     yyaxis left
    l8= plot(gr,phiCc,'-o','Color',color{8},'LineWidth',1, 'markerfacecolor', color{8}, 'MarkerSize', 8 );

    l3= plot(gr,phiEc,'-s','Color',color{2},'LineWidth',1, 'markerfacecolor', color{2}, 'MarkerSize', 8);
    l7= plot(gr,phi_max,'-x','Color',color{7}, 'LineWidth',1, 'markerfacecolor', color{7}, 'MarkerSize', 8);

    ylabel('Proteome fraction', 'Fontsize', 18);
    lgd = legend([l1;l2;l3;l4;l5;l8;l6;l7] ,'\phi C_{m}','\phi R','\phi E_{c}', '\phi E_{r}', '\phi E_{m}','\phi C_{c}', '\phi^{m}','\phi_{max}');
    set(lgd, 'Box', 'off', 'Location', 'northwest','Fontsize', 18);
    axis([0 inf 0 0.5])
    hold off;
    
    filename = ['MAFBA_C_lim_sec_g_', str,'.png'];
    saveas(fig,filename,'png');
    
end

