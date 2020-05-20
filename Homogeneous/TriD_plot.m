% load('AGZ_matrix_lb_no_wCc.mat');
load('AGY_matrix_lb_no_wCc.mat');
[z,g,v] = size(data);
data_short = data(:,:,1:3);
gr_b = data_short(:,1:36,3);
% gr_a = reshape(gr_b,[50,36]); % Z matrix
gr_a = reshape(gr_b,[60,36]); % Y matrix
gr_a(gr_a<0) = nan; 

gr_t = gr_a';

u_v = linspace(0,0.35, 36);
% gamma_v = [linspace(0.001,0.35,25), linspace(0.35,1,25) ];
gamma_v = [linspace(0.001,1,60) ];

% u_v = linspace(0,0.35, 60); % n vector -> X
% gamma_v = linspace(0,1,60); % m vector -> Y

y = u_v;
x = gamma_v;

[X,Y] = meshgrid(x,y);


fig=figure(); hold on; 
set(fig, 'Position',[10 10 800 600]);   
title({'Acetate threshold profile vs \gamma ratio'})
surf(X,Y,gr_t);
xlabel('\gamma ratio', 'Fontsize', 20);
zlabel('acetate threshold (\lambda_{ac})', 'Fontsize', 20);
ylabel('\phi Z', 'Fontsize', 20);
%% 
% color = {sscanf('65AF65','%2x%2x%2x',[1 3])/255, sscanf('C74848','%2x%2x%2x',[1 3])/255, sscanf('9065D1','%2x%2x%2x',[1 3])/255,...
% sscanf('266DD7','%2x%2x%2x',[1 3])/255, sscanf('F2A1E4','%2x%2x%2x',[1 3])/255, sscanf('D2C725','%2x%2x%2x',[1 3])/255,...
% sscanf('B2B2B2','%2x%2x%2x',[1 3])/255,sscanf('EC8327','%2x%2x%2x',[1 3])/255};
% fig=figure(); hold on; 
% set(fig,'defaultAxesColorOrder',[color{4}; color{8}], 'Position',[10 10 800 600]);   
% title({'Acetate threshold profile vs \gamma ratio'})
% l4= plot(x,gr_t(1,:),'-d','Color',color{1}, 'LineWidth',1, 'markerfacecolor', color{1}, 'MarkerSize', 5);
% 
% xlabel('\gamma ratio', 'Fontsize', 14);
% ylabel('acetate threshold (\lambda_{ac})', 'Color', color{4}, 'Fontsize', 14);
% hold off;
% 
% %% 
% color = {sscanf('65AF65','%2x%2x%2x',[1 3])/255, sscanf('C74848','%2x%2x%2x',[1 3])/255, sscanf('9065D1','%2x%2x%2x',[1 3])/255,...
% sscanf('266DD7','%2x%2x%2x',[1 3])/255, sscanf('F2A1E4','%2x%2x%2x',[1 3])/255, sscanf('D2C725','%2x%2x%2x',[1 3])/255,...
% sscanf('B2B2B2','%2x%2x%2x',[1 3])/255,sscanf('EC8327','%2x%2x%2x',[1 3])/255};
% fig=figure(); hold on; 
% set(fig,'defaultAxesColorOrder',[color{1}; color{1}], 'Position',[10 10 800 600]);   
% title({'Acetate threshold profile vs \phi Y'})
% l5= plot(y,gr_t(:,end),'-s','Color',color{8}, 'LineWidth',1, 'markerfacecolor', color{8}, 'MarkerSize', 5);
% 
% xlabel('\phi Z', 'Fontsize', 14);
% ylabel('acetate threshold (\lambda_{ac})', 'Fontsize', 14);
% hold off;

%% 
% 
% 
% fig=figure(); hold on; 
% set(fig, 'Position',[10 10 800 600]);   
% title({'Acetate threshold profile vs \gamma ratio and \phi Z'})
% 
% yyaxis left
% set(gca,'FontSize',12)
% l1= plot(gr,glc,'-o','Color',color{1},'LineWidth',1, 'markerfacecolor', color{1}, 'MarkerSize', 5 );
% 
% l3= plot(gr,akgdh,'-^','Color',color{3},'LineWidth',1, 'markerfacecolor', color{3}, 'MarkerSize', 5);
% l4= plot(gr,memb,'-d','Color',color{4}, 'LineWidth',1, 'markerfacecolor', color{4}, 'MarkerSize', 5);
% 
% l5= plot(gr,mals,'-v','Color',color{5}, 'LineWidth',1, 'markerfacecolor', color{5}, 'MarkerSize', 5);
% l6= plot(gr,edd,'->','Color',color{6}, 'LineWidth',1, 'markerfacecolor', color{6}, 'MarkerSize', 5);
% l7= plot(gr,co2,'-<','Color',color{7}, 'LineWidth',1, 'markerfacecolor', color{7}, 'MarkerSize', 5);
% l2= plot(gr,ac,'-s','Color',color{2},'LineWidth',1, 'markerfacecolor', color{2}, 'MarkerSize', 7);
% 
% xlabel('Growth rate (h^{-1})', 'Fontsize', 14);
% ylabel('Flux (mmol/g_{DW}h)', 'Color', color{4}, 'Fontsize', 14);
% 
% yyaxis right
% l8= plot(gr,-wvec,'-+','Color',color{8},'LineWidth',1, 'markerfacecolor', color{8}, 'MarkerSize', 5 );
% ylabel('low boundary of EX_{glc} (mmol/g_{DW}h)', 'Color', color{8}, 'Fontsize', 14);
% lgd = legend([l1;l2;l3;l4;l5;l6;l7;l8] ,'Glucose uptake','Acetate excretion','AKGDH (TCA)','ATPS4r flux (Membrane)',...
%     'Glyoxylate shunt ', 'ED pathway','CO_{2} excretion','lb of EX_{glc}');
% set(lgd, 'Box', 'on', 'Location', 'northwest','Fontsize', 10);
% hold off;