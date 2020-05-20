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
glc_r=find(strcmp(model_iJR_M3.rxns,'EX_glc_e_'));
bm_r=model_iJR_M3.protGroup(3).rxns;
akgdh_r=find(strcmp(model_iJR_M3.rxns,'AKGDH'));
ace_r=find(strcmp(model_iJR_M3.rxns,'EX_ac_e_'));
mals_r=find(strcmp(model_iJR_M3.rxns,'MALS'));
edd_r=find(strcmp(model_iJR_M3.rxns,'EDD'));
memb_r=find(strcmp(model_iJR_M3.rxns,'ATPS4r'));
co2_r=find(strcmp(model_iJR_M3.rxns,'EX_co2_e_'));

B(1) = 0.4741 ; B(2) = -0.0435;
u_v = linspace(0,0.35, 36);
u = 1;
% gamma_v = [linspace(0.001,0.35,25), linspace(0.35,1,25) ];
gamma_v = [linspace(0.001,1,60) ];

data = deal(zeros(length(u_v),length(gamma_v),7));


for u=1:length(u_v)
    
z_r = find(strcmp(model_iJR_M3.rxns,'RxnZ'));
model_iJR_M3.lb(z_r) = 1;
model_iJR_M3.ub(z_r) = 1;
model_iJR_M3.w(z_r)=u_v(1); % percentage of fraction Z

y_r = find(strcmp(model_iJR_M3.rxns,'RxnY'));
model_iJR_M3.lb(y_r) = 1;
model_iJR_M3.ub(y_r) = 1;
model_iJR_M3.w(y_r)=u_v(u); % percentage of fraction Y

    for f=1:length(gamma_v)
        str = num2str(gamma_v(f));
        gamma_string = ['\gamma = ', str ];
        
        fprintf("%1.1i - %1.1i \n",u,f);

        % Computation will take some seconds
        [gr,glc,ac,akgdh,mals,edd,memb,co2]=deal(zeros(Npoints,1));
        for i=1:Npoints
            model_iJR_M3.lb(ex_glc_r)=-wvec(i);
%             y = B(1)/wvec(i) + B(2);
%             if y < 0
%                y = 0; 
%             end
%             model_iJR_M3=setWeights(model_iJR_M3,7,y);

            sol_M3(i,f) = MAFBA_OptCbModel_cplex_v4(model_iJR_M3, 'gamma', gamma_v(f), 'cSense_m', 'U');

            gr(i)=sol_M3(i,f).x(bm_r);
            glc(i)=-sol_M3(i,f).x(glc_r); % note the minus sign
            ac(i) = sol_M3(i,f).x(ace_r);
            akgdh(i)=sol_M3(i,f).x(akgdh_r);
            mals(i)=sol_M3(i,f).x(mals_r);
            edd(i)=sol_M3(i,f).x(edd_r);
            memb(i)=sol_M3(i,f).x(memb_r);
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
        
        % I need another way to detect the acetate threshold, maybe when
        % the diference in the change of memb flux is negative.        
        data(f,u,1) = u_v(u); % Z%
        data(f,u,2) = gamma_v(f);
        data(f,u,3) = T.gr(ac_t);
        data(f,u,4) = T.ac(ac_t); 
        data(f,u,5) = T.memb(ac_t);
        data(f,u,6) = T.glc(ac_t);
        data(f,u,7) = sol_M3(i,f).phiM;

%        
    end
    
end

filename = 'AGY_matrix_lb_no_wCc.mat';
save(filename, 'data');
