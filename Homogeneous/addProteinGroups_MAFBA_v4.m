function model=addProteinGroups_MAFBA_v4(model,optC)
% model=addProteinGroups_MAFBA(model,optC)
%
% This function adds the 'protGroup' field to a COBRA model including the
% membrane allocation groups.
% Inputs:
%    model : COBRA metabolic model
%    optC  : Chooses the content of the phiC group. Here we report some
%            examples:
%            {'default'} (equivalent to 'glc')
%            'glc'
%            'lcts'
%            'fru'
%            'akg'
%            'g6p'
%            'g6p+glcn'
%           The full list can be found in the setPhiC() subfunction.
% Output:
%    model : same model as input, with the additional 'protGroup' field
%
%

if ~exist('optC','var')
    optC='glc';
end
phiEopt='subSys';

% Reactions
Biomass=find(model.c);

if ~isfield(model,'protGroup')
    
    protGroup(3).name='phiR';
    protGroup(3).phi0=0.066;
    protGroup(3).rxns=Biomass;
     
    protGroup(1).name='phiCm';
    protGroup(1).phi0=0;
    [model,rxns_Cm,rxns_Cc] = setPhiC(model,optC);
    protGroup(1).rxns=rxns_Cm;
    protGroup(7).name='phiCc';
    protGroup(7).phi0=0;
    protGroup(7).rxns=rxns_Cc;
    protGroup(7).rxns=removeRxnsOverlap(protGroup(7).rxns,protGroup(1).rxns);
    
    
    % Defines all protGroup entries
    protGroup(2).name='phiEc';
    protGroup(4).name='phiEr';
    protGroup(5).name='phiEm';
    protGroup(2).phi0=0;
    protGroup(4).phi0=0;
    protGroup(5).phi0=0;
    [protGroup(2).rxns , protGroup(4).rxns, protGroup(5).rxns ]=set_PhiE_list(model,phiEopt);
    protGroup(2).rxns=removeRxnsOverlap(protGroup(2).rxns,protGroup(1).rxns);
    protGroup(2).rxns=removeRxnsOverlap(protGroup(2).rxns,protGroup(7).rxns);
    protGroup(2).rxns=removeRxnsOverlap(protGroup(2).rxns,protGroup(3).rxns);
    protGroup(4).rxns=removeRxnsOverlap(protGroup(4).rxns,protGroup(1).rxns);
    protGroup(4).rxns=removeRxnsOverlap(protGroup(4).rxns,protGroup(2).rxns);
    protGroup(4).rxns=removeRxnsOverlap(protGroup(4).rxns,protGroup(3).rxns);
    protGroup(4).rxns=removeRxnsOverlap(protGroup(4).rxns,protGroup(7).rxns);
    protGroup(5).rxns=removeRxnsOverlap(protGroup(5).rxns,protGroup(1).rxns);
    protGroup(5).rxns=removeRxnsOverlap(protGroup(5).rxns,protGroup(2).rxns);
    protGroup(5).rxns=removeRxnsOverlap(protGroup(5).rxns,protGroup(3).rxns);
    protGroup(5).rxns=removeRxnsOverlap(protGroup(5).rxns,protGroup(7).rxns);
    protGroup(5).rxns=removeRxnsOverlap(protGroup(5).rxns,protGroup(4).rxns);
    
    
    
    protGroup(6).name='phiQ and others'; % weights=0
    protGroup(6).phi0=0;
    protGroup(6).rxns=(1:length(model.rxns))';
    protGroup(6).rxns=removeRxnsOverlap(protGroup(6).rxns,protGroup(1).rxns);
    protGroup(6).rxns=removeRxnsOverlap(protGroup(6).rxns,protGroup(2).rxns);
    protGroup(6).rxns=removeRxnsOverlap(protGroup(6).rxns,protGroup(3).rxns);
    protGroup(6).rxns=removeRxnsOverlap(protGroup(6).rxns,protGroup(4).rxns);
    protGroup(6).rxns=removeRxnsOverlap(protGroup(6).rxns,protGroup(5).rxns);
    protGroup(6).rxns=removeRxnsOverlap(protGroup(6).rxns,protGroup(7).rxns);
    
    model.protGroup=protGroup;
else
    % only updates the first prot group
    [model,rxns_Cm,rxns_Cc] = setPhiC(model,optC);
    model.protGroup(1).rxns=rxns_Cm;
    model.protGroup(7).rxns=rxns_Cc;
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION - Set the PhiC rxns vector and bounds %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [model,rxns_Cm, rxns_Cc] = setPhiC(model,optC)
% This function also updates the bound for the exchange reactions

EX_glc=find(strcmp(model.rxns,'EX_glc_e_')); % glc transport
EX_biomass=find(strcmp(model.rxns,'EX_biomass_c')); % glc transport

% Upper glycolysis
if strcmp(optC,'default') || strcmp(optC,'glc')
    rxns_Cm=find(ismember(model.rxns,{'GLCpts','GLCt2'} ));
    rxns_Cc=EX_biomass;
elseif strcmp(optC,'g6p')
    EX_g6p=find(strcmp(model.rxns,'EX_g6p_e_')); % g6p transport
    t_g6p=find(strcmp(model.rxns,'G6Pt6_2')); % g6p transport
    rxns_Cm=t_g6p;
    rxns_Cc=EX_g6p;
    model.lb(EX_g6p)=-1000;
    model.lb(EX_glc)=0; 
elseif strcmp(optC,'man6p')
    EX_man6p=find(strcmp(model.rxns,'EX_man6p_e_')); % man6p transport
    t_man6p=find(strcmp(model.rxns,'MAN6Pt6_2')); % man6p transport
    rxns_Cc=EX_man6p;
    rxns_Cm=t_man6p;
    model.lb(EX_man6p)=-1000;
    model.lb(EX_glc)=0; 
elseif strcmp(optC,'mnl')
    EX_mnl=find(strcmp(model.rxns,'EX_mnl_e_')); % mannitol transport
    t_mnl=find(strcmp(model.rxns,'MNLpts')); % mannitol transport
    rxns_Cm=t_mnl;
    rxns_Cc=EX_mnl;
    model.lb(EX_mnl)=-1000;
    model.lb(EX_glc)=0;
elseif strcmp(optC,'lcts')
    EX_lcts=find(strcmp(model.rxns,'EX_lcts_e_') ); % lactose transport
    t_lcts=find(strcmp(model.rxns,'LCTSt') ); % lactose transport
    rxns_Cc=EX_lcts;
    rxns_Cm=t_lcts;
    model.lb(EX_lcts)=-1000;
    model.lb(EX_glc)=0;
elseif strcmp(optC,'sucr')
    EX_sucr=find(strcmp(model.rxns,'EX_sucr_e_') ); % sucrose transport
    t_sucr=find(strcmp(model.rxns,'SUCpts') ); % lactose transport
    rxns_Cc=EX_sucr;
    rxns_Cm=t_sucr;
    model.lb(EX_sucr)=-1000;
    model.lb(EX_glc)=0;
elseif strcmp(optC,'malt')
    EX_malt=find(strcmp(model.rxns,'EX_malt_e_') ); % maltose transport
    t_malt=find(ismember(model.rxns,{'MALTabc','MALTpts'} )); % maltose transport
    rxns_Cc=EX_malt;
    rxns_Cm=t_malt;
    model.lb(EX_malt)=-1000;
    model.lb(EX_glc)=0;
elseif strcmp(optC,'fru')
    EX_fru=find(strcmp(model.rxns,'EX_fru_e_')); % fructose transport
    t_fru=find(ismember(model.rxns,{'FRUpts','FRUpts2'} )); % fructose transport
    rxns_Cc=EX_fru;     rxns_Cm=t_fru;
    model.lb(EX_fru)=-1000;
    model.lb(EX_glc)=0;
elseif strcmp(optC,'sbt')
    EX_sbt=find(strcmp(model.rxns,'EX_sbt_D_e_')); % sorbitol
    t_sbt=find(strcmp(model.rxns,'SBTpts')); % sorbitol
    rxns_Cc=EX_sbt;
    rxns_Cm=t_sbt;
    model.lb(EX_sbt)=-1000;
    model.lb(EX_glc)=0;
elseif strcmp(optC,'man')
    EX_man=find(strcmp(model.rxns,'EX_man_e_')); % mannose
    t_man=find(strcmp(model.rxns,'MANpts')); % mannose
    rxns_Cc=EX_man;
    rxns_Cm=t_man;
    model.lb(EX_man)=-1000;
    model.lb(EX_glc)=0;
elseif strcmp(optC,'glcn')
    EX_glcn=find(strcmp(model.rxns,'EX_glcn_e_')); % gluconate
    t_glcn=find(strcmp(model.rxns,'GLCNt2r')); % gluconate
    rxns_Cc=EX_glcn;
    rxns_Cm=t_glcn;
    model.lb(EX_glcn)=-1000;
    model.lb(EX_glc)=0;
elseif strcmp(optC,'arab')
    EX_arab=find(strcmp(model.rxns,'EX_arab_L_e_')); % arabinose
    t_arab=find(ismember(model.rxns,{'ARBabc','ARBt2r'} )); % arabinose
    rxns_Cc=EX_arab;
    rxns_Cm=t_arab;
    model.lb(EX_arab)=-1000;
    model.lb(EX_glc)=0;
elseif strcmp(optC,'rib')
    EX_rib=find(strcmp(model.rxns,'EX_rib_D_e_')); % ribose
    t_rib=find(strcmp(model.rxns,'RIBabc')); % ribose
    rxns_Cc=EX_rib;
    rxns_Cm=t_rib;
    model.lb(EX_rib)=-1000;
    model.lb(EX_glc)=0;
elseif strcmp(optC,'xyl')
    EX_xyl=find(strcmp(model.rxns,'EX_xyl_D_e_')); % xylose
    t_xyl=find(ismember(model.rxns,{'XYLabc','XYLt2'} )); % xylose
    rxns_Cc=EX_xyl;
    rxns_Cm=t_xyl;
    model.lb(EX_xyl)=-1000;
    model.lb(EX_glc)=0;
elseif strcmp(optC,'gal')
    EX_gal=find(strcmp(model.rxns,'EX_gal_e_')); % galactose
    t_gal=find(ismember(model.rxns,{'GALabc','GALt2'} )); % galactose
    rxns_Cc=EX_gal;
    rxns_Cm=t_gal;
    model.lb(EX_gal)=-1000;
    model.lb(EX_glc)=0;


elseif strcmp(optC,'g6p+glcn')
    EX_g6p=find(strcmp(model.rxns,'EX_g6p_e_')); % g6p transport
    t_g6p=find(strcmp(model.rxns,'G6Pt6_2')); % g6p transport
    EX_glcn=find(strcmp(model.rxns,'EX_glcn_e_')); % gluconate
    t_glcn=find(strcmp(model.rxns,'GLCNt2r')); % gluconate
    rxns_Cc=[EX_glcn,EX_g6p];
    rxns_Cm=[t_g6p,t_glcn];
    model.lb(EX_glcn)=-1000;
    model.lb(EX_g6p)=-1000;
    model.lb(EX_glc)=0;
elseif strcmp(optC,'man6p+glcn')
    EX_man6p=find(strcmp(model.rxns,'EX_man6p_e_')); % man6p transport
    t_man6p=find(strcmp(model.rxns,'MAN6Pt6_2')); % man6p transport
    EX_glcn=find(strcmp(model.rxns,'EX_glcn_e_')); % gluconate
    t_glcn=find(strcmp(model.rxns,'GLCNt2r')); % gluconate
    rxns_Cc=[EX_glcn,EX_man6p];
    rxns_Cm=[t_man6p,t_glcn];
    model.lb(EX_glcn)=-1000;
    model.lb(EX_man6p)=-1000;
    model.lb(EX_glc)=0; 
    
% Mid glycolysis
elseif strcmp(optC,'glcr')
    EX_glcr=find(strcmp(model.rxns,'EX_glcr_e_')); % glucorate
    t_glcr=find(strcmp(model.rxns,'GLCRt2r')); % glucorate
    rxns_Cc=EX_glcr;
    rxns_Cm=t_glcr;
    model.lb(EX_glcr)=-1000;
    model.lb(EX_glc)=0;
elseif strcmp(optC,'dha')
    EX_dha=find(strcmp(model.rxns,'EX_dha_e_')); % dihydroxyacetone
    t_dha=find(strcmp(model.rxns,'DHAt')); % dihydroxyacetone
    rxns_Cc=EX_dha;
    rxns_Cm=t_dha;
    model.lb(EX_dha)=-1000;
    model.lb(EX_glc)=0;
elseif strcmp(optC,'glyc')
    EX_glyc=find(strcmp(model.rxns,'EX_glyc_e_')); % glycerol
    t_glyc=find(strcmp(model.rxns,'GLYCt')); % glycerol
    rxns_Cc=EX_glyc;
    rxns_Cm=t_glyc;
    model.lb(EX_glyc)=-1000;
    model.lb(EX_glc)=0;

%TODO: Continue adding more reactions
    
    % TCA C-lim
elseif strcmp(optC,'pyr')
    EX_pyr=find(strcmp(model.rxns,'EX_pyr_e_')); % pyruvate
    rxns_Cc=EX_pyr;
    model.lb(EX_pyr)=-1000;
    model.lb(EX_glc)=0; 
elseif strcmp(optC,'akg')
    EX_akg=find(strcmp(model.rxns,'EX_akg_e_')); % akg
    rxns_Cc=EX_akg;
    model.lb(EX_akg)=-1000;
    model.lb(EX_glc)=0; 
elseif strcmp(optC,'cit')
    EX_cit=find(strcmp(model.rxns,'EX_cit_e_')); % citrate
    rxns_Cc=EX_cit;
    model.lb(EX_cit)=-1000;
    model.lb(EX_glc)=0; 
elseif strcmp(optC,'succ')
    EX_succ=find(strcmp(model.rxns,'EX_succ_e_')); % succinate
    rxns_Cc=EX_succ;
    model.lb(EX_succ)=-1000;
    model.lb(EX_glc)=0;
elseif strcmp(optC,'fum')
    EX_fum=find(strcmp(model.rxns,'EX_fum_e_')); % fumarate
    rxns_Cc=EX_fum;
    model.lb(EX_fum)=-1000;
    model.lb(EX_glc)=0; 
elseif strcmp(optC,'mal')
    EX_mal=find(strcmp(model.rxns,'EX_mal_L_e_') | strcmp(model.rxns,'EX_mal-L_e_')); % malate
    rxns_Cc=EX_mal;
    model.lb(EX_mal)=-1000;
    model.lb(EX_glc)=0; 
    
% other C-lim
elseif strcmp(optC,'lac')
    EX_lac=find(strcmp(model.rxns,'EX_lac_D_e_') | strcmp(model.rxns,'EX_lac_D_e_')); % lactate transport
    rxns_Cc=EX_lac;
    model.lb(EX_lac)=-1000;
    model.lb(EX_glc)=0; 
elseif strcmp(optC,'ac')
    EX_ac=find(strcmp(model.rxns,'EX_ac_e_')); % acetate
    rxns_Cc=EX_ac;
    model.lb(EX_ac)=-1000;
    model.lb(EX_glc)=0;
elseif strcmp(optC,'for')
    EX_for=find(strcmp(model.rxns,'EX_for_e_')); % formate
    rxns_Cc=EX_for;
    model.lb(EX_for)=-1000;
    model.lb(EX_glc)=0;
elseif strcmp(optC,'etoh')
    EX_etoh=find(strcmp(model.rxns,'EX_etoh_e_')); % ethanol
    rxns_Cc=EX_etoh;
    model.lb(EX_etoh)=-1000;
    model.lb(EX_glc)=0;
else
    disp('Error in ''addProteinGroupsToModel'': ''optC'' is not valid. Exiting.');
    return;
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION - Set the PhiE rxns vector  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [phiEcVec, phiErVec, phiEmVec]=set_PhiE_list(model,phiEopt)

nRxns=length(model.rxns);
phiEcVec=[];

% find the constraint vector
if strcmpi(phiEopt, 'All')    
    for n=1:nRxns
        if n~=nGrowth
            phiEcVec=cat(1,phiEcVec,n);
        end
    end
% elseif strcmpi(phiEopt, 'subSys_w_ACM')
%     % This option works with iJR904, iAF1260 and iJO1366
%     subSysList={'Transport Extracellular'; ...
%                 'Transport Inner Membrane'; ...
%                 'Transport Outer Membrane'; ...
%                 'Transport Outer Membrane Porin'; ...
%                 'Putative Transporters'; ...
%                 'Exchange'; ...
%                 'Others'; ...
%                 'Alternate Carbon Metabolism';...
%                 ''};
elseif strcmpi(phiEopt, 'subSys')
    % This option works with iJR904, iAF1260 and iJO1366
    subSysL_q={'Transport Extracellular'; ...
                'Transport Inner Membrane'; ...
                'Transport Outer Membrane'; ...
                'Transport Outer Membrane Porin'; ...
                'Putative Transporters'; ...
                'Exchange'; ...
                'Others'; ...
                'Alternate Carbon Metabolism';...
                ''};
    subSysL_m = {'Transport Extracellular'; ...
                'Transport Inner Membrane'; ...
                'Putative Transporters'};
    tprt_diff = {'CO2t', 'DHAt', 'FORt', 'GLYALDt', 'H2Ot', 'O2t', 'UREAt'};        
    
    M_subSysList = cell([length(model.subSystems) 1]);
    for k1 = 1:length(model.subSystems)
        M_subSysList(k1) =  model.subSystems{k1};
    end
    M_subSysList = M_subSysList(1:k1);
    if contains(model.modelID,'iJR904','IgnoreCase',true)
        Resp_rxnlist = { 'ATPS4r'; 'CYTBD'; 'CYTBO3'; 'FDH2';'FDH3'; 
            'NADH5'; 'NADH6'; 'NADH7'; 'NADH8'; 'NADH9'};
    else
        Resp_rxnlist = { 'ATPS4rpp'; 'CYTBD2pp'; 'CYTBDpp'; 'CYTBO3_4pp';...
            'FDH4pp';'FDH5pp'; 'NADH16pp';'NADH5';'NADH17pp';...
            'NADH18pp';'NADH9'};
    end
    phiErVec = find(ismember(model.rxns, Resp_rxnlist));
    phiEcVec=find(~ismember(M_subSysList,subSysL_q));
    phiEcVec= setdiff(phiEcVec,phiErVec);
    phiEmVec= find(ismember(M_subSysList,subSysL_m));
    difVec = find(ismember(model.rxns, tprt_diff)); 
    phiEmVec= setdiff(phiEmVec,difVec);   % Remove difussion transport
%     fprintf('%i reactions inserted in phiE.\n',length(phiEcVec));
else
    disp('No valid option in set_PhiE_list.');
end
end

function vec3=removeRxnsOverlap(vec1,vec2)
% Subtracts from vec1 all entries  which also appear in vec2
vec3=vec1(~ismember(vec1,vec2));
end