function model=addProteinGroupsToModel(model,optC)
% model=addProteinGroupsToModel(model,optC)
%
% This function adds the 'protGroup' field to a COBRA model.
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
    protGroup(1).name='phiC';
    protGroup(1).phi0=0;
    [model,rxns_C] = setPhiC(model,optC);
    protGroup(1).rxns=rxns_C;
    
    % Defines all protGroup entries
    protGroup(2).name='phiE';
    protGroup(2).phi0=0.45;
    protGroup(2).rxns=set_PhiE_list(model,phiEopt);
    protGroup(2).rxns=removeRxnsOverlap(protGroup(2).rxns,protGroup(1).rxns);
    
    protGroup(3).name='phiR';
    protGroup(3).phi0=0.066;
    protGroup(3).rxns=Biomass;
    protGroup(3).rxns=removeRxnsOverlap(protGroup(3).rxns,protGroup(1).rxns);
    protGroup(3).rxns=removeRxnsOverlap(protGroup(3).rxns,protGroup(2).rxns);
    
    protGroup(4).name='phiQ and others'; % weights=0
    protGroup(4).phi0=0;
    protGroup(4).rxns=(1:length(model.rxns))';
    protGroup(4).rxns=removeRxnsOverlap(protGroup(4).rxns,protGroup(1).rxns);
    protGroup(4).rxns=removeRxnsOverlap(protGroup(4).rxns,protGroup(2).rxns);
    protGroup(4).rxns=removeRxnsOverlap(protGroup(4).rxns,protGroup(3).rxns);
    
    model.protGroup=protGroup;
else
    % only updates the first prot group
    [model,rxns_C] = setPhiC(model,optC);
    model.protGroup(1).rxns=rxns_C;
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION - Set the PhiC rxns vector and bounds %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [model,rxns_C] = setPhiC(model,optC)
% This function also updates the bound for the exchange reactions

EX_glc=find(strcmp(model.rxns,'EX_glc_e_')); % glc transport

% Upper glycolysis
if strcmp(optC,'default') || strcmp(optC,'glc')
    rxns_C=EX_glc;
elseif strcmp(optC,'g6p')
    EX_g6p=find(strcmp(model.rxns,'EX_g6p_e_')); % g6p transport
    rxns_C=EX_g6p;
    model.lb(EX_g6p)=-1000;
    model.lb(EX_glc)=0; 
elseif strcmp(optC,'man6p')
    EX_man6p=find(strcmp(model.rxns,'EX_man6p_e_')); % g6p transport
    rxns_C=EX_man6p;
    model.lb(EX_man6p)=-1000;
    model.lb(EX_glc)=0; 
elseif strcmp(optC,'mnl')
    EX_mnl=find(strcmp(model.rxns,'EX_mnl_e_')); % mannitol transport
    rxns_C=EX_mnl;
    model.lb(EX_mnl)=-1000;
    model.lb(EX_glc)=0;
elseif strcmp(optC,'lcts')
    EX_lcts=find(strcmp(model.rxns,'EX_lcts_e_') ); % lactose transport
    rxns_C=EX_lcts;
    model.lb(EX_lcts)=-1000;
    model.lb(EX_glc)=0;
elseif strcmp(optC,'sucr')
    EX_sucr=find(strcmp(model.rxns,'EX_sucr_e_') ); % sucrose transport
    rxns_C=EX_sucr;
    model.lb(EX_sucr)=-1000;
    model.lb(EX_glc)=0;
elseif strcmp(optC,'malt')
    EX_malt=find(strcmp(model.rxns,'EX_malt_e_') ); % maltose transport
    rxns_C=EX_malt;
    model.lb(EX_malt)=-1000;
    model.lb(EX_glc)=0;
elseif strcmp(optC,'fru')
    EX_fru=find(strcmp(model.rxns,'EX_fru_e_')); % fructose transport
    rxns_C=EX_fru;
    model.lb(EX_fru)=-1000;
    model.lb(EX_glc)=0;
elseif strcmp(optC,'sbt')
    EX_sbt=find(strcmp(model.rxns,'EX_sbt_D_e_')); % sorbitol
    rxns_C=EX_sbt;
    model.lb(EX_sbt)=-1000;
    model.lb(EX_glc)=0;
elseif strcmp(optC,'man')
    EX_man=find(strcmp(model.rxns,'EX_man_e_')); % mannose
    rxns_C=EX_man;
    model.lb(EX_man)=-1000;
    model.lb(EX_glc)=0;
elseif strcmp(optC,'glcn')
    EX_glcn=find(strcmp(model.rxns,'EX_glcn_e_')); % gluconate
    rxns_C=EX_glcn;
    model.lb(EX_glcn)=-1000;
    model.lb(EX_glc)=0;
elseif strcmp(optC,'arab')
    EX_arab=find(strcmp(model.rxns,'EX_arab_L_e_')); % arabinose
    rxns_C=EX_arab;
    model.lb(EX_arab)=-1000;
    model.lb(EX_glc)=0;
elseif strcmp(optC,'rib')
    EX_rib=find(strcmp(model.rxns,'EX_rib_D_e_')); % ribulose
    rxns_C=EX_rib;
    model.lb(EX_rib)=-1000;
    model.lb(EX_glc)=0;
elseif strcmp(optC,'xyl')
    EX_xyl=find(strcmp(model.rxns,'EX_xyl_D_e_')); % xylose
    rxns_C=EX_xyl;
    model.lb(EX_xyl)=-1000;
    model.lb(EX_glc)=0;
elseif strcmp(optC,'gal')
    EX_gal=find(strcmp(model.rxns,'EX_gal_e_')); % galactose
    rxns_C=EX_gal;
    model.lb(EX_gal)=-1000;
    model.lb(EX_glc)=0;

elseif strcmp(optC,'g6p+glcn')
    EX_g6p=find(strcmp(model.rxns,'EX_g6p_e_')); % g6p transport
    EX_glcn=find(strcmp(model.rxns,'EX_glcn_e_')); % gluconate
    rxns_C=[EX_glcn,EX_g6p];
    model.lb(EX_glcn)=-1000;
    model.lb(EX_g6p)=-1000;
    model.lb(EX_glc)=0;
elseif strcmp(optC,'man6p+glcn')
    EX_man6p=find(strcmp(model.rxns,'EX_man6p_e_')); % man6p transport
    EX_glcn=find(strcmp(model.rxns,'EX_glcn_e_')); % gluconate
    rxns_C=[EX_glcn,EX_man6p];
    model.lb(EX_glcn)=-1000;
    model.lb(EX_man6p)=-1000;
    model.lb(EX_glc)=0; 
    
% Mid glycolysis
elseif strcmp(optC,'glcr')
    EX_glcr=find(strcmp(model.rxns,'EX_glcr_e_')); % glucorate
    rxns_C=EX_glcr;
    model.lb(EX_glcr)=-1000;
    model.lb(EX_glc)=0;
elseif strcmp(optC,'dha')
    EX_dha=find(strcmp(model.rxns,'EX_dha_e_')); % dihydroxyacetone
    rxns_C=EX_dha;
    model.lb(EX_dha)=-1000;
    model.lb(EX_glc)=0;
elseif strcmp(optC,'glyc')
    EX_glyc=find(strcmp(model.rxns,'EX_glyc_e_')); % glycerol
    rxns_C=EX_glyc;
    model.lb(EX_glyc)=-1000;
    model.lb(EX_glc)=0;

    
    % TCA C-lim
elseif strcmp(optC,'pyr')
    EX_pyr=find(strcmp(model.rxns,'EX_pyr_e_')); % pyruvate
    rxns_C=EX_pyr;
    model.lb(EX_pyr)=-1000;
    model.lb(EX_glc)=0; 
elseif strcmp(optC,'akg')
    EX_akg=find(strcmp(model.rxns,'EX_akg_e_')); % akg
    rxns_C=EX_akg;
    model.lb(EX_akg)=-1000;
    model.lb(EX_glc)=0; 
elseif strcmp(optC,'cit')
    EX_cit=find(strcmp(model.rxns,'EX_cit_e_')); % citrate
    rxns_C=EX_cit;
    model.lb(EX_cit)=-1000;
    model.lb(EX_glc)=0; 
elseif strcmp(optC,'succ')
    EX_succ=find(strcmp(model.rxns,'EX_succ_e_')); % succinate
    rxns_C=EX_succ;
    model.lb(EX_succ)=-1000;
    model.lb(EX_glc)=0;
elseif strcmp(optC,'fum')
    EX_fum=find(strcmp(model.rxns,'EX_fum_e_')); % fumarate
    rxns_C=EX_fum;
    model.lb(EX_fum)=-1000;
    model.lb(EX_glc)=0; 
elseif strcmp(optC,'mal')
    EX_mal=find(strcmp(model.rxns,'EX_mal_L_e_') | strcmp(model.rxns,'EX_mal-L_e_')); % malate
    rxns_C=EX_mal;
    model.lb(EX_mal)=-1000;
    model.lb(EX_glc)=0; 
    
% other C-lim
elseif strcmp(optC,'glcn')
    EX_glcn=find(strcmp(model.rxns,'EX_glcn_e_')); % lactate transport
    rxns_C=EX_glcn;
    model.lb(EX_glcn)=-1000;
    model.lb(EX_glc)=0; 
elseif strcmp(optC,'lac')
    EX_lac=find(strcmp(model.rxns,'EX_lac_D_e_') | strcmp(model.rxns,'EX_lac_D_e_')); % lactate transport
    rxns_C=EX_lac;
    model.lb(EX_lac)=-1000;
    model.lb(EX_glc)=0; 
elseif strcmp(optC,'ac')
    EX_ac=find(strcmp(model.rxns,'EX_ac_e_')); % acetate
    rxns_C=EX_ac;
    model.lb(EX_ac)=-1000;
    model.lb(EX_glc)=0;
elseif strcmp(optC,'for')
    EX_for=find(strcmp(model.rxns,'EX_for_e_')); % formate
    rxns_C=EX_for;
    model.lb(EX_for)=-1000;
    model.lb(EX_glc)=0;
elseif strcmp(optC,'etoh')
    EX_etoh=find(strcmp(model.rxns,'EX_etoh_e_')); % ethanol
    rxns_C=EX_etoh;
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
function phiEVec=set_PhiE_list(model,phiEopt)

nRxns=length(model.rxns);

phiEVec=[];
% find the constraint vector
if strcmpi(phiEopt, 'All')    
    for n=1:nRxns
        if n~=nGrowth
            phiEVec=cat(1,phiEVec,n);
        end
    end
elseif strcmpi(phiEopt, 'subSys_w_ACM')
    % This option works with iJR904, iAF1260 and iJO1366
    subSysList={'Transport Extracellular'; ...
                'Transport Inner Membrane'; ...
                'Transport Outer Membrane'; ...
                'Transport Outer Membrane Porin'; ...
                'Putative Transporters'; ...
                'Exchange'; ...
                'Others'; ...
                'Alternate Carbon Metabolism';...
                ''};
elseif strcmpi(phiEopt, 'subSys')
    % This option works with iJR904, iAF1260 and iJO1366
    subSysList={'Transport Extracellular'; ...
                'Transport Inner Membrane'; ...
                'Transport Outer Membrane'; ...
                'Transport Outer Membrane Porin'; ...
                'Putative Transporters'; ...
                'Exchange'; ...
                'Others'; ...
                'Alternate Carbon Metabolism';...
                ''};
    M_subSysList ={};
    for k1 = 1:length(model.subSystems)
        M_subSysList = [M_subSysList; model.subSystems{k1}];
    end
    phiEVec=find(~ismember(M_subSysList,subSysList));
    %fprintf('%i reactions inserted in phiE.\n',length(phiEVec));
else
    disp('No valid option in set_PhiE_list.');
end
end

function vec3=removeRxnsOverlap(vec1,vec2)
% Subtracts from vec1 all entries  which also appear in vec2
vec3=vec1(~ismember(vec1,vec2));
end