%% MTN_AllInOne_FSA_v4_RegimeStructured_3INPUT.m
% Regime-structured Factorial Sensitivity Analysis (FSA) for chronic loading.
%
% Updated to the 3-input exclusive loading architecture:
%   Hypo, NL, HL
%
% Regimes / domains analyzed:
%   1) Hypo domain            : Hypo scan, NL fixed low, HL fixed low
%   2) Hypo-Normal transition : Hypo x NL factorial, HL fixed low
%   3) Normal domain          : NL scan, Hypo fixed low, HL fixed low
%   4) Normal-Hyper transition: NL x HL factorial, Hypo fixed low
%   5) Hyper domain           : HL scan, Hypo fixed low, NL fixed hyper-baseline low
%
% Notes:
% - Uses the same chronic regime philosophy as the updated baseline/rescue codes
% - Robust clamping logic works even if a clamped value is 0
% - Transition domains produce ANOVA/Sobol-like decompositions
% - Pure domains produce response-amplitude analyses (range and endpoint shift)
% - Exports CSV, Excel, MAT and figure outputs per domain
%
% Dependencies:
%   CreateMatrices_new.m
%   ODESysFun.m
%
% Outputs:
%   ./outputs_FSA_RegimeStructured
%   ./Figures_FSA_RegimeStructured

clc; close all; clear;

%% ------------------------- USER SETTINGS --------------------------------
netfile    = 'MT_PRIMARY4_1.xlsx';
out_xlsx   = 'Mechanotransduction_Results_FSA_RegimeStructured.xlsx';
fig_dir    = fullfile(pwd,'Figures_FSA_RegimeStructured');    ensure_dir(fig_dir);
out_dir    = fullfile(pwd,'outputs_FSA_RegimeStructured');    ensure_dir(out_dir);

% Domain-specific factor levels

% Hypo is biologically restricted to < 0.2
Hypo_levels_domain = [0.10 0.15 0.20 0.25 0.30];

% Hypo -> Normal transition:
% Hypo decreases while NL increases
Hypo_levels_trans  = [0.30 0.25 0.20 0.15 0.10];
NL_levels_trans_HN = [0.10 0.15 0.20 0.25 0.30];

% Normal domain
NL_levels_domain   = [0.10 0.15 0.20 0.25 0.30];

% Normal -> Hyper transition:
% NL decreases while HL increases
NL_levels_trans_NH = [0.30 0.35 0.20 0.15 0.10];
HL_levels_trans    = [0.10 0.15 0.20 0.25 0.30];

% Hyper domain
HL_levels_domain   = [0.10 0.15 0.20 0.25 0.30];

% Optional biological invalid-combination filters for the 2-factor transition domains
HN_Hypo_high_thresh = 0.18;   % Hypo-Normal: exclude relatively high Hypo with high NL
HN_NL_high_thresh   = 0.30;

NH_NL_high_thresh   = 0.30;   % Normal-Hyper: exclude high NL with high HL
NH_HL_high_thresh   = 0.30;

% Simulation controls
Nrep = 100;
tspan = [0 100];
rng(1);

USE_PAR = true;
PAR_WORKERS = [];

RelTol  = 1e-8;
AbsTol  = 1e-10;
MaxStep = 0.5;

% Optional response surfaces for these nodes (2-factor domains only)
SURFACE_NODES = {'SOX9','NF-κB','ROS'};

% Plot controls
TOPN = 50;
YLIM_ACTIVATION = [0 1];

%% ---------------------- GROUP CATEGORIES --------------------------------
group_categories = {
    {'Cytokines, chemokines, proteases & others', ...
        {'TNF','IL6','IL1β','IL8','CCL2','CXCL1','CXCL3','ADAMTS4/5','MMP1','MMP3','MMP13'}};

    {'Growth factors', ...
        {'TGFβ','VEGF','IGF1','BMP2','CCN2','GDF5','FGF2','FGF18','Wnt3a','Wnt5a'}};

    {'Transcription Factors', ...
        {'CREB','HIF-1α','HIF-2α','NF-κB','AP-1','FOXO','SOX9','NFAT','RUNX2', ...
         'YAP/TAZ','MRTF-A','NRF2','HSF1','TonEBP','ELK1','PPARγ','CITED2'}};

    {'Ion channels & related', ...
        {'Ca2+os','Ca2+su','CaMKII','PKC-E','PKC-M','PLCγ-M','PLCγ-E','CaN','IP3','PLA2','AQP1','AQP5'}};

    {'Metabolic & related', ...
        {'LKB1','NAD+','AMPK','mTORC1','mTORC2','SIRT1','PI3K-M','PI3K-E','PIP3-M','PIP3-E', ...
         'PDK1-M','PDK1-E','AKT1-M','AKT1-E','GSK3B','ULK1','PTEN','PLD2','PGE2','COX-2', ...
         'CAT','GPX1','SOD1','SOD2','HO-1','PHD2','VHL','Rheb'}};

    {'ECM anabolism & phenotype markers', ...
        {'COL2A1','COL1A1','COL10A1','ACAN','TIMP3'}};

    {'Mechanical stimuli & its receptors', ...
        {'Hypo','NL','HL','α5β1-FN','α5β1-Fs','αVβ3','αVβ6','SDC4-E','SDC4-M','TRPV4','PIEZO1','NutD','MitD'}};

    {'Cell survival, apoptosis & mitophagy/DNA-damage', ...
        {'Bcl2','BAX','CASP3','CASP9','BNIP3','GADD45','DRP1','MOMP'}};

    {'Oxidative-stress defense & proteostasis', ...
        {'HO-1','GPX1','SOD1','SOD2','CAT','HSP70','HSP27','ROS'}};

    {'MAPK & stress-activated kinases', ...
        {'RAS-M','RAS-E','RAF-M','RAF-E','MEK-M','MEK-E','ERK-M','ERK-E','MKK3/6','MKK4/7','JNK','p38','RSK','TAK1'}};

    {'Rho GTPases, cytoskeletal & Hippo regulators', ...
        {'RhoA-M','RhoA-E','RAC1-M','RAC1-E','CDC42','ROCK-M','ROCK-E','PAK1','PKN1','FAK-M','FAK-E','MST1/2','LATS1/2'}};
};

%% ------------------------- LOAD NETWORK ---------------------------------
[Mact, Minh, NodeNames, NumOfNodes] = CreateMatrices_new(netfile);
NodeNames = string(NodeNames(:));
NodeLabels = ensure_unique_labels(cellstr(NodeNames)); %#ok<NASGU>

idxHypo = find(NodeNames=="Hypo",1);
idxNL   = find(NodeNames=="NL",1);
idxHL   = find(NodeNames=="HL",1);
assert(~isempty(idxHypo) && ~isempty(idxNL) && ~isempty(idxHL), 'Hypo/NL/HL nodes not found.');

clamp_idx = [idxHypo; idxNL; idxHL];

varNames = matlab.lang.makeValidName(cellstr(NodeNames), 'ReplacementStyle','delete');
varNames = matlab.lang.makeUniqueStrings(varNames, {}, namelengthmax); %#ok<NASGU>

assignin('base','Mact', Mact);
assignin('base','Minh', Minh);
assignin('base','NodeNames', NodeNames);
assignin('base','NumOfNodes', NumOfNodes);

%% ------------------------- ODE OPTIONS ----------------------------------
S = (Mact~=0) | (Minh~=0);
S = S | eye(size(S));
ode_opts = odeset('RelTol',RelTol,'AbsTol',AbsTol,'MaxStep',MaxStep, ...
                  'NonNegative',1:NumOfNodes,'JPattern',S);
solver_handle = @ode45;

%% ------------------------- PARPOOL --------------------------------------
BASE_VARS = [];
if USE_PAR
    try
        p = gcp('nocreate');
        if isempty(p)
            if isempty(PAR_WORKERS)
                p = parpool('local');
            else
                p = parpool('local', PAR_WORKERS);
            end
        end
        try
            ofn = which('ODESysFun');
            if ~isempty(ofn), addAttachedFiles(p, {ofn}); end
        catch
        end
        BASE_VARS = parallel.pool.Constant(@() initWorkerBaseVars(Mact,Minh,NodeNames,NumOfNodes));
    catch ME
        warning('Parallel unavailable (%s). Running serially.', ME.message);
        USE_PAR = false;
    end
end

simulate_replicates = @(HypoVal,NLVal,HLVal) run_replicates_clamped3( ...
    HypoVal,NLVal,HLVal,clamp_idx,NumOfNodes,tspan,Nrep,solver_handle,ode_opts,BASE_VARS);

%% ------------------------- DOMAIN CONFIG --------------------------------
Domain = struct('name',{},'type',{},'factor1_name',{},'factor1_levels',{}, ...
                'factor2_name',{},'factor2_levels',{}, ...
                'fixedHypo',{},'fixedNL',{},'fixedHL',{}, ...
                'valid_mask_fun',{});

% 1) Hypo domain
Domain(end+1) = struct( ...
    'name','Hypo_Domain', ...
    'type','1D', ...
    'factor1_name','Hypo', ...
    'factor1_levels',Hypo_levels_domain, ...
    'factor2_name','', ...
    'factor2_levels',[], ...
    'fixedHypo',NaN, ...
    'fixedNL',0.01, ...
    'fixedHL',0.01, ...
    'valid_mask_fun',[]);

% 2) Hypo -> Normal transition
Domain(end+1) = struct( ...
    'name','Hypo_to_Normal', ...
    'type','2D', ...
    'factor1_name','Hypo', ...
    'factor1_levels',Hypo_levels_trans, ...
    'factor2_name','NL', ...
    'factor2_levels',NL_levels_trans_HN, ...
    'fixedHypo',NaN, ...
    'fixedNL',NaN, ...
    'fixedHL',0.01, ...
    'valid_mask_fun',@(A,B) ~((A > HN_Hypo_high_thresh) & (B > HN_NL_high_thresh)));

% 3) Normal domain
Domain(end+1) = struct( ...
    'name','Normal_Domain', ...
    'type','1D', ...
    'factor1_name','NL', ...
    'factor1_levels',NL_levels_domain, ...
    'factor2_name','', ...
    'factor2_levels',[], ...
    'fixedHypo',0.01, ...
    'fixedNL',NaN, ...
    'fixedHL',0.01, ...
    'valid_mask_fun',[]);

% 4) Normal -> Hyper transition
Domain(end+1) = struct( ...
    'name','Normal_to_Hyper', ...
    'type','2D', ...
    'factor1_name','NL', ...
    'factor1_levels',NL_levels_trans_NH, ...
    'factor2_name','HL', ...
    'factor2_levels',HL_levels_trans, ...
    'fixedHypo',0.01, ...
    'fixedNL',NaN, ...
    'fixedHL',NaN, ...
    'valid_mask_fun',@(A,B) ~((A > NH_NL_high_thresh) & (B > NH_HL_high_thresh)));

% 5) Hyper domain
Domain(end+1) = struct( ...
    'name','Hyper_Domain', ...
    'type','1D', ...
    'factor1_name','HL', ...
    'factor1_levels',HL_levels_domain, ...
    'factor2_name','', ...
    'factor2_levels',[], ...
    'fixedHypo',0.01, ...
    'fixedNL',0.10, ...
    'fixedHL',NaN, ...
    'valid_mask_fun',[]);

%% ------------------------- RUN DOMAINS ----------------------------------
SummaryMaster = table( ...
    strings(0,1), ...  % Node
    strings(0,1), ...  % Domain
    strings(0,1), ...  % Type
    strings(0,1), ...  % Factor1
    strings(0,1), ...  % Factor2
    zeros(0,1), ...    % ResponseRange
    zeros(0,1), ...    % EndpointShift
    zeros(0,1), ...    % S_Factor1
    zeros(0,1), ...    % S_Factor2
    zeros(0,1), ...    % S_Interaction
    'VariableNames', {'Node','Domain','Type','Factor1','Factor2', ...
                      'ResponseRange','EndpointShift','S_Factor1','S_Factor2','S_Interaction'});

for d = 1:numel(Domain)
    D = Domain(d);
    fprintf('\n=== Running domain: %s ===\n', D.name);

    dom_fig_dir = fullfile(fig_dir, D.name); ensure_dir(dom_fig_dir);
    dom_out_dir = fullfile(out_dir, D.name); ensure_dir(dom_out_dir);

    if strcmpi(D.type,'1D')
        RES = run_domain_1D(D, simulate_replicates, NodeNames, NumOfNodes, Nrep);
        save(fullfile(dom_out_dir,[D.name '_Results.mat']), 'RES','-v7.3');

        writetable(RES.ComboTable, fullfile(dom_out_dir,[D.name '_ComboMeans.csv']));
        writetable(RES.ResponseTable, fullfile(dom_out_dir,[D.name '_ResponseSummary.csv']));
        try
            writetable(RES.ComboTable, out_xlsx, 'Sheet', safe_sheet_name([D.name '_Combos']));
            writetable(RES.ResponseTable, out_xlsx, 'Sheet', safe_sheet_name([D.name '_Summary']));
        catch ME
            warning('Excel export failed for %s (%s).', D.name, ME.message);
        end

        % Plots for 1D domains
        plot_top_range_bar_1D(RES, NodeNames, TOPN, dom_fig_dir);
        plot_top_heatmap_1D(RES, NodeNames, TOPN, dom_fig_dir);
        plot_group_heatmaps_1D(RES, NodeNames, group_categories, dom_fig_dir);

        % Add to master summary
        Tmaster = table( ...
            NodeNames(:), ...
            repmat(string(D.name),NumOfNodes,1), ...
            repmat("1D",NumOfNodes,1), ...
            repmat(string(D.factor1_name),NumOfNodes,1), ...
            repmat("",NumOfNodes,1), ...
            RES.ResponseRange(:), ...
            RES.EndpointShift(:), ...
            nan(NumOfNodes,1), ...
            nan(NumOfNodes,1), ...
            nan(NumOfNodes,1), ...
            'VariableNames', {'Node','Domain','Type','Factor1','Factor2', ...
                              'ResponseRange','EndpointShift','S_Factor1','S_Factor2','S_Interaction'});
        SummaryMaster = [SummaryMaster; Tmaster]; %#ok<AGROW>

    elseif strcmpi(D.type,'2D')
        RES = run_domain_2D(D, simulate_replicates, NodeNames, NumOfNodes, Nrep);
        save(fullfile(dom_out_dir,[D.name '_Results.mat']), 'RES','-v7.3');

        writetable(RES.ComboTable, fullfile(dom_out_dir,[D.name '_ComboMeans.csv']));
        writetable(RES.SensTable_A, fullfile(dom_out_dir,[D.name '_S_' D.factor1_name '.csv']));
        writetable(RES.SensTable_B, fullfile(dom_out_dir,[D.name '_S_' D.factor2_name '.csv']));
        writetable(RES.SensTable_I, fullfile(dom_out_dir,[D.name '_S_Interaction.csv']));
        try
            writetable(RES.ComboTable, out_xlsx, 'Sheet', safe_sheet_name([D.name '_Combos']));
            writetable(RES.SensTable_A, out_xlsx, 'Sheet', safe_sheet_name([D.name '_S_' D.factor1_name]));
            writetable(RES.SensTable_B, out_xlsx, 'Sheet', safe_sheet_name([D.name '_S_' D.factor2_name]));
            writetable(RES.SensTable_I, out_xlsx, 'Sheet', safe_sheet_name([D.name '_S_Int']));
        catch ME
            warning('Excel export failed for %s (%s).', D.name, ME.message);
        end

        % Plots for 2D domains
        plot_top_heatmap_2D(RES, NodeNames, TOPN, dom_fig_dir);
        plot_group_heatmaps_2D(RES, NodeNames, group_categories, dom_fig_dir);
        plot_surface_nodes_2D(RES, NodeNames, SURFACE_NODES, dom_fig_dir);

        Tmaster = table( ...
            NodeNames(:), ...
            repmat(string(D.name),NumOfNodes,1), ...
            repmat("2D",NumOfNodes,1), ...
            repmat(string(D.factor1_name),NumOfNodes,1), ...
            repmat(string(D.factor2_name),NumOfNodes,1), ...
            nan(NumOfNodes,1), ...
            nan(NumOfNodes,1), ...
            RES.S_A(:), ...
            RES.S_B(:), ...
            RES.S_Int(:), ...
            'VariableNames', {'Node','Domain','Type','Factor1','Factor2', ...
                              'ResponseRange','EndpointShift','S_Factor1','S_Factor2','S_Interaction'});
        SummaryMaster = [SummaryMaster; Tmaster]; %#ok<AGROW>
    else
        error('Unknown domain type for %s.', D.name);
    end
end

writetable(SummaryMaster, fullfile(out_dir,'FSA_MasterSummary.csv'));
try
    writetable(SummaryMaster, out_xlsx, 'Sheet', safe_sheet_name('FSA_MasterSummary'));
catch
end

disp('Regime-structured FSA complete: all results and figures saved.');

%% ========================= HELPER FUNCTIONS =============================

function token = initWorkerBaseVars(Mact,Minh,NodeNames,NumOfNodes)
assignin('base','Mact',Mact);
assignin('base','Minh',Minh);
assignin('base','NodeNames',NodeNames);
assignin('base','NumOfNodes',NumOfNodes);
token = true;
end

function [SS_all, SS_mu, SS_sd] = run_replicates_clamped3(HypoVal,NLVal,HLVal,clamp_idx,N,tspan,Nrep,solver,opts,BASE)
SS_all = nan(Nrep,N);
clamp_vals = [HypoVal; NLVal; HLVal];
clamp_mask = false(N,1);
clamp_mask(clamp_idx) = true;

if ~isempty(gcp('nocreate'))
    parfor r = 1:Nrep
        if ~isempty(BASE); BASE.Value; end
        Xi = rand(N,1);
        Xi(clamp_idx) = clamp_vals;
        [~,X] = solver(@(t,x) ODESysFun_Clamped_safe3(t,x,clamp_idx,clamp_vals,clamp_mask), tspan, Xi, opts);
        SS_all(r,:) = X(end,:);
    end
else
    for r = 1:Nrep
        Xi = rand(N,1);
        Xi(clamp_idx) = clamp_vals;
        [~,X] = solver(@(t,x) ODESysFun_Clamped_safe3(t,x,clamp_idx,clamp_vals,clamp_mask), tspan, Xi, opts);
        SS_all(r,:) = X(end,:);
    end
end
SS_mu = mean(SS_all,1,'omitnan');
SS_sd = std(SS_all,0,1,'omitnan');
end

function dxdt = ODESysFun_Clamped_safe3(~,x,clamp_idx,clamp_vals,clamp_mask)
x(clamp_idx) = clamp_vals;
dxdt = ODESysFun(0,x);
dxdt(clamp_mask) = 0;
end

function RES = run_domain_1D(D, simulate_replicates, NodeNames, NumOfNodes, Nrep)
L = D.factor1_levels(:);
K = numel(L);

ComboMeans = nan(K, NumOfNodes);
ComboSDs   = nan(K, NumOfNodes);

fprintf('1D domain "%s": %d levels x %d reps = %d ODEs\n', D.name, K, Nrep, K*Nrep);

for k = 1:K
    [HypoVal, NLVal, HLVal] = build_input_1D(D, L(k));
    [~,SS_mu,SS_sd] = simulate_replicates(HypoVal,NLVal,HLVal);
    ComboMeans(k,:) = SS_mu;
    ComboSDs(k,:)   = SS_sd;
end

ResponseRange = max(ComboMeans,[],1) - min(ComboMeans,[],1);
EndpointShift = ComboMeans(end,:) - ComboMeans(1,:);
ResponseStd   = std(ComboMeans,0,1,'omitnan');

ComboTable = table();
ComboTable.Level = L;
ComboTable.Hypo  = nan(K,1);
ComboTable.NL    = nan(K,1);
ComboTable.HL    = nan(K,1);
for k = 1:K
    [HypoVal, NLVal, HLVal] = build_input_1D(D, L(k));
    ComboTable.Hypo(k) = HypoVal;
    ComboTable.NL(k)   = NLVal;
    ComboTable.HL(k)   = HLVal;
end
for n = 1:NumOfNodes
    ComboTable.(matlab.lang.makeValidName(char(NodeNames(n)))) = ComboMeans(:,n);
end

ResponseTable = table(NodeNames(:), ResponseRange(:), EndpointShift(:), ResponseStd(:), ...
    'VariableNames', {'Node','ResponseRange','EndpointShift','AcrossLevelStd'});

RES = struct();
RES.DomainName    = D.name;
RES.Type          = D.type;
RES.FactorName    = D.factor1_name;
RES.Levels        = L;
RES.ComboMeans    = ComboMeans;
RES.ComboSDs      = ComboSDs;
RES.ResponseRange = ResponseRange(:);
RES.EndpointShift = EndpointShift(:);
RES.ResponseStd   = ResponseStd(:);
RES.ComboTable    = ComboTable;
RES.ResponseTable = ResponseTable;
end

function [HypoVal, NLVal, HLVal] = build_input_1D(D, level)
HypoVal = D.fixedHypo;
NLVal   = D.fixedNL;
HLVal   = D.fixedHL;

switch D.factor1_name
    case 'Hypo'
        HypoVal = level;
    case 'NL'
        NLVal   = level;
    case 'HL'
        HLVal   = level;
    otherwise
        error('Unknown factor1_name "%s" in 1D domain.', D.factor1_name);
end
end

function RES = run_domain_2D(D, simulate_replicates, NodeNames, NumOfNodes, Nrep)
A_levels = D.factor1_levels(:);
B_levels = D.factor2_levels(:);

[Agrid, Bgrid] = ndgrid(A_levels, B_levels);
Avec = Agrid(:);
Bvec = Bgrid(:);

if isempty(D.valid_mask_fun)
    valid = true(size(Avec));
else
    valid = D.valid_mask_fun(Avec, Bvec);
end

Avec = Avec(valid);
Bvec = Bvec(valid);
K = numel(Avec);

fprintf('2D domain "%s": %d valid combos x %d reps = %d ODEs\n', D.name, K, Nrep, K*Nrep);

ComboMeans = nan(K, NumOfNodes);
ComboSDs   = nan(K, NumOfNodes);

for k = 1:K
    [HypoVal, NLVal, HLVal] = build_input_2D(D, Avec(k), Bvec(k));
    [~,SS_mu,SS_sd] = simulate_replicates(HypoVal,NLVal,HLVal);
    ComboMeans(k,:) = SS_mu;
    ComboSDs(k,:)   = SS_sd;
end

[S_A, S_B, S_Int] = sensitivity_decompose_2D(Avec, Bvec, A_levels, B_levels, ComboMeans);

ComboTable = table();
ComboTable.(char(D.factor1_name)) = Avec;
ComboTable.(char(D.factor2_name)) = Bvec;
ComboTable.Hypo = nan(K,1);
ComboTable.NL   = nan(K,1);
ComboTable.HL   = nan(K,1);
for k = 1:K
    [HypoVal, NLVal, HLVal] = build_input_2D(D, Avec(k), Bvec(k));
    ComboTable.Hypo(k) = HypoVal;
    ComboTable.NL(k)   = NLVal;
    ComboTable.HL(k)   = HLVal;
end
for n = 1:NumOfNodes
    ComboTable.(matlab.lang.makeValidName(char(NodeNames(n)))) = ComboMeans(:,n);
end

SensTable_A = table(NodeNames(:), S_A(:), 'VariableNames', {'Node',['S_' char(D.factor1_name)]});
SensTable_B = table(NodeNames(:), S_B(:), 'VariableNames', {'Node',['S_' char(D.factor2_name)]});
SensTable_I = table(NodeNames(:), S_Int(:), 'VariableNames', {'Node','S_Interaction'});

RES = struct();
RES.DomainName = D.name;
RES.Type       = D.type;
RES.FactorA    = D.factor1_name;
RES.FactorB    = D.factor2_name;
RES.A_levels   = A_levels;
RES.B_levels   = B_levels;
RES.Avec       = Avec;
RES.Bvec       = Bvec;
RES.ComboMeans = ComboMeans;
RES.ComboSDs   = ComboSDs;
RES.S_A        = S_A(:);
RES.S_B        = S_B(:);
RES.S_Int      = S_Int(:);
RES.ComboTable = ComboTable;
RES.SensTable_A= SensTable_A;
RES.SensTable_B= SensTable_B;
RES.SensTable_I= SensTable_I;
end

function [HypoVal, NLVal, HLVal] = build_input_2D(D, A, B)
HypoVal = D.fixedHypo;
NLVal   = D.fixedNL;
HLVal   = D.fixedHL;

switch D.factor1_name
    case 'Hypo'
        HypoVal = A;
    case 'NL'
        NLVal   = A;
    case 'HL'
        HLVal   = A;
    otherwise
        error('Unknown factor1_name "%s".', D.factor1_name);
end

switch D.factor2_name
    case 'Hypo'
        HypoVal = B;
    case 'NL'
        NLVal   = B;
    case 'HL'
        HLVal   = B;
    otherwise
        error('Unknown factor2_name "%s".', D.factor2_name);
end
end

function [S_A, S_B, S_Int] = sensitivity_decompose_2D(Avec, Bvec, A_levels, B_levels, ComboMeans)
[~, A_idx] = ismembertol(Avec, A_levels, 1e-12);
[~, B_idx] = ismembertol(Bvec, B_levels, 1e-12);

K = numel(Avec);
w_A = accumarray(A_idx, 1, [numel(A_levels) 1]) / K;
w_B = accumarray(B_idx, 1, [numel(B_levels) 1]) / K;

NumOfNodes = size(ComboMeans,2);
S_A = zeros(NumOfNodes,1);
S_B = zeros(NumOfNodes,1);
S_Int = zeros(NumOfNodes,1);

A_group = cell(numel(A_levels),1);
B_group = cell(numel(B_levels),1);
for i = 1:numel(A_levels), A_group{i} = find(A_idx==i); end
for j = 1:numel(B_levels), B_group{j} = find(B_idx==j); end

for n = 1:NumOfNodes
    y = ComboMeans(:,n);
    mu = mean(y,'omitnan');
    Vt = var(y,1,'omitnan');
    if Vt <= eps, continue; end

    fA = nan(numel(A_levels),1);
    for i = 1:numel(A_levels)
        idx = A_group{i};
        if ~isempty(idx), fA(i) = mean(y(idx),'omitnan'); end
    end

    fB = nan(numel(B_levels),1);
    for j = 1:numel(B_levels)
        idx = B_group{j};
        if ~isempty(idx), fB(j) = mean(y(idx),'omitnan'); end
    end

    VA = nansum(w_A .* (fA - mu).^2);
    VB = nansum(w_B .* (fB - mu).^2);
    VI = max(0, Vt - VA - VB);

    S_A(n)   = VA / Vt;
    S_B(n)   = VB / Vt;
    S_Int(n) = VI / Vt;
end

sumS = S_A + S_B + S_Int;
mask = sumS > 0;
S_A(mask)   = S_A(mask)   ./ sumS(mask);
S_B(mask)   = S_B(mask)   ./ sumS(mask);
S_Int(mask) = S_Int(mask) ./ sumS(mask);
end

function plot_top_range_bar_1D(RES, NodeNames, TOPN, dom_fig_dir)
[~,ord] = sort(RES.ResponseRange,'descend');
topN = min(TOPN, numel(ord));
idx = ord(1:topN);

f = figure('Color','w','Position',[100 100 1400 800]);
ax = axes('Parent',f);
bar(ax, RES.ResponseRange(idx), 'FaceColor',[0.2 0.5 0.85], 'EdgeColor','none');
grid(ax,'on'); box(ax,'on');

set(ax, ...
    'XTick',1:topN, ...
    'XTickLabel',cellstr(NodeNames(idx)), ...
    'XTickLabelRotation',45, ...
    'FontWeight','bold', ...
    'FontSize',16, ...
    'TickLabelInterpreter','none', ...
    'LineWidth',1.2);

ylabel(ax,'Response range across levels','FontWeight','bold','FontSize',18);
title(ax,sprintf('%s: Top-%d nodes by response range', pretty_label(RES.DomainName), topN), ...
    'FontWeight','bold','FontSize',20,'Interpreter','none');

exportgraphics(f, fullfile(dom_fig_dir,[RES.DomainName '_TopRangeBar.png']), 'Resolution',600);
close(f);
end

function plot_top_heatmap_1D(RES, NodeNames, TOPN, dom_fig_dir)
[~,ord] = sort(RES.ResponseRange,'descend');
topN = min(TOPN, numel(ord));
idx = ord(1:topN);

Z = RES.ComboMeans(:,idx)'; % nodes x levels
xlab = cellstr(compose('%.2f', RES.Levels(:)));

titleStr = sprintf('%s: Top-%d node responses', pretty_label(RES.DomainName), topN);
xlabStr  = char(RES.FactorName);

draw_custom_heatmap( ...
    Z, ...
    xlab, ...
    cellstr(NodeNames(idx)), ...
    titleStr, ...
    xlabStr, ...
    'Top responsive nodes', ...
    fullfile(dom_fig_dir,[RES.DomainName '_TopHeatmap.png']), ...
    'FigPos',[100 100 1000 1200], ...
    'TitleSize',22, ...
    'LabelSize',20, ...
    'TickSize',18, ...
    'CellTextSize',18);
end

function plot_group_heatmaps_1D(RES, NodeNames, group_categories, dom_fig_dir)
for g = 1:length(group_categories)
    group_name = group_categories{g}{1};
    node_list  = group_categories{g}{2};
    idx = map_custom_names(node_list, cellstr(NodeNames));
    idx = idx(idx>0);
    if isempty(idx), continue; end

    Z = RES.ComboMeans(:,idx)'; % nodes x levels
    xlab = cellstr(compose('%.2f', RES.Levels(:)));

    titleStr = sprintf('%s: %s', pretty_label(RES.DomainName), group_name);
    xlabStr  = char(RES.FactorName);

    draw_custom_heatmap( ...
        Z, ...
        xlab, ...
        cellstr(NodeNames(idx)), ...
        titleStr, ...
        xlabStr, ...
        'Nodes', ...
        fullfile(dom_fig_dir, sprintf('%s_Group_%s.png', ...
            RES.DomainName, regexprep(group_name,'[^A-Za-z0-9_]','_'))), ...
        'FigPos',[150 100 1000 750], ...
        'TitleSize',20, ...
        'LabelSize',18, ...
        'TickSize',16, ...
        'CellTextSize',16);
end
end

function plot_top_heatmap_2D(RES, NodeNames, TOPN, dom_fig_dir)
maxS = max([RES.S_A RES.S_B RES.S_Int],[],2);
[~,ord] = sort(maxS,'descend');
topN = min(TOPN, numel(ord));
idx = ord(1:topN);

SensMat = [RES.S_A(idx) RES.S_B(idx) RES.S_Int(idx)];
xlabs = {tex_s_label(RES.FactorA), tex_s_label(RES.FactorB), 'S_{INT}'};

titleStr = sprintf('%s: Top-%d by max sensitivity', pretty_label(RES.DomainName), topN);

draw_custom_heatmap( ...
    SensMat, ...
    xlabs, ...
    cellstr(NodeNames(idx)), ...
    titleStr, ...
    'Sensitivity component', ...
    sprintf('Top-%d nodes', topN), ...
    fullfile(dom_fig_dir,[RES.DomainName '_TopSensitivityHeatmap.png']), ...
    'FigPos',[100 100 1000 1200], ...
    'TitleSize',22, ...
    'LabelSize',20, ...
    'TickSize',18, ...
    'CellTextSize',18);
end

function plot_group_heatmaps_2D(RES, NodeNames, group_categories, dom_fig_dir)
for g = 1:length(group_categories)
    group_name = group_categories{g}{1};
    node_list  = group_categories{g}{2};
    idx = map_custom_names(node_list, cellstr(NodeNames));
    idx = idx(idx>0);
    if isempty(idx), continue; end

    SensMat = [RES.S_A(idx) RES.S_B(idx) RES.S_Int(idx)];
    xlabs = {tex_s_label(RES.FactorA), tex_s_label(RES.FactorB), 'S_{INT}'};

    titleStr = sprintf('%s: %s', pretty_label(RES.DomainName), group_name);

    draw_custom_heatmap( ...
        SensMat, ...
        xlabs, ...
        cellstr(NodeNames(idx)), ...
        titleStr, ...
        'Sensitivity component', ...
        'Nodes', ...
        fullfile(dom_fig_dir, sprintf('%s_Group_%s.png', ...
            RES.DomainName, regexprep(group_name,'[^A-Za-z0-9_]','_'))), ...
        'FigPos',[150 100 1000 750], ...
        'TitleSize',20, ...
        'LabelSize',18, ...
        'TickSize',16, ...
        'CellTextSize',16);
end
end

function draw_custom_heatmap(Z, xlabels, ylabels, titleStr, xlabStr, ylabStr, outFile, varargin)
    p = inputParser;
    addParameter(p,'FigPos',[100 100 900 700]);
    addParameter(p,'FontSize',16);
    addParameter(p,'TitleSize',20);
    addParameter(p,'LabelSize',18);
    addParameter(p,'TickSize',16);
    addParameter(p,'CellTextSize',16);
    parse(p,varargin{:});
    P = p.Results;

    f = figure('Color','w','Position',P.FigPos);
    ax = axes('Parent',f);

    imagesc(ax, Z);
    colormap(ax, parula);
    cb = colorbar(ax);
    cb.FontSize = P.TickSize;
    cb.FontWeight = 'bold';

    axis(ax,'tight');

    set(ax, ...
        'XTick',1:numel(xlabels), ...
        'XTickLabel',xlabels, ...
        'YTick',1:numel(ylabels), ...
        'YTickLabel',ylabels, ...
        'TickLabelInterpreter','tex', ...
        'FontWeight','bold', ...
        'FontSize',P.TickSize, ...
        'LineWidth',1.2);

    xtickangle(ax,0);

    xlabel(ax, xlabStr, 'Interpreter','tex', 'FontWeight','bold', 'FontSize',P.LabelSize);
    ylabel(ax, ylabStr, 'Interpreter','tex', 'FontWeight','bold', 'FontSize',P.LabelSize);
    title(ax, titleStr, 'Interpreter','none', 'FontWeight','bold', 'FontSize',P.TitleSize);

    hold(ax,'on');
    [nr,nc] = size(Z);

    for i = 0.5:1:(nr+0.5)
        plot(ax,[0.5 nc+0.5],[i i],'k-','LineWidth',0.5);
    end
    for j = 0.5:1:(nc+0.5)
        plot(ax,[j j],[0.5 nr+0.5],'k-','LineWidth',0.5);
    end

    clim = caxis(ax);
    midv = mean(clim);

    for i = 1:nr
        for j = 1:nc
            val = Z(i,j);
            if isnan(val)
                txt = 'NaN';
            else
                txt = sprintf('%.4g', val);
            end

            if val > midv
                txtColor = 'k';
            else
                txtColor = 'w';
            end

            text(ax, j, i, txt, ...
                'HorizontalAlignment','center', ...
                'VerticalAlignment','middle', ...
                'FontWeight','bold', ...
                'FontSize',P.CellTextSize, ...
                'Color',txtColor, ...
                'Interpreter','none');
        end
    end

    set(ax,'YDir','normal');

    exportgraphics(f, outFile, 'Resolution', 600);
    close(f);
end

function s = tex_s_label(factorName)
    factorName = char(string(factorName));
    switch factorName
        case 'NL'
            s = 'S_{NL}';
        case 'HL'
            s = 'S_{HL}';
        case 'Hypo'
            s = 'S_{Hypo}';
        otherwise
            s = ['S_{' factorName '}'];
    end
end

function s = pretty_label(s)
    s = char(string(s));
    s = strrep(s, '_', ' ');
end
function plot_surface_nodes_2D(RES, NodeNames, SURFACE_NODES, dom_fig_dir)
for s = 1:numel(SURFACE_NODES)
    node_name = SURFACE_NODES{s};
    node_idx = find(strcmp(NodeNames, node_name), 1);
    if isempty(node_idx)
        warning('Surface node "%s" not found. Skipping in %s.', node_name, RES.DomainName);
        continue;
    end

    Z = nan(numel(RES.A_levels), numel(RES.B_levels));
    for k = 1:numel(RES.Avec)
        ii = find(abs(RES.A_levels - RES.Avec(k)) < 1e-12, 1);
        jj = find(abs(RES.B_levels - RES.Bvec(k)) < 1e-12, 1);
        if ~isempty(ii) && ~isempty(jj)
            Z(ii,jj) = RES.ComboMeans(k,node_idx);
        end
    end

    [X,Y] = meshgrid(RES.B_levels, RES.A_levels);
    figure('Color','w','Position',[150 120 700 520]);
    surf(X,Y,Z,'EdgeColor','none');
    view(135,30);
    xlabel(char(RES.FactorB));
    ylabel(char(RES.FactorA));
    zlabel(sprintf('%s steady-state', node_name), 'Interpreter','none');
    title(sprintf('%s: %s response surface', pretty_label(RES.DomainName), node_name), 'Interpreter','none');
    colorbar; grid on;
    exportgraphics(gcf, fullfile(dom_fig_dir, sprintf('%s_Surface_%s.png', ...
        RES.DomainName, regexprep(node_name,'[^A-Za-z0-9_]','_'))), 'Resolution',300);
    close(gcf);
end
end

function idx = map_custom_names(custom_names, NodeNames)
if isempty(custom_names), idx = []; return; end
if isstring(custom_names), custom_names = cellstr(custom_names); end
if isstring(NodeNames), NodeNames = cellstr(NodeNames); end

custom_names = custom_names(~cellfun(@isempty, custom_names));
if isempty(custom_names), idx = []; return; end

NN_norm = normalize_for_match(NodeNames);
idx = [];
for k = 1:numel(custom_names)
    c = custom_names{k};
    c_fix = fix_common_synonyms(c);
    c_norm = normalize_for_match(c_fix);
    hit = find(strcmpi(c_norm, NN_norm), 1);
    if isempty(hit)
        hit = find(strcmpi(c_fix, NodeNames), 1);
    end
    if ~isempty(hit)
        idx(end+1) = hit; %#ok<AGROW>
    end
end
[~, ia] = unique(idx, 'stable');
idx = idx(sort(ia));
end

function out = normalize_for_match(strs)
if ischar(strs), strs = {strs}; end
out = cell(size(strs));
for i = 1:numel(strs)
    s = lower(strs{i});
    s = strrep(s, 'α', 'a');
    s = strrep(s, 'β', 'b');
    s = strrep(s, 'γ', 'g');
    s = strrep(s, 'κ', 'k');
    s = regexprep(s, '\s+', '');
    s = strrep(s, '–', '-');
    s = regexprep(s, '[^a-z0-9/\+\-_]', '');
    out{i} = s;
end
end

function s = fix_common_synonyms(s)
switch s
    case 'TPV4', s = 'TRPV4';
    case 'GSK3B', s = 'GSK3β';
    case {'PPAR𝛾','PPAR-gamma'}, s = 'PPARγ';
    case 'α5β1_FN', s = 'α5β1-FN';
    case 'α5β1_Fs', s = 'α5β1-Fs';
    case {'IL-1β','IL1beta','IL-1b','IL1b'}, s = 'IL1β';
end
end

function ensure_dir(d)
if exist(d,'dir')~=7
    mkdir(d);
end
end

function labels = ensure_unique_labels(base)
labels = base(:);
[~, ~, gid] = unique(labels, 'stable');
counts = accumarray(gid, 1);
if any(counts > 1)
    seen = zeros(size(counts));
    for i = 1:numel(labels)
        g = gid(i);
        seen(g) = seen(g) + 1;
        if counts(g) > 1
            labels{i} = sprintf('%s (%d)', labels{i}, seen(g));
        end
    end
end
end

function sh = safe_sheet_name(s)
s = regexprep(s,'[^A-Za-z0-9_]','_');
if numel(s) > 31
    s = s(1:31);
end
sh = s;
end

