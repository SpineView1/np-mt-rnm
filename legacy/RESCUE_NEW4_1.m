%% NP_MorPerturb_From_Hyper_ALLINONE_PAIRED_PARALLEL_v3_3INPUT.m
% Three-regime baseline simulation using mutually exclusive clamped inputs:
%   Hypo, NL, HL
%
% Regime logic:
%   Hypo loading   = Hypo high, NL low, HL very low
%   Normal loading = Hypo very low, NL high, HL very low
%   Hyper loading  = Hypo very low, NL low, HL high
%
% Perturbations start from the Hyper steady state.
%
% PAIRED-SEQUENTIAL VERSION:
% For each replicate r:
%   random initial -> Hyper steady state -> Perturbed steady state
%
% Delta is computed replicate-wise:
%   Delta_r = Perturbed_r - Hyper_r
%
% Parallelization:
%   - parpool starts workers
%   - parfor distributes independent replicates across workers

clc; close all; clear;

%% ===================== GLOBAL: SAVE ONLY (NO DISPLAY) ====================
set(groot,'DefaultFigureVisible','off');

%% ===================== GLOBAL STYLE ======================================
set(groot,'DefaultAxesTickLabelInterpreter','tex');
set(groot,'DefaultTextInterpreter','tex');
set(groot,'DefaultLegendInterpreter','tex');

set(groot,'DefaultAxesFontSize',24);
set(groot,'DefaultTextFontSize',24);
set(groot,'DefaultAxesLineWidth',1.8);
set(groot,'DefaultLineLineWidth',2.2);

set(groot,'DefaultAxesFontWeight','bold');
set(groot,'DefaultTextFontWeight','bold');

try, set(groot,'DefaultAxesTitleFontWeight','bold'); end %#ok<TRYNC>
try, set(groot,'DefaultLegendFontSize',22); end %#ok<TRYNC>

FS_TICK  = 24;
FS_LABEL = 28;
FS_TITLE = 30;
FS_STAR  = 20;

%% ===================== PARALLEL SETUP ====================================
USE_PARALLEL = true;
N_WORKERS    = [];   % [] = MATLAB default

if USE_PARALLEL
    ensure_parpool(N_WORKERS);
end

%% --------------- I/O & model --------------------------------------------
netfile  = 'MT_PRIMARY4_1.xlsx';
ODIR     = fullfile(pwd,'outputs_RESCUER_NEW4_1','HyperToNormal');  ensure_dir(ODIR);
FDIR     = fullfile(ODIR,'Figures_HyperToNormal');                ensure_dir(FDIR);
OUT_XLSX = fullfile(ODIR,'Sensitivity_FromHyper4_1.xlsx');

[Mact, Minh, NodeNames, NumOfNodes] = CreateMatrices_new(netfile);
NodeNames = string(NodeNames(:));

iHypo = find(NodeNames=="Hypo",1);
iNL   = find(NodeNames=="NL",1);
iHL   = find(NodeNames=="HL",1);
assert(~isempty(iHypo) && ~isempty(iNL) && ~isempty(iHL), 'Hypo/NL/HL nodes not found.');

ode_opts  = odeset('RelTol',1e-8,'AbsTol',1e-10,'MaxStep',0.5,'NonNegative',1:NumOfNodes);
try
    JP = (Mact~=0) | (Minh~=0) | eye(NumOfNodes)>0;
    ode_opts = odeset(ode_opts,'JPattern',JP);
catch
end
solver_handle = @ode45;

BASE_SEED = 1;

%% --------------- Baseline clamps & windows -------------------------------
% Use the same three-input chronic loading regime scheme as the baseline code
Hypo_hypo = 0.20; NL_hypo = 0.01; HL_hypo = 0.01;   % Hypo regime
Hypo_norm = 0.01; NL_norm = 0.80; HL_norm = 0.01;   % Normal regime
Hypo_hype = 0.01; NL_hype = 0.01; HL_hype = 0.80;   % Hyper regime

tspan_baseline = [0 100];
tspan_pert     = [0 100];

N_REP_BASE = 100;
N_PERM     = 1000;
ALPHA_Q    = 0.05;

%% --------------- Panels & markers for composite metrics ------------------
ANABOLIC_TX = string(["ACAN","COL2A1","TIMP3","PPARγ","GDF5","IGF1","SOX9","Bcl2"]);
CATABOLIC_TX= string(["MMP1","MMP3","MMP13","ADAMTS4/5","IL1β","IL6","IL8","TNF","iNOS","COX-2"]);

aIdx = map_custom_names(ANABOLIC_TX, NodeNames); aIdx = aIdx(aIdx>0);
cIdx = map_custom_names(CATABOLIC_TX, NodeNames); cIdx = cIdx(cIdx>0);

nNodes = NumOfNodes;

%% --------------- Group categories for readouts ---------------------------
group_categories = {
    {'Cytokines, chemokines, proteases & others', ...
        {'TNF','IL6','IL1β','IL8','CCL2','CXCL1','CXCL3','ADAMTS4/5','MMP1','MMP3','MMP13'}};

    {'Growth factors', ...
        {'TGFβ','VEGF','IGF1','BMP2','CCN2','GDF5','FGF2','FGF18','Wnt3a','Wnt5a'}};

    {'Transcription Factors', ...
        {'CREB','HIF-1α','HIF-2α','NF-κB','AP-1','FOXO','SOX9','ATF2','NFAT','RUNX2', ...
         'YAP/TAZ','MRTF-A','NRF2','HSF1','TonEBP','ELK1','PPARγ','CITED2'}};

    {'Ion channels & related', ...
        {'Ca2+os','Ca2+su','CaMKII','PKC-E','PKC-M','PLCγ-M','PLCγ-E','CaN','IP3','PLA2'}};

    {'Metabolic & related', ...
        {'LKB1','NAD+','AMPK','mTORC1','mTORC2','SIRT1','PI3K-M','PI3K-E','PIP3-M','PIP3-E', ...
         'PDK1-M','PDK1-E','AKT1-M','AKT1-E','GSK3B','ULK1','PTEN','PLD2','PGE2','COX-2', ...
         'PHD2','VHL','Rheb'}};

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

nGroups = size(group_categories,1);
READOUT = struct('name',[],'nodes',[],'idx',[]);

for g = 1:nGroups
    pair = group_categories{g};
    READOUT(g).name  = pair{1};
    READOUT(g).nodes = string(pair{2});
    READOUT(g).idx   = map_custom_names(READOUT(g).nodes, NodeNames);
    READOUT(g).idx   = READOUT(g).idx(READOUT(g).idx>0);
end

%% --------------- Candidate perturbations ---------------------------------
CATAB_KD_LIST  = string(["RhoA-E","PIEZO1","PI3K-E","FAK-E","ROS"]);
ANAB_UP_LIST   = string(["SOX9","PPARγ","HIF-1α","NRF2","IκBα"]);

CATAB_KD_LEVEL = 0;
ANAB_UP_LEVEL  = 1;

clamp_specs = struct('names',{},'values',{});

for k = 1:numel(CATAB_KD_LIST)
    clamp_specs(end+1) = struct('names', {{char(CATAB_KD_LIST(k))}}, ...
                                'values', CATAB_KD_LEVEL); %#ok<AGROW>
end

for u = 1:numel(ANAB_UP_LIST)
    clamp_specs(end+1) = struct('names', {{char(ANAB_UP_LIST(u))}}, ...
                                'values', ANAB_UP_LEVEL); %#ok<AGROW>
end

for u = 1:numel(ANAB_UP_LIST)
    for k = 1:numel(CATAB_KD_LIST)
        clamp_specs(end+1) = struct('names', {{char(ANAB_UP_LIST(u)), char(CATAB_KD_LIST(k))}}, ...
                                    'values', [ANAB_UP_LEVEL, CATAB_KD_LEVEL]); %#ok<AGROW>
    end
end

[clamp_specs, ~, ClampVars, ShortClampLabels] = normalize_clamps(clamp_specs, NodeNames);
ClampVars = matlab.lang.makeUniqueStrings(string(ClampVars), {}, namelengthmax);
nC = numel(clamp_specs);

%% --------------- Sim 1–3 baselines --------------------------------------
[SS_hyp_all,   SS_hyp_mu]   = run_reps_baseline_parallel3( ...
    Hypo_hypo, NL_hypo, HL_hypo, ...
    iHypo, iNL, iHL, NumOfNodes, tspan_baseline, N_REP_BASE, ...
    solver_handle, ode_opts, USE_PARALLEL, BASE_SEED+1000, Mact, Minh);

[SS_norm_all,  SS_norm_mu]  = run_reps_baseline_parallel3( ...
    Hypo_norm, NL_norm, HL_norm, ...
    iHypo, iNL, iHL, NumOfNodes, tspan_baseline, N_REP_BASE, ...
    solver_handle, ode_opts, USE_PARALLEL, BASE_SEED+2000, Mact, Minh); %#ok<ASGLU>

[SS_hyper_all, SS_hyper_mu] = run_reps_baseline_parallel3( ...
    Hypo_hype, NL_hype, HL_hype, ...
    iHypo, iNL, iHL, NumOfNodes, tspan_baseline, N_REP_BASE, ...
    solver_handle, ode_opts, USE_PARALLEL, BASE_SEED+3000, Mact, Minh);

BaselinesTbl = table(NodeNames(:), SS_hyp_mu(:), SS_norm_mu(:), SS_hyper_mu(:), ...
    'VariableNames', {'Node','Hypo','Normal','Hyper'});
writetable(BaselinesTbl, OUT_XLSX, 'Sheet','Baselines');

%% --------------- Baseline bar graph --------------------------------------
f_base = figure('Color','w','Position',[80 80 1800 900], 'Renderer','painters', 'Visible','off');
ax = axes('Parent', f_base);
hold(ax,'on'); grid(ax,'on'); box(ax,'on');

X = (1:NumOfNodes)';
Y = [SS_hyp_mu, SS_norm_mu, SS_hyper_mu];

hb = bar(ax, X, Y, 'grouped');
hb(1).FaceColor = [0.20 0.40 0.80];
hb(2).FaceColor = [0.25 0.70 0.40];
hb(3).FaceColor = [0.80 0.30 0.30];

set(ax,'XTick',X, ...
       'XTickLabel',NodeNames, ...
       'XTickLabelRotation',45, ...
       'TickLabelInterpreter','none', ...
       'FontWeight','bold', ...
       'FontSize',FS_TICK);

ylabel(ax,'Steady-state activation','FontWeight','bold','FontSize',FS_LABEL);
title(ax,'Baseline','FontWeight','bold','FontSize',FS_TITLE);
legend(ax,{'Hypo','Normal','Hyper'},'Location','northeastoutside','FontSize',FS_TICK,'FontWeight','bold');

export_png(f_base, fullfile(FDIR,'baseline_bars.png'));

%% --------------- PAIRED SEQUENTIAL PERTURBATIONS FROM HYPER -------------
baseClamp = nan(NumOfNodes,1);
baseClamp(iHypo) = Hypo_hype;
baseClamp(iNL)   = NL_hype;
baseClamp(iHL)   = HL_hype;

DELTA_mean_all = zeros(nNodes, nC);
DELTA_sd_all   = zeros(nNodes, nC);
Pmat_all       = ones(nNodes, nC);
Dmat_all       = zeros(nNodes, nC);
pertMean_all   = zeros(nNodes, nC);
FINAL          = zeros(NumOfNodes, nC);

BaseRep_all  = SS_hyper_all;
baseMean_all = mean(SS_hyper_all,1,'omitnan')';

[~, ~, Bal_base] = panel_metrics(SS_hyper_all, SS_hyper_mu, aIdx, cIdx);
sub    = unique([aIdx; cIdx]);
d_base = rowdist(SS_hyper_all(:,sub), SS_norm_mu(sub)');

Results = table('Size',[0 11], 'VariableTypes', ...
  {'string','double','double','double','double','double','double','double','double','double','double'}, ...
  'VariableNames', {'Node','ClampValue','Anab_mu','Catab_mu','Balance','Dist_to_Normal', ...
                    'd_Balance','p_Balance','q_Balance','d_Dist','p_Dist'});

for c = 1:nC
    clampVec = baseClamp;
    cv = nan(NumOfNodes,1);
    cv(clamp_specs(c).idx) = clamp_specs(c).values;
    mask = ~isnan(cv);
    clampVec(mask) = cv(mask);

    [PertRep, SS_mu] = run_reps_from_replicate_states_parallel( ...
        SS_hyper_all, clampVec, tspan_pert, solver_handle, ode_opts, USE_PARALLEL, Mact, Minh);

    FINAL(:,c) = mean(PertRep,1,'omitnan')';
    pertMean_all(:,c) = mean(PertRep,1,'omitnan')';

    DeltaRep = PertRep - BaseRep_all;

    for m = 1:nNodes
        dvec = DeltaRep(:,m);
        dvec = dvec(~isnan(dvec));

        x = PertRep(:,m);     x = x(~isnan(x));
        y = BaseRep_all(:,m); y = y(~isnan(y));

        if isempty(dvec) || isempty(x) || isempty(y)
            DELTA_mean_all(m,c) = 0;
            DELTA_sd_all(m,c)   = 0;
            Dmat_all(m,c)       = 0;
            Pmat_all(m,c)       = 1;
        else
            DELTA_mean_all(m,c) = mean(dvec,'omitnan');
            DELTA_sd_all(m,c)   = std(dvec,0,'omitnan');
            Dmat_all(m,c)       = cohen_d(x,y);
            Pmat_all(m,c)       = perm_p_two_sided(x,y,N_PERM);
        end
    end

    [~, ~, Bal_pert] = panel_metrics(PertRep, SS_mu, aIdx, cIdx);
    d_pert = rowdist(PertRep(:,sub), SS_norm_mu(sub)');

    dBal  = cohen_d(Bal_pert, Bal_base);
    pBal  = perm_p_two_sided(Bal_pert, Bal_base, N_PERM);
    dDist = mean(d_pert,'omitnan') - mean(d_base,'omitnan');
    pDist = perm_p_two_sided(d_pert, d_base, N_PERM);

    Results = [Results; {
        clamp_specs(c).label_node, clamp_specs(c).label_value, ...
        mean(SS_mu(aIdx),'omitnan'), mean(SS_mu(cIdx),'omitnan'), ...
        mean(Bal_pert,'omitnan'), mean(d_pert,'omitnan'), ...
        dBal, pBal, NaN, dDist, pDist
    }]; %#ok<AGROW>
end

Results.q_Balance = bh_fdr(Results.p_Balance);
Qmat_all = reshape(bh_fdr(Pmat_all(:)), size(Pmat_all));

%% --------------- Save tables ---------------------------------------------
writetable(Results, fullfile(ODIR,'Perturbation_Summary.csv'));

F = array2table(FINAL, 'VariableNames', cellstr(ClampVars), 'RowNames', cellstr(NodeNames));
writetable(F, OUT_XLSX, 'Sheet','FinalStates','WriteRowNames',true);

D = array2table(DELTA_mean_all, 'VariableNames', cellstr(ClampVars), 'RowNames', cellstr(NodeNames));
writetable(D, OUT_XLSX, 'Sheet','Sensitivity_Delta','WriteRowNames',true);

K  = nNodes * nC;
MarkerCol    = repmat(NodeNames(:), nC, 1);
ClampCol     = reshape(repmat(string(ShortClampLabels(:)).', nNodes, 1), [], 1);
PertMeanCol  = reshape(pertMean_all,    K, 1);
BaseMeanCol  = repmat(baseMean_all(:),  nC, 1);
DeltaMeanCol = reshape(DELTA_mean_all,  K, 1);
DeltaSDCol   = reshape(DELTA_sd_all,    K, 1);
CohenDCol    = reshape(Dmat_all,        K, 1);
PvalCol      = reshape(Pmat_all,        K, 1);
QvalCol      = reshape(Qmat_all,        K, 1);

MarkerStats = table( ...
    MarkerCol, ...
    ClampCol, ...
    PertMeanCol, ...
    BaseMeanCol, ...
    DeltaMeanCol, ...
    DeltaSDCol, ...
    CohenDCol, ...
    PvalCol, ...
    QvalCol, ...
    'VariableNames', {'Marker','Clamp','PertMean','BaseMean','DeltaMean','DeltaSD','Cohen_d','p_value','q_value'} ...
);
writetable(MarkerStats, OUT_XLSX, 'Sheet','MarkerStats');

%% --------------- Global composite figures --------------------------------
[vals, ix] = maxk(Results.d_Balance, min(10,height(Results)));
balance_labels = Results.Node(ix);
bar_top(vals, balance_labels, 'Balance', fullfile(FDIR,'top_rescuers_balance.png'));

[valsD, ixD] = maxk(-Results.d_Dist, min(10,height(Results)));
dist_labels = Results.Node(ixD);
bar_top(valsD, dist_labels, 'Distance to Normal', fullfile(FDIR,'top_rescuers_distance.png'));

%% --------------- Group-specific figures ----------------------------------
for g = 1:nGroups
    rows = READOUT(g).idx;
    if isempty(rows), continue; end

    gname  = READOUT(g).name;
    gnodes = NodeNames(rows);

    ttlG = gname;
    fnG  = fullfile(FDIR, sprintf('heatmap_%s.png', sanitize_name(gname)));
    heatmap_with_sig(DELTA_mean_all(rows,:), Qmat_all(rows,:), ShortClampLabels, gnodes, ttlG, fnG, ALPHA_Q, FS_TICK, FS_LABEL, FS_TITLE, FS_STAR);

    scores = mean(abs(DELTA_mean_all(rows,:)), 1, 'omitnan');

    if all(~isfinite(scores) | scores==0)
        warning('Group "%s": no finite Delta values. Skipping group plots.', gname);
        continue;
    end

    [valsG, ixG] = maxk(scores, min(10, numel(scores)));
    labelsG = ShortClampLabels(ixG);

    ttlBarG = gname;
    fnBarG  = fullfile(FDIR, sprintf('top_rescuers_%s.png', sanitize_name(gname)));
    bar_top(valsG, labelsG, ttlBarG, fnBarG, FS_TICK, FS_LABEL, FS_TITLE);

    for c = 1:nC
        meansG = DELTA_mean_all(rows, c);
        sdsG   = DELTA_sd_all(rows,   c);
        qcolG  = Qmat_all(rows,       c);

        ttlBarSD = sprintf('%s: %s', ...
                           gname, format_clamp_title(clamp_specs(c), ANAB_UP_LIST, CATAB_KD_LIST));

        fnBarSD  = fullfile(FDIR, sprintf('bar_withSD_%s_%s.png', ...
                    sanitize_name(gname), sanitize_name(ShortClampLabels{c})));

        bar_with_sd_and_sig(gnodes, meansG, sdsG, qcolG, ttlBarSD, fnBarSD, ALPHA_Q, FS_TICK, FS_LABEL, FS_TITLE, FS_STAR);
    end
end

fprintf('Saved paired sequential parallel results to %s\n', ODIR);

%% ===================== HELPERS ===========================================

function ensure_parpool(nWorkers)
try
    p = gcp('nocreate');
    if isempty(p)
        if isempty(nWorkers)
            parpool;
        else
            parpool(nWorkers);
        end
    else
        if ~isempty(nWorkers) && p.NumWorkers ~= nWorkers
            delete(p);
            parpool(nWorkers);
        end
    end
catch ME
    warning('Could not start parpool. Continuing without creating a new pool. Reason: %s', ME.message);
end
end

function [SS_all, SS_mu] = run_reps_baseline_parallel3(HypoVal,NLVal,HLVal,iHypo,iNL,iHL,N,tspan,nrep,solver,opts,useParallel,baseSeed,Mact,Minh)
SS_all = nan(nrep,N);

clampVec = nan(N,1);
clampVec(iHypo) = HypoVal;
clampVec(iNL)   = NLVal;
clampVec(iHL)   = HLVal;

if useParallel
    parfor r = 1:nrep
        rng(baseSeed + r, 'twister');
        x0 = rand(N,1);
        x0(iHypo) = HypoVal;
        x0(iNL)   = NLVal;
        x0(iHL)   = HLVal;
        [~,X] = solver(@(t,x) ODE_ClampMask_NaN(t,x,clampVec,Mact,Minh,N), tspan, x0, opts);
        SS_all(r,:) = X(end,:);
    end
else
    for r = 1:nrep
        rng(baseSeed + r, 'twister');
        x0 = rand(N,1);
        x0(iHypo) = HypoVal;
        x0(iNL)   = NLVal;
        x0(iHL)   = HLVal;
        [~,X] = solver(@(t,x) ODE_ClampMask_NaN(t,x,clampVec,Mact,Minh,N), tspan, x0, opts);
        SS_all(r,:) = X(end,:);
    end
end

SS_mu = mean(SS_all,1,'omitnan')';
end

function [SS_all, SS_mu] = run_reps_from_replicate_states_parallel(X0_all, clampVec, tspan, solver, opts, useParallel, Mact, Minh)
nrep = size(X0_all,1);
N    = size(X0_all,2);
SS_all = nan(nrep,N);

if useParallel
    parfor r = 1:nrep
        x0 = X0_all(r,:)';
        mask = ~isnan(clampVec);
        x0(mask) = clampVec(mask);
        [~,X] = solver(@(t,x) ODE_ClampMask_NaN(t,x,clampVec,Mact,Minh,N), tspan, x0, opts);
        SS_all(r,:) = X(end,:);
    end
else
    for r = 1:nrep
        x0 = X0_all(r,:)';
        mask = ~isnan(clampVec);
        x0(mask) = clampVec(mask);
        [~,X] = solver(@(t,x) ODE_ClampMask_NaN(t,x,clampVec,Mact,Minh,N), tspan, x0, opts);
        SS_all(r,:) = X(end,:);
    end
end

SS_mu = mean(SS_all,1,'omitnan')';
end

function dxdt = ODE_ClampMask(~,x,mask,vals,Mact,Minh,NumOfNodes)
x(mask) = vals(mask);
dxdt = ODESysFunS(0,x,Mact,Minh,NumOfNodes);
dxdt(mask) = 0;
end

function dxdt = ODE_ClampMask_NaN(~, x, clampVec, Mact, Minh, NumOfNodes)
mask = ~isnan(clampVec);
x(mask) = clampVec(mask);
dxdt = ODESysFunS(0, x, Mact, Minh, NumOfNodes);
dxdt(mask) = 0;
end

function [mA, mC, balance] = panel_metrics(SS_all, ~, aIdx, cIdx)
aVals = SS_all(:, aIdx);
cVals = SS_all(:, cIdx);
mA = mean(aVals,2,'omitnan');
mC = mean(cVals,2,'omitnan');
zA = (mA - mean(mA,'omitnan')) / max(std(mA,0,'omitnan'), eps);
zC = (mC - mean(mC,'omitnan')) / max(std(mC,0,'omitnan'), eps);
balance = zA - zC;
end

function idx = map_custom_names(names, NodeNames)
if isstring(names), names = cellstr(names); end
NN = normalize_for_match(NodeNames);
idx = zeros(numel(names),1);
for i=1:numel(names)
    c = normalize_for_match(string(names{i}));
    hit = find(strcmpi(c,NN),1);
    if ~isempty(hit), idx(i)=hit; end
end
end

function out = normalize_for_match(strs)
if isstring(strs), strs = cellstr(strs); end
if ischar(strs), strs = {strs}; end
out = cell(size(strs));
for i=1:numel(strs)
    s = lower(strs{i});
    s = strrep(s,'α','a'); s = strrep(s,'β','b'); s = strrep(s,'γ','g'); s = strrep(s,'κ','k');
    s = strrep(s,'–','-'); s = strrep(s,'—','-');
    s = regexprep(s,'\s+','');
    s = regexprep(s,'[^a-z0-9/\+\-_]','');
    out{i} = s;
end
end

function [specs, labels, vars, shortLabels] = normalize_clamps(specs, NodeNames)
labels = cell(1,numel(specs));
vars = labels;
shortLabels = labels;

for k = 1:numel(specs)
    if isfield(specs(k),'name') && ~isfield(specs(k),'names')
        specs(k).names = {specs(k).name};
    end
    if ~iscell(specs(k).names)
        specs(k).names = cellstr(specs(k).names);
    end
    if isscalar(specs(k).values)
        specs(k).values = repmat(specs(k).values,1,numel(specs(k).names));
    end

    idx = map_custom_names(string(specs(k).names), NodeNames);
    idx = idx(idx>0);

    specs(k).idx = idx(:)';
    specs(k).values = specs(k).values(1:numel(idx));

    if numel(idx)==1
        specs(k).label_node  = NodeNames(idx);
        specs(k).label_value = specs(k).values(1);
        labels{k}      = sprintf('%s=%.3g', specs(k).label_node, specs(k).label_value);
        shortLabels{k} = char(specs(k).label_node);
    else
        parts = cell(1,numel(idx));
        for i=1:numel(idx)
            parts{i} = sprintf('%s=%.3g', NodeNames(idx(i)), specs(k).values(i));
        end
        shortLabels{k} = char(strjoin(NodeNames(idx), '-'));
        specs(k).label_node  = string(shortLabels{k});
        specs(k).label_value = NaN;
        labels{k} = strjoin(parts,' + ');
    end

    vars{k} = regexprep(shortLabels{k},'[^A-Za-z0-9_]','_');
    if isempty(vars{k}), vars{k} = sprintf('Clamp_%d', k); end
end
end

function d = cohen_d(x,y)
x=x(:); y=y(:); x=x(~isnan(x)); y=y(~isnan(y));
nx=numel(x); ny=numel(y);
if nx<2||ny<2, d=0; return; end
s = sqrt(((nx-1)*var(x,0)+(ny-1)*var(y,0))/max(1,(nx+ny-2)));
if s<eps, d=0; else, d=(mean(x)-mean(y))/s; end
end

function p = perm_p_two_sided(a,b,nperm)
a=a(:); b=b(:); a=a(~isnan(a)); b=b(~isnan(b));
na=numel(a); nb=numel(b);
if na==0||nb==0, p=1; return; end
if nperm <= 0, p=1; return; end
obs=mean(a)-mean(b); pool=[a;b]; n=numel(pool); cnt=0;
for k=1:nperm
    idx=randperm(n); A=pool(idx(1:na)); B=pool(idx(na+1:end));
    if abs(mean(A)-mean(B))>=abs(obs), cnt=cnt+1; end
end
p=(cnt+1)/(nperm+1);
end

function q = bh_fdr(p)
p=p(:); [ps,ix]=sort(p); m=numel(p); qtmp=ps.*(m./(1:m))'; qtmp=min(1,qtmp);
for i=m-1:-1:1, qtmp(i)=min(qtmp(i),qtmp(i+1)); end
q=zeros(size(p)); q(ix)=qtmp;
end

function heatmap_with_sig(DELTA, Q, xlabels, ylabels, ttl, outpng, alpha_q, FS_TICK, FS_LABEL, FS_TITLE, FS_STAR)
f = figure('Color','w','Position',[50 50 2200 1100], 'Renderer','painters', 'Visible','off');
ax = axes('Parent', f);

imagesc(ax, DELTA);
colormap(ax, parula);
cb = colorbar(ax);
cb.FontSize = FS_TICK;
cb.FontWeight = 'bold';

hold(ax,'on');

set(ax,'XTick',1:numel(xlabels), ...
       'XTickLabel',xlabels, ...
       'XTickLabelRotation',45, ...
       'YTick',1:numel(ylabels), ...
       'YTickLabel',ylabels, ...
       'TickLabelInterpreter','none', ...
       'FontWeight','bold', ...
       'FontSize',FS_TICK);

vmax=max(abs(DELTA),[],'all');
if vmax<1e-9, vmax=1e-9; end
caxis(ax, [-vmax vmax]);

xlabel(ax,'Perturbations','FontWeight','bold','FontSize',FS_LABEL);
ylabel(ax,'Nodes','FontWeight','bold','FontSize',FS_LABEL);
title(ax, ttl, 'Interpreter','none','FontWeight','bold','FontSize',FS_TITLE);
box(ax,'on');

[nM,nC]=size(DELTA);
for j=1:nC
    for i=1:nM
        q=Q(i,j);
        if q<0.001, txt='***';
        elseif q<0.01, txt='**';
        elseif q<alpha_q, txt='*';
        else, txt=''; end
        if ~isempty(txt)
            text(ax, j, i, txt, 'Color','k', 'FontSize',FS_STAR, ...
                'HorizontalAlignment','center', 'FontWeight','bold');
        end
    end
end

drawnow;
export_png(f,outpng);
end

function bar_top(vals, labels, ttl, outpng, FS_TICK, FS_LABEL, FS_TITLE)
if nargin < 5
    FS_TICK  = 24;
    FS_LABEL = 28;
    FS_TITLE = 30;
end

f = figure('Color','w', 'Position',[100 100 1600 800], 'Renderer','painters', 'Visible','off');
ax = axes('Parent', f);

bar(ax, vals);
grid(ax, 'on'); box(ax, 'on');

labels = cellstr(string(labels));

set(ax, 'XTick', 1:numel(vals), ...
        'XTickLabel', labels, ...
        'XTickLabelRotation', 45, ...
        'TickLabelInterpreter', 'none', ...
        'FontWeight','bold', ...
        'FontSize',FS_TICK);

ylabel(ax, 'Mean |\Delta|', 'Interpreter','tex', 'FontWeight','bold', 'FontSize',FS_LABEL);
title(ax, ttl, 'Interpreter','none', 'FontWeight','bold', 'FontSize',FS_TITLE);

drawnow;
export_png(f, outpng);
end

function bar_with_sd_and_sig(labels, means, sds, qcol, ttl, outpng, alpha_q, FS_TICK, FS_LABEL, FS_TITLE, FS_STAR)
f = figure('Color','w','Position',[90 90 1700 850], 'Renderer','painters', 'Visible','off');
ax = axes('Parent', f); hold(ax,'on');

bar(ax, means, 'FaceColor',[0.25 0.45 0.85]);
grid(ax,'on'); box(ax,'on');

er = errorbar(ax, 1:numel(means), means, sds, sds, '.k');
er.LineStyle = 'none';
er.LineWidth = 1.8;

set(ax,'XTick',1:numel(labels), ...
       'XTickLabel',labels, ...
       'XTickLabelRotation',45, ...
       'TickLabelInterpreter','none', ...
       'FontWeight','bold', ...
       'FontSize',FS_TICK);

ylabel(ax,'\Delta (mean \pm SD)', 'Interpreter','tex', 'FontWeight','bold', 'FontSize',FS_LABEL);
title(ax, ttl, 'Interpreter','tex', 'FontWeight','bold', 'FontSize',FS_TITLE);

yoff = max([sds(:); 1e-6])*0.5 + 0.02;
for i = 1:numel(means)
    q = qcol(i);
    if q < 0.001
        txt = '***';
    elseif q < 0.01
        txt = '**';
    elseif q < alpha_q
        txt = '*';
    else
        txt = '';
    end
    if ~isempty(txt)
        text(ax, i, means(i) + sds(i) + yoff, txt, ...
            'HorizontalAlignment','center', ...
            'FontWeight','bold', ...
            'FontSize',FS_STAR, ...
            'Interpreter','none');
    end
end

drawnow;
export_png(f,outpng);
end

function out = format_clamp_title(spec, ANAB_UP_LIST, CATAB_KD_LIST)
names = string(spec.names);
parts = strings(0);

for i = 1:numel(names)
    if any(names(i) == ANAB_UP_LIST)
        parts(end+1) = "\uparrow " + strip_suffix_for_title(names(i)); %#ok<AGROW>
    end
end

for i = 1:numel(names)
    if any(names(i) == CATAB_KD_LIST)
        parts(end+1) = "\downarrow " + strip_suffix_for_title(names(i)); %#ok<AGROW>
    end
end

out = strjoin(parts, ', ');
end

function nm = strip_suffix_for_title(nm)
nm = string(nm);
nm = regexprep(nm,'-(E|M)$','');
end

function d = rowdist(A, centerRow)
D = bsxfun(@minus, A, centerRow); %#ok<BSXFUN>
d = sqrt(sum(D.^2,2));
end

function ensure_dir(d)
if ~isempty(d) && exist(d,'dir')~=7, mkdir(d); end
end

function s = sanitize_name(s)
s = regexprep(s,'[^A-Za-z0-9_]','_');
end

function export_png(h, filename)
try
    if ~isgraphics(h)
        warning('export_png: invalid graphics handle for %s', filename);
        return;
    end

    [p,~,ext] = fileparts(filename);
    if ~isempty(p) && exist(p,'dir')~=7
        mkdir(p);
    end
    if isempty(ext)
        filename = [filename '.png'];
    end

    if isprop(h,'Renderer')
        h.Renderer = 'painters';
    end
    set(h,'InvertHardcopy','off');

    drawnow;

    try
        print(h, filename, '-dpng', '-r600');
    catch
        try
            exportgraphics(h, filename, 'Resolution', 600);
        catch ME2
            warning('Could not export figure %s: %s', filename, ME2.message);
        end
    end

catch ME
    warning('export_png failed for %s: %s', filename, ME.message);
end

try
    if isgraphics(h)
        close(h);
    end
end
end