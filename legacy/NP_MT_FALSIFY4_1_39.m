%% MTN_AllInOne_EXCLUSIVE_BASELINES_ONLY_v4_FALSIFY_3INPUT_39.m
% Corrected version:
% 1) Uses the same exclusive 3-input chronic loading scheme as the main baseline code
% 2) Grouped bar plots sized for quarter-A4 PowerPoint composition
% 3) Large node sets automatically split across multiple figures
% 4) All nodes plotted correctly in all-node grouped bar plots
% 5) Robust clamping logic works even when a clamped value is 0
% 6) Falsification uses zero-based 95% bootstrap CI directional criteria
%    (no +/-0.02 tolerance bounds).
% 7) Baseline grouped bars have asymmetric 95% bootstrap CIs on mean activation.
% 8) Baseline three-condition and pairwise tests are nonparametric, with BH-FDR.
%    These analyze distributions across simulation initializations, not biology.

clc; close all; clear;

%% ------------------------- USER SETTINGS -------------------------------
netfile   = 'MT_PRIMARY4_1.xlsx';
out_xlsx  = 'Mechanotransduction_Results_39.xlsx';
fig_dir   = fullfile(pwd,'figures_FALSIF4_1_39');  ensure_dir(fig_dir);
out_dir   = fullfile(pwd,'outputs_FALSIF4_1_39');  ensure_dir(out_dir);

observed_names = {'Hypo','NL','HL','α5β1-FN','α5β1-Fs','αVβ6','αVβ3','TRPV4','PIEZO1'};

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
         'CAT','GPX1','SOD1','SOD2','HO-1','PHD2','VHL','Rheb','MitD'}};
    {'ECM anabolism & phenotype markers', ...
        {'COL2A1','COL1A1','COL10A1','ACAN','TIMP3'}};
    {'Mechanical stimuli & its receptors', ...
        {'Hypo','NL','HL','α5β1-FN','α5β1-Fs','αVβ3','αVβ6','SDC4-E','SDC4-M','TRPV4','PIEZO1'}};
    {'Cell survival, apoptosis & mitophagy/DNA-damage', ...
        {'Bcl2','BAX','CASP3','CASP9','BNIP3','GADD45','DRP1','MOMP'}};
    {'Oxidative-stress defense & proteostasis', ...
        {'HO-1','GPX1','SOD1','SOD2','CAT','HSP70','HSP27','ROS'}};
    {'MAPK & stress-activated kinases', ...
        {'RAS-M','RAS-E','RAF-M','RAF-E','MEK-M','MEK-E','ERK-M','ERK-E','MKK3/6','MKK4/7','JNK','p38','RSK','TAK1'}};
    % MST1/2 remains in the general all-node category, but NOT the benchmark panel.
    {'Rho GTPases, cytoskeletal & Hippo regulators', ...
        {'RhoA-M','RhoA-E','RAC1-M','RAC1-E','CDC42','ROCK-M','ROCK-E','PAK1','PKN1','FAK-M','FAK-E','MST1/2','LATS1/2'}};
};

tspan_baseline = [0 100];
rng(1);

% Three-input exclusive regime definitions
Hypo_hypo = 0.20; NL_hypo = 0.01; HL_hypo = 0.01;
Hypo_norm = 0.01; NL_norm = 0.80; HL_norm = 0.01;
% NOTE: NL_hype is 0.10 in this source file (Methods text states 0.01).
% Left unchanged to preserve the original simulation configuration.
Hypo_hype = 0.01; NL_hype = 0.10; HL_hype = 0.80;

YLIM_ACTIVATION = [0 1];

USE_PAR      = true;
PAR_WORKERS  = [];

% Revised experimental registry B01-B06 (literature-derived NL > HL).
% The legacy internal label 'anabolic' means an upward benchmark response;
% it does not imply that every marker is biologically anabolic.
ANABOLIC_LIST = {...
    'ACAN','COL2A1','TIMP3','Bcl2','SOD2','SOX9'};

% Revised experimental registry B07-B39 (literature-derived HL > NL).
% The legacy internal label 'catabolic' means a downward benchmark response.
% B38: GSK3B replaces MST1/2; this benchmark refers to active GSK3beta
% kinase (not its inhibitory Ser9 phosphorylation or total abundance).
% Reference: Li et al., Cell Physiol Biochem 2017;44:229-239,
% doi:10.1159/000484649 (2% versus 20% cyclic compression in rat NP cells).
CATABOLIC_LIST = {...
    'COL1A1','COL10A1','MMP3','MMP13','ADAMTS4/5',...
    'TNF','IL1β','IL6','IL8','CCL2','CXCL1','PGE2','COX-2',...
    'TGFβ','Wnt3a','BAX','CASP3','CASP9','HSF1','MRTF-A',...
    'YAP/TAZ','NF-κB','β-catenin','p38','JNK','ROS',...
    'PIEZO1','DRP1','PKN1','VEGF','CCN2',...
    'GSK3B','RUNX2'};
% Zero represents no difference in mean activation (Normal - Hyper).
% A benchmark passes only when its 95% bootstrap CI lies entirely
% on the experimentally expected side of zero.
FALS_ZERO  = 0;
FALS_ALPHA = 0.05;
FALS_NBOOT = 10000;

% Separate baseline-plot uncertainty from the directional-effect falsification.
% 95% bootstrap percentile intervals around the condition-specific MEAN activation.
BASE_NBOOT = 10000;
BASE_ALPHA = 0.05;
BASE_STATS_ALPHA = 0.05;   % BH-adjusted q-value threshold for reporting


%% ------------------------- LOAD NETWORK --------------------------------
[Mact, Minh, NodeNames, NumOfNodes] = CreateMatrices_new(netfile);
NodeLabels = ensure_unique_labels(NodeNames);

% Fail before simulation if any revised benchmark is absent or maps ambiguously.
validate_benchmark_registry(ANABOLIC_LIST, CATABOLIC_LIST, NodeNames);

idxHypo = find(ismember(NodeNames, {'Hypo'}));
idxNL   = find(ismember(NodeNames, {'NL'}));
idxHL   = find(ismember(NodeNames, {'HL'}));

if isempty(idxHypo)
    error('Node "Hypo" was not found in NodeNames. Add it to the topology file.');
end
if isempty(idxNL)
    error('Node "NL" was not found in NodeNames.');
end
if isempty(idxHL)
    error('Node "HL" was not found in NodeNames.');
end

clamp_idx = [idxHypo(:); idxNL(:); idxHL(:)];

observed = map_custom_names(observed_names, NodeNames);
if ~isempty(observed), observed_names = NodeNames(observed); end

varNames = matlab.lang.makeValidName(NodeNames, 'ReplacementStyle','delete');
varNames = matlab.lang.makeUniqueStrings(varNames, {}, namelengthmax);

assignin('base','Mact', Mact);
assignin('base','Minh', Minh);
assignin('base','NodeNames', NodeNames);
assignin('base','NumOfNodes', NumOfNodes);

%% ----------------- ODE OPTIONS -----------------------------------------
ode_opts = odeset('RelTol',1e-8,'AbsTol',1e-10,'MaxStep',0.5,'NonNegative',1:NumOfNodes);
USE_STIFF    = false;
USE_JPATTERN = true;
solver_handle = @ode45;
if USE_STIFF, solver_handle = @ode15s; end
if USE_JPATTERN
    S = (Mact~=0) | (Minh~=0);
    S = S | eye(size(S));
    ode_opts = odeset(ode_opts, 'JPattern', S);
end

%% ----------------- PARPOOL ---------------------------------------------
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
        warning('Parallel pool unavailable (%s). Running serially.', ME.message);
        USE_PAR = false;
    end
end

%% ----------------- SIM HELPERS -----------------------------------------
simulate_replicates = @(HypoVal,NLVal,HLVal,tspan,numiter) run_replicates_clamped3(...
    HypoVal,NLVal,HLVal,clamp_idx,NumOfNodes,tspan,numiter,solver_handle,ode_opts,BASE_VARS);

%% ------------------ PART A: THREE EXCLUSIVE BASELINES ------------------
fprintf('Running baselines...\n');

[SS_hypo_all, SS_hypo_mu, SS_hypo_sd] = simulate_replicates(Hypo_hypo, NL_hypo, HL_hypo, tspan_baseline, 100);
[SS_norm_all, SS_norm_mu, SS_norm_sd] = simulate_replicates(Hypo_norm, NL_norm, HL_norm, tspan_baseline, 100);
[SS_hype_all, SS_hype_mu, SS_hype_sd] = simulate_replicates(Hypo_hype, NL_hype, HL_hype, tspan_baseline, 100);

BaselinesTbl = table(NodeNames(:), SS_hypo_mu(:), SS_norm_mu(:), SS_hype_mu(:), ...
  'VariableNames', {'Node','Hypo','Normal','Hyper'});
writetable(BaselinesTbl, out_xlsx, 'Sheet','Baselines');

VarNameMap = table(cellstr(varNames(:)), NodeNames(:), ...
    'VariableNames', {'TableVarName','OriginalNodeName'});
writetable(VarNameMap, out_xlsx, 'Sheet','VarNameMap');

ClampDiagnostic = table( ...
    ["Hypo";"Normal";"Hyper"], ...
    [SS_hypo_mu(idxHypo); SS_norm_mu(idxHypo); SS_hype_mu(idxHypo)], ...
    [SS_hypo_mu(idxNL);   SS_norm_mu(idxNL);   SS_hype_mu(idxNL)], ...
    [SS_hypo_mu(idxHL);   SS_norm_mu(idxHL);   SS_hype_mu(idxHL)], ...
    'VariableNames', {'Condition','Hypo_mean','NL_mean','HL_mean'});
writetable(ClampDiagnostic, out_xlsx, 'Sheet','ClampDiagnostic');

save(fullfile(out_dir,'Baselines_Replicates.mat'), ...
     'NodeNames','SS_hypo_all','SS_norm_all','SS_hype_all', ...
     'SS_hypo_mu','SS_norm_mu','SS_hype_mu', ...
     'SS_hypo_sd','SS_norm_sd','SS_hype_sd','-v7.3');

RepLong = build_replicates_long(NodeNames, SS_hypo_all, SS_norm_all, SS_hype_all);
writetable(RepLong, fullfile(out_dir,'Replicates_Long.csv'));

BaseSum = table(NodeNames(:), ...
    SS_hypo_mu(:), SS_hypo_sd(:), ...
    SS_norm_mu(:), SS_norm_sd(:), ...
    SS_hype_mu(:), SS_hype_sd(:), ...
    'VariableNames',{'Node','Hypo_Mean','Hypo_SD','Normal_Mean','Normal_SD','Hyper_Mean','Hyper_SD'});
writetable(BaseSum, fullfile(out_dir,'Baseline_Summary.csv'));
writetable(BaseSum, out_xlsx, 'Sheet','Baseline_Summary');

%% ------------- PART B: Pairwise differences ----------------------------
diff_Norm_vs_Hypo  = SS_norm_mu - SS_hypo_mu;
diff_Hyper_vs_Norm = SS_hype_mu - SS_norm_mu;
diff_Hyper_vs_Hypo = SS_hype_mu - SS_hypo_mu;

pooled_sd = @(A,B) sqrt((var(A,0,1,'omitnan')*(size(A,1)-1) + var(B,0,1,'omitnan')*(size(B,1)-1)) ...
                         ./ max(1,(size(A,1)+size(B,1)-2)));

sd_Norm_vs_Hypo = pooled_sd(SS_norm_all, SS_hypo_all);
sd_Hyper_vs_Norm = pooled_sd(SS_hype_all, SS_norm_all);
sd_Hyper_vs_Hypo = pooled_sd(SS_hype_all, SS_hypo_all);

d_Norm_vs_Hypo  = safe_divide(diff_Norm_vs_Hypo, sd_Norm_vs_Hypo);
d_Hyper_vs_Norm = safe_divide(diff_Hyper_vs_Norm, sd_Hyper_vs_Norm);
d_Hyper_vs_Hypo = safe_divide(diff_Hyper_vs_Hypo, sd_Hyper_vs_Hypo);

PairwiseAll = table(NodeNames(:), ...
    diff_Norm_vs_Hypo(:), sd_Norm_vs_Hypo(:), d_Norm_vs_Hypo(:), ...
    diff_Hyper_vs_Norm(:), sd_Hyper_vs_Norm(:), d_Hyper_vs_Norm(:), ...
    diff_Hyper_vs_Hypo(:), sd_Hyper_vs_Hypo(:), d_Hyper_vs_Hypo(:), ...
    'VariableNames',{'Node', ...
      'Delta_Norm_vs_Hypo','SDpooled_Norm_vs_Hypo','d_Norm_vs_Hypo', ...
      'Delta_Hyper_vs_Normal','SDpooled_Hyper_vs_Normal','d_Hyper_vs_Normal', ...
      'Delta_Hyper_vs_Hypo','SDpooled_Hyper_vs_Hypo','d_Hyper_vs_Hypo'});

writetable(PairwiseAll, out_xlsx, 'Sheet','Pairwise_All');
writetable(PairwiseAll, fullfile(out_dir,'Pairwise_All.csv'));

PairwiseTbl = table(NodeNames(:), diff_Norm_vs_Hypo(:), diff_Hyper_vs_Hypo(:), ...
  d_Norm_vs_Hypo(:), d_Hyper_vs_Hypo(:), ...
  'VariableNames',{'Node','Diff_Norm_vs_Hypo','Diff_Hyper_vs_Hypo','d_Norm_vs_Hypo','d_Hyper_vs_Hypo'});
writetable(PairwiseTbl, out_xlsx, 'Sheet','Pairwise_vs_Hypo');

%% ======================= PART B2: FALSIFICATION ========================
idxAnab = map_custom_names(ANABOLIC_LIST, NodeNames);
idxCat  = map_custom_names(CATABOLIC_LIST, NodeNames);
assert(numel(idxAnab)==6 && numel(idxCat)==33 && ...
    isempty(intersect(idxAnab,idxCat)), ...
    'Revised benchmark mapping must contain 6 distinct NL-higher and 33 distinct HL-higher nodes.');

Classes = strings(NumOfNodes,1); Classes(:) = "neutral";
Classes(idxAnab) = "anabolic";
Classes(idxCat)  = "catabolic";

MeanL = SS_hypo_mu(:);
MeanN = SS_norm_mu(:);
MeanH = SS_hype_mu(:);

DiffNH = MeanN - MeanH;
[ciLow_NH, ciHigh_NH] = bootstrap_diff_CI(SS_norm_all, SS_hype_all, FALS_NBOOT, FALS_ALPHA);

Pass_NH = false(NumOfNodes,1);
for i = 1:NumOfNodes
    if Classes(i)=="anabolic"
        Pass_NH(i) = (ciLow_NH(i) > FALS_ZERO);
    elseif Classes(i)=="catabolic"
        Pass_NH(i) = (ciHigh_NH(i) < FALS_ZERO);
    else
        Pass_NH(i) = true;
    end
end

FalsTbl_NH = table(string(NodeNames(:)), Classes, MeanN, MeanH, DiffNH, ...
    ciLow_NH(:), ciHigh_NH(:), Pass_NH, ...
    'VariableNames', {'Node','Class','Mean_Normal','Mean_Hyper', ...
                      'Delta_NminusH','CI_low','CI_high','Pass'});

maskTest = (Classes=="anabolic") | (Classes=="catabolic");
FalsTbl_NH_Export = FalsTbl_NH(maskTest,:);

sum_total_NH = sum(maskTest);
sum_pass_NH  = sum(FalsTbl_NH.Pass(maskTest));
sum_fail_NH  = sum(~FalsTbl_NH.Pass(maskTest));
sum_anab_NH  = [sum(Classes=="anabolic"), sum(Pass_NH & Classes=="anabolic")];
sum_cat_NH   = [sum(Classes=="catabolic"), sum(Pass_NH & Classes=="catabolic")];

FalsSummary_NH = table(...
    sum_total_NH, sum_pass_NH, sum_fail_NH, ...
    sum_anab_NH(1), sum_anab_NH(2), safe_divide(sum_anab_NH(2), max(1,sum_anab_NH(1))), ...
    sum_cat_NH(1),  sum_cat_NH(2),  safe_divide(sum_cat_NH(2),  max(1,sum_cat_NH(1))), ...
    FALS_ZERO, FALS_ALPHA, FALS_NBOOT, ...
    'VariableNames',{'TotalRules','Passed','Failed', ...
    'Anabolic_N','Anabolic_Passed','Anabolic_PassRate', ...
    'Catabolic_N','Catabolic_Passed','Catabolic_PassRate', ...
    'ZeroReference','Alpha','NBoot'});

writetable(FalsTbl_NH_Export, fullfile(out_dir,'Falsification_Normal_vs_Hyper.csv'));
try, writetable(FalsTbl_NH_Export, out_xlsx, 'Sheet','Falsification_NvH'); catch, end
writetable(FalsSummary_NH, fullfile(out_dir,'Falsification_Summary_NvH.csv'));
try, writetable(FalsSummary_NH, out_xlsx, 'Sheet','Falsification_Summary_NvH'); catch, end

DiffNL = MeanN - MeanL;
[ciLow_NL, ciHigh_NL] = bootstrap_diff_CI(SS_norm_all, SS_hypo_all, FALS_NBOOT, FALS_ALPHA);

Pass_NL = false(NumOfNodes,1);
for i = 1:NumOfNodes
    if Classes(i)=="anabolic"
        Pass_NL(i) = (ciLow_NL(i) > FALS_ZERO);
    elseif Classes(i)=="catabolic"
        Pass_NL(i) = (ciHigh_NL(i) < FALS_ZERO);
    else
        Pass_NL(i) = true;
    end
end

FalsTbl_NL = table(string(NodeNames(:)), Classes, MeanL, MeanN, DiffNL, ...
    ciLow_NL(:), ciHigh_NL(:), Pass_NL, ...
    'VariableNames', {'Node','Class','Mean_Hypo','Mean_Normal', ...
                      'Delta_NminusL','CI_low','CI_high','Pass'});

maskTest_NL = maskTest;
FalsTbl_NL_Export = FalsTbl_NL(maskTest_NL,:);

sum_total_NL = sum(maskTest_NL);
sum_pass_NL  = sum(FalsTbl_NL.Pass(maskTest_NL));
sum_fail_NL  = sum(~FalsTbl_NL.Pass(maskTest_NL));
sum_anab_NL  = [sum(Classes=="anabolic"), sum(Pass_NL & Classes=="anabolic")];
sum_cat_NL   = [sum(Classes=="catabolic"), sum(Pass_NL & Classes=="catabolic")];

FalsSummary_NL = table(...
    sum_total_NL, sum_pass_NL, sum_fail_NL, ...
    sum_anab_NL(1), sum_anab_NL(2), safe_divide(sum_anab_NL(2), max(1,sum_anab_NL(1))), ...
    sum_cat_NL(1),  sum_cat_NL(2),  safe_divide(sum_cat_NL(2),  max(1,sum_cat_NL(1))), ...
    FALS_ZERO, FALS_ALPHA, FALS_NBOOT, ...
    'VariableNames',{'TotalRules','Passed','Failed', ...
    'Anabolic_N','Anabolic_Passed','Anabolic_PassRate', ...
    'Catabolic_N','Catabolic_Passed','Catabolic_PassRate', ...
    'ZeroReference','Alpha','NBoot'});

writetable(FalsTbl_NL_Export, fullfile(out_dir,'Falsification_Normal_vs_Hypo.csv'));
try, writetable(FalsTbl_NL_Export, out_xlsx, 'Sheet','Falsification_NvL'); catch, end
writetable(FalsSummary_NL, fullfile(out_dir,'Falsification_Summary_NvL.csv'));
try, writetable(FalsSummary_NL, out_xlsx, 'Sheet','Falsification_Summary_NvL'); catch, end

%% ===================== MANUSCRIPT-QUALITY VISUALIZATION ================
COL_HYPO  = [0 0 1];
COL_NORM  = [0 0.60 0];
COL_HYPER = [1 0 0];
COL_ANAB  = [0.15 0.65 0.15];
COL_CATA  = [0.85 0.25 0.25];
COL_NEUT  = [0.55 0.55 0.55];

FS = get_pub_style();

%% Forest plot: Normal vs Hyper (tested nodes on x-axis)
% Vertical 95% bootstrap confidence intervals around each mean difference.
% Zero is the only reference: no +/- tolerance threshold lines.
testedIdx_NH = find(maskTest);
[~, order_NH] = sort(DiffNH(testedIdx_NH), 'descend');
ordIdx_NH = testedIdx_NH(order_NH);

if ~isempty(ordIdx_NH)
    f = create_pub_figure('forest');
    ax = axes('Parent',f); hold(ax,'on');
    x = 1:numel(ordIdx_NH);
    capHalfWidth = 0.16;  % cap width in node-position units

    for k = 1:numel(ordIdx_NH)
        i = ordIdx_NH(k);
        % Vertical CI whisker and horizontal end caps.
        line(ax,[x(k) x(k)], [ciLow_NH(i) ciHigh_NH(i)], ...
            'LineWidth',1.5,'Color',[0.55 0.55 0.55]);
        line(ax,x(k)+[-capHalfWidth capHalfWidth], [ciLow_NH(i) ciLow_NH(i)], ...
            'LineWidth',1.2,'Color',[0.55 0.55 0.55]);
        line(ax,x(k)+[-capHalfWidth capHalfWidth], [ciHigh_NH(i) ciHigh_NH(i)], ...
            'LineWidth',1.2,'Color',[0.55 0.55 0.55]);
    end
    for k = 1:numel(ordIdx_NH)
        i = ordIdx_NH(k);
        col = COL_ANAB;
        if Classes(i)=="catabolic", col = COL_CATA; end
        mk = 'o';
        if ~Pass_NH(i), mk = 's'; end
        plot(ax,x(k),DiffNH(i),mk,'MarkerFaceColor',col,'MarkerEdgeColor','k', ...
            'MarkerSize',6.0,'LineWidth',0.9);
    end

    % Horizontal line at zero because mean difference is on the y-axis.
    yline(ax,FALS_ZERO,'k-','LineWidth',1.1);

    labels = make_readable_labels(NodeNames(ordIdx_NH),28,false);
    set(ax,'XTick',x,'XTickLabel',labels,'XTickLabelRotation',90);
    xlim(ax,[0.4 numel(ordIdx_NH)+0.6]);
    xlabel(ax,'Tested nodes','FontSize',FS.label,'FontWeight','bold');
    ylabel(ax,'Mean difference (Normal - Hyper)','FontSize',FS.label,'FontWeight','bold');
    title(ax,sprintf('Falsification forest plot (Normal - Hyper, %.0f%% CI)',(1-FALS_ALPHA)*100), ...
        'FontSize',FS.title,'FontWeight','bold');
    apply_pub_axes(ax,FS,'forest');
    ax.TickLabelInterpreter = 'none';

    % Dashed grid lines for the forest plot
    ax.XGrid = 'on';          % Vertical grid lines at each node
    ax.YGrid = 'on';          % Horizontal grid lines
    ax.GridLineStyle = '--';  % Dashed lines
    ax.GridAlpha = 0.20;      % Subtle grid appearance
    savefig_safe(f,fullfile(out_dir,'falsification_forest_NormalMinusHyper'));
    close(f);
end

%% Pass/fail strip: Normal vs Hyper
vals_NH = nan(numel(NodeNames),1);
vals_NH(Classes=="anabolic")  = (Pass_NH(Classes=="anabolic")*2 - 1);
vals_NH(Classes=="catabolic") = (Pass_NH(Classes=="catabolic")*2 - 1);

[~, oNH] = sort(DiffNH(testedIdx_NH), 'descend');
ord_NH = testedIdx_NH(oNH);
if ~isempty(ord_NH)
    f = create_pub_figure('quarterwide');
    ax = axes('Parent',f);
    imagesc(ax, vals_NH(ord_NH)');
    colormap(ax, [0.85 0.25 0.25; 0.96 0.96 0.96; 0.20 0.65 0.20]);
    caxis(ax,[-1 1]);
    labels = make_readable_labels(NodeNames(ord_NH), 18, false);
    idxShow = choose_tick_idx(numel(ord_NH), 18);
    set(ax,'YTick',1,'YTickLabel',{'Pass/Fail'}, ...
        'XTick',idxShow,'XTickLabel',labels(idxShow), 'XTickLabelRotation',45);
    title(ax,'Pass/fail strip: Normal vs Hyper','FontSize',FS.title,'FontWeight','bold');
    apply_pub_axes(ax,FS,'strip');
    ax.TickLabelInterpreter = 'none';
    savefig_safe(f, fullfile(out_dir,'falsification_passfail_strip_NvH'));
    close(f);
end

%% Summary pass-rate bars
f = create_pub_figure('quarter');
ax = axes('Parent',f);
cats = categorical({'NL-higher','HL-higher','Overall'});
vals = [sum(Pass_NH & Classes=="anabolic") / max(1,sum(Classes=="anabolic")), ...
        sum(Pass_NH & Classes=="catabolic") / max(1,sum(Classes=="catabolic")), ...
        sum(Pass_NH(maskTest)) / max(1,sum(maskTest))];
b = bar(ax,cats,vals,0.65,'FaceColor','flat','EdgeColor','none');
b.CData = [COL_ANAB; COL_CATA; 0.35 0.35 0.35];
ylim(ax,[0 1]);
ylabel(ax,'Pass rate','FontSize',FS.label,'FontWeight','bold');
title(ax,'Falsification pass rates (Normal vs Hyper)','FontSize',FS.title,'FontWeight','bold');
apply_pub_axes(ax,FS,'bar');
savefig_safe(f, fullfile(out_dir,'falsification_pass_rates_NvH'));
close(f);

%% Robustness map: Normal vs Hyper
f = create_pub_figure('quarter');
ax = axes('Parent',f); hold(ax,'on');
box(ax,'on');

mA = (Classes=="anabolic");
mC = (Classes=="catabolic");
mN = ~(mA|mC);
mAP_NH = mA & Pass_NH; mAF_NH = mA & ~Pass_NH;
mCP_NH = mC & Pass_NH; mCF_NH = mC & ~Pass_NH;

scatter(ax,MeanH(mAP_NH), MeanN(mAP_NH), 55, COL_ANAB, 'o', 'filled', 'MarkerEdgeColor','k', 'LineWidth',0.8);
scatter(ax,MeanH(mAF_NH), MeanN(mAF_NH), 75, COL_ANAB, 's', 'filled', 'MarkerEdgeColor','k', 'LineWidth',0.9);
scatter(ax,MeanH(mCP_NH), MeanN(mCP_NH), 55, COL_CATA, 'o', 'filled', 'MarkerEdgeColor','k', 'LineWidth',0.8);
scatter(ax,MeanH(mCF_NH), MeanN(mCF_NH), 75, COL_CATA, 's', 'filled', 'MarkerEdgeColor','k', 'LineWidth',0.9);
scatter(ax,MeanH(mN), MeanN(mN), 28, COL_NEUT, 'o', 'filled', 'MarkerEdgeColor',[0.3 0.3 0.3], 'MarkerFaceAlpha',0.75);

lims = [min([MeanH;MeanN]) max([MeanH;MeanN])];
if isempty(lims) || any(isnan(lims)), lims = [0 1]; end
lims = expand_limits(lims, 0.04);
plot(ax,lims,lims,'k--','LineWidth',1.0);

xlabel(ax,'Hyper steady-state activation','FontSize',FS.label,'FontWeight','bold');
ylabel(ax,'Normal steady-state activation','FontSize',FS.label,'FontWeight','bold');
title(ax,'Robustness map: Normal vs Hyper','FontSize',FS.title,'FontWeight','bold');
axis(ax,'square'); xlim(ax,lims); ylim(ax,lims);
legend(ax,{'NL-higher (pass)','NL-higher (fail)','HL-higher (pass)','HL-higher (fail)','Neutral'}, ...
       'Location','southoutside','Orientation','horizontal','Box','off', ...
       'FontSize',FS.legend,'FontWeight','bold');
apply_pub_axes(ax,FS,'default');
savefig_safe(f, fullfile(out_dir,'robustness_map_Normal_vs_Hyper_enhanced'));
close(f);

%% Directional effects bar: Normal vs Hyper
[~, ord] = sort(DiffNH(testedIdx_NH), 'descend');
ordIdx = testedIdx_NH(ord);

if ~isempty(ordIdx)

    maxPerFig = 50;   % 🔥 adjust if needed (30–50 is good)
    N = numel(ordIdx);
    nPages = ceil(N / maxPerFig);

    for pg = 1:nPages

        idx1 = (pg-1)*maxPerFig + 1;
        idx2 = min(pg*maxPerFig, N);
        idxp = idx1:idx2;

        % Use the same wide layout and vertical node labels as the
        % falsification forest plot to avoid overlap in the directional bar plot.
        f = create_pub_figure('forest');
        ax = axes('Parent',f); hold(ax,'on');

        % Data for this page
        vals = DiffNH(ordIdx(idxp));

        bh = bar(ax, vals, 0.92, 'FaceColor','flat','EdgeColor','none');

        % Colors
        C = zeros(numel(idxp),3);
        C(Classes(ordIdx(idxp))=="anabolic",:) = repmat(COL_ANAB, sum(Classes(ordIdx(idxp))=="anabolic"),1);
        C(Classes(ordIdx(idxp))=="catabolic",:) = repmat(COL_CATA, sum(Classes(ordIdx(idxp))=="catabolic"),1);
        bh.CData = C;

        % Reference line at zero only.
        yline(ax,FALS_ZERO,'k-','LineWidth',1.0);

        % Match the node-label treatment used in the bootstrap 95% CI forest plot.
        labels = make_readable_labels(NodeNames(ordIdx(idxp)), 28, false);
        x = 1:numel(idxp);
        set(ax, ...
            'XTick',x, ...
            'XTickLabel',labels, ...
            'XTickLabelRotation',90);
        xlim(ax,[0.4 numel(idxp)+0.6]);

        ylabel(ax,'\Delta (Normal - Hyper)','FontSize',FS.label,'FontWeight','bold');
        xlabel(ax,'Tested nodes','FontSize',FS.label,'FontWeight','bold');

        if nPages > 1
            title_here = sprintf('Directional effects (Normal - Hyper) (%d/%d)', pg, nPages);
            savepath_here = fullfile(out_dir, sprintf('directional_effects_bar_NvH_part_%02d', pg));
        else
            title_here = 'Directional effects by node (Normal - Hyper)';
            savepath_here = fullfile(out_dir,'directional_effects_bar_NvH');
        end

        title(ax, title_here, 'FontSize',FS.title,'FontWeight','bold');

        apply_pub_axes(ax,FS,'forest');
        ax.TickLabelInterpreter = 'none';
        ax.XGrid = 'on';
        ax.YGrid = 'on';
        ax.GridLineStyle = '--';
        ax.GridAlpha = 0.20;

        savefig_safe(f, savepath_here);
        close(f);

    end
end

%% Stability map: CV Normal vs Hyper
f = create_pub_figure('quarter');
ax = axes('Parent',f); hold(ax,'on');
sdN = std(SS_norm_all,0,1,'omitnan');
sdH = std(SS_hype_all,0,1,'omitnan');
cvN = sdN ./ max(MeanN', realmin);
cvH = sdH ./ max(MeanH', realmin);

scatter(ax,cvH(Classes=="anabolic"),  cvN(Classes=="anabolic"),  50, COL_ANAB, 'o', 'filled', 'MarkerEdgeColor','k', 'LineWidth',0.8, 'MarkerFaceAlpha',0.9);
scatter(ax,cvH(Classes=="catabolic"), cvN(Classes=="catabolic"), 50, COL_CATA, 'o', 'filled', 'MarkerEdgeColor','k', 'LineWidth',0.8, 'MarkerFaceAlpha',0.9);
scatter(ax,cvH(~(Classes=="anabolic" | Classes=="catabolic")), ...
           cvN(~(Classes=="anabolic" | Classes=="catabolic")), ...
           25, COL_NEUT, 'o', 'filled', 'MarkerEdgeColor',[0.3 0.3 0.3], 'MarkerFaceAlpha',0.7);

limsCV = [0 max([cvH cvN],[],'all','omitnan')];
if isempty(limsCV) || any(isnan(limsCV)), limsCV = [0 0.01]; end
limsCV = expand_limits(limsCV, 0.05);
plot(ax,limsCV,limsCV,'k--','LineWidth',1.0);

xlabel(ax,'CV (Hyper)','FontSize',FS.label,'FontWeight','bold');
ylabel(ax,'CV (Normal)','FontSize',FS.label,'FontWeight','bold');
title(ax,'Stability map: CV (Normal vs Hyper)','FontSize',FS.title,'FontWeight','bold');
axis(ax,'square'); xlim(ax,limsCV); ylim(ax,limsCV);
legend(ax,{'NL-higher','HL-higher','Neutral'}, 'Location','southoutside', ...
    'Orientation','horizontal','Box','off','FontSize',FS.legend,'FontWeight','bold');
apply_pub_axes(ax,FS,'default');
savefig_safe(f, fullfile(out_dir,'stability_map_CV_NvH_enhanced'));
close(f);

%% -------- BASELINE UNCERTAINTY AND NONPARAMETRIC COMPARISONS ----------
% Use the SAME replicate outcomes as the existing baseline plots.
% Bootstrap mean CIs allow asymmetric bounds for skewed distributions.
% These CIs describe uncertainty in estimated condition-specific MEANS, not
% the spread of individual outcomes or biological experimental variation.
fprintf('Computing baseline mean 95%% bootstrap confidence intervals...\n');
[CI_hypo_lo, CI_hypo_hi] = bootstrap_mean_CI(SS_hypo_all, BASE_NBOOT, BASE_ALPHA);
[CI_norm_lo, CI_norm_hi] = bootstrap_mean_CI(SS_norm_all, BASE_NBOOT, BASE_ALPHA);
[CI_hype_lo, CI_hype_hi] = bootstrap_mean_CI(SS_hype_all, BASE_NBOOT, BASE_ALPHA);

CI_LOW_ALL  = [CI_hypo_lo(:), CI_norm_lo(:), CI_hype_lo(:)];
CI_HIGH_ALL = [CI_hypo_hi(:), CI_norm_hi(:), CI_hype_hi(:)];
MEANS_ALL   = [SS_hypo_mu(:), SS_norm_mu(:), SS_hype_mu(:)];

BaselineCI = table(string(NodeNames(:)), ...
    SS_hypo_mu(:), CI_hypo_lo(:), CI_hypo_hi(:), ...
    SS_norm_mu(:), CI_norm_lo(:), CI_norm_hi(:), ...
    SS_hype_mu(:), CI_hype_lo(:), CI_hype_hi(:), ...
    'VariableNames',{'Node','Hypo_Mean','Hypo_CI_Low','Hypo_CI_High', ...
    'NL_Mean','NL_CI_Low','NL_CI_High', ...
    'HL_Mean','HL_CI_Low','HL_CI_High'});
writetable(BaselineCI, fullfile(out_dir,'Baseline_Means_95CI.csv'));
try, writetable(BaselineCI,out_xlsx,'Sheet','Baseline_Means_95CI'); catch, end

%% -------- SUPPLEMENTARY FIGURES A AND B (ORIGINAL v12 VISUAL STYLE) -----
% Uses this run's 39-benchmark falsification and all-node baseline CI tables directly.
plot_falsification_AB_v12(FalsTbl_NH_Export, BaselineCI, out_dir);

% Independent randomized initializations under each condition justify
% unpaired, rank-based comparisons of the simulated outcome distributions.
% KW assesses all 3 conditions; rank-sum provides pairwise comparisons.
% BH-FDR q-values are calculated across ALL nodes within each test family.
fprintf('Computing rank-based baseline comparisons and BH-FDR...\n');
BaselineStats = baseline_nonparametric_stats(SS_hypo_all, SS_norm_all, ...
    SS_hype_all, NodeNames, clamp_idx);
writetable(BaselineStats, fullfile(out_dir,'Baseline_Nonparametric_Statistics.csv'));
try, writetable(BaselineStats,out_xlsx,'Sheet','Baseline_Statistics'); catch, end

nKW = sum(BaselineStats.q_KW_3group < BASE_STATS_ALPHA);
nNH = sum(BaselineStats.q_NL_HL < BASE_STATS_ALPHA);
nEligible = sum(isfinite(BaselineStats.p_KW_3group));
fprintf('Baseline: KW q<%.3f in %d/%d nodes; NL-HL rank-sum q<%.3f in %d/%d nodes.\n', ...
    BASE_STATS_ALPHA,nKW,nEligible, ...
    BASE_STATS_ALPHA,nNH,nEligible);
summaryPath = fullfile(out_dir,'Baseline_Statistical_Summary.txt');
fid = fopen(summaryPath,'w');
if fid ~= -1
    fprintf(fid,'Baseline comparisons of RNM outcomes across randomized initial conditions\n');
    fprintf(fid,'Mean activation error bars: percentile bootstrap %.1f%% CIs, %d resamples/condition\n', ...
        100*(1-BASE_ALPHA),BASE_NBOOT);
    fprintf(fid,'Kruskal-Wallis: three conditions (Hypo, NL, HL) per node.\n');
    fprintf(fid,'Pairwise Wilcoxon rank-sum: Hypo-NL, Hypo-HL, NL-HL per node.\n');
    fprintf(fid,'BH-FDR correction across nodes within each comparison family.\n');
    fprintf(fid,'Network nodes: %d; eligible tested nodes: %d (fixed loading inputs excluded).\n', ...
        height(BaselineStats),nEligible);
    fprintf(fid,'KW: q<%.3f for %d nodes.\n',BASE_STATS_ALPHA,nKW);
    fprintf(fid,'NL-HL: q<%.3f for %d nodes.\n',BASE_STATS_ALPHA,nNH);
    fprintf(fid,'These p/q values describe model-output variation among initializations, NOT biological replicates.\n');
    fclose(fid);
end

%% -------- PART C: Grouped bars for observed, with 95% CIs -------------
if ~isempty(observed)
  plot_grouped_bars_with_CI( ...
      observed_names, MEANS_ALL(observed,:), ...
      CI_LOW_ALL(observed,:), CI_HIGH_ALL(observed,:), ...
      {'Hypo','Normal','Hyper'}, [COL_HYPO; COL_NORM; COL_HYPER], ...
      'Observed nodes across exclusive load states', ...
      fullfile(fig_dir,'baseline_observed_grouped'), ...
      YLIM_ACTIVATION, true, 8);
end

%% -------- PART D: Category plots with grouped bars and 95% CIs ---------
for i = 1:length(group_categories)
  category_title = group_categories{i}{1};
  cat_nodes = map_custom_names(group_categories{i}{2}, NodeNames);
  if isempty(cat_nodes)
      warning('Category "%s" empty or unmatched; skipping.', category_title);
      continue;
  end
  plot_grouped_bars_with_CI( ...
      NodeLabels(cat_nodes), MEANS_ALL(cat_nodes,:), ...
      CI_LOW_ALL(cat_nodes,:), CI_HIGH_ALL(cat_nodes,:), ...
      {'Hypo','Normal','Hyper'}, [COL_HYPO; COL_NORM; COL_HYPER], ...
      category_title, fullfile(fig_dir, sprintf('group_%s', sanitize_name(category_title))), ...
      YLIM_ACTIVATION, false, 50);
end

%% -------- PART F: All-node comparisons with 95% bootstrap CIs ----------
% Split dense plots into pages so the error bars remain visible.
plot_grouped_bars_with_CI(NodeLabels, MEANS_ALL(:,[1,2]), ...
    CI_LOW_ALL(:,[1,2]), CI_HIGH_ALL(:,[1,2]), ...
    {'Hypo','Normal'}, [COL_HYPO; COL_NORM], ...
    'All nodes: Hypo vs Normal', fullfile(fig_dir,'allnodes_Hypo_vs_Normal'), ...
    YLIM_ACTIVATION, false, 40);

plot_grouped_bars_with_CI(NodeLabels, MEANS_ALL(:,[2,3]), ...
    CI_LOW_ALL(:,[2,3]), CI_HIGH_ALL(:,[2,3]), ...
    {'Normal','Hyper'}, [COL_NORM; COL_HYPER], ...
    'All nodes: Normal vs Hyper', fullfile(fig_dir,'allnodes_Normal_vs_Hyper'), ...
    YLIM_ACTIVATION, false, 40);

disp('Done. Corrected quarter-A4 figures and exports were written.');

%% ========================= HELPER FUNCTIONS ============================

function FS = get_pub_style()
    FS.title  = 11;
    FS.label  = 9.5;
    FS.tick   = 8;
    FS.legend = 8;
    FS.lw     = 1.0;
    FS.axlw   = 0.9;
end

function f = create_pub_figure(mode)
    switch lower(mode)
        case 'quarter'
            pos = [1 1 13.5 9.0];
        case 'quarterwide'
            pos = [1 1 14.0 9.0];
        case 'square'
            pos = [1 1 10.5 10.0];
        case 'forest'
            pos = [1 1 26.0 15.0];
        otherwise
            pos = [1 1 13.5 9.0];
    end

    f = figure( ...
        'Color','w', ...
        'Units','centimeters', ...
        'Position',pos, ...
        'PaperUnits','centimeters', ...
        'PaperPosition',[0 0 pos(3) pos(4)], ...
        'PaperSize',[pos(3) pos(4)], ...
        'Renderer','painters');
end

function apply_pub_axes(ax,FS,kind)
    if nargin<3, kind = 'default'; end

    ax.FontName = 'Arial';
    ax.FontSize = FS.tick;
    ax.FontWeight = 'bold';
    ax.LineWidth = FS.axlw;
    ax.Box = 'off';
    ax.TickDir = 'out';
    ax.TickLength = [0.010 0.010];
    ax.Layer = 'top';
    ax.XColor = [0.1 0.1 0.1];
    ax.YColor = [0.1 0.1 0.1];

    grid(ax,'on');
    ax.GridColor = [0.80 0.80 0.80];
    ax.GridAlpha = 0.35;

    switch lower(kind)
        case 'bar'
            ax.Position = [0.10 0.28 0.86 0.58];
        case 'bar_rot45'
            ax.Position = [0.10 0.34 0.86 0.52];
        case 'bar_rot60'
            ax.Position = [0.10 0.40 0.86 0.46];
        case 'forest'
            ax.Position = [0.08 0.37 0.90 0.55];
        case 'strip'
            ax.Position = [0.08 0.33 0.88 0.50];
        otherwise
            ax.Position = [0.12 0.14 0.82 0.78];
    end
end

function labels = make_readable_labels(labelsIn, maxChars, wrapLabels)
    labels = string(labelsIn);

    for ii = 1:numel(labels)
        s = char(labels(ii));

        if wrapLabels
            s = regexprep(s, '([\-_/])', '$1 ');
            words = strsplit(s, ' ');
            out = "";
            line = "";

            for k = 1:numel(words)
                candidate = strtrim(line + " " + string(words{k}));
                if strlength(candidate) > maxChars && strlength(line) > 0
                    out = out + line + newline;
                    line = string(words{k});
                else
                    line = candidate;
                end
            end
            s = char(out + line);
        end

        labels(ii) = string(strtrim(s));
    end
end

function idx = choose_tick_idx(n, maxLabels)
    if n <= maxLabels
        idx = 1:n;
    else
        idx = unique(round(linspace(1, n, maxLabels)));
    end
end

function lims = expand_limits(lims, frac)
    d = lims(2) - lims(1);
    if d <= 0, d = max(abs(lims(1)),1); end
    pad = frac * d;
    lims = [lims(1)-pad, lims(2)+pad];
end

function token = initWorkerBaseVars(Mact,Minh,NodeNames,NumOfNodes)
    assignin('base','Mact', Mact);
    assignin('base','Minh', Minh);
    assignin('base','NodeNames', NodeNames);
    assignin('base','NumOfNodes', NumOfNodes);
    token = true;
end

function idx = map_custom_names(custom_names, NodeNames)
    if isempty(custom_names), idx = []; return; end
    custom_names = custom_names(~cellfun(@isempty, custom_names));
    if isempty(custom_names), idx = []; return; end
    NN_norm = normalize_for_match(NodeNames);
    idx = [];
    missing = {};
    for k = 1:numel(custom_names)
        c = custom_names{k};
        c_fix = fix_common_synonyms(c);
        c_norm = normalize_for_match(c_fix);
        hit = find(strcmpi(c_norm, NN_norm), 1);
        if isempty(hit), hit = find(strcmpi(c_fix, NodeNames), 1); end
        if isempty(hit), missing{end+1} = c; else, idx(end+1) = hit; end %#ok<AGROW>
    end
    if ~isempty(missing)
        warning('Custom items not found: %s', strjoin(missing, ', '));
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

function RepLong = build_replicates_long(NodeNames, H, N, Hy)
    n = numel(NodeNames);
    R_H = size(H,1); R_N = size(N,1); R_Hy = size(Hy,1);
    Node_col = [repmat(NodeNames(:), R_H, 1); repmat(NodeNames(:), R_N, 1); repmat(NodeNames(:), R_Hy, 1)];
    Cond_col = [repmat("Hypo", R_H*n, 1); repmat("Normal", R_N*n, 1); repmat("Hyper", R_Hy*n, 1)];
    Repl_col = [repelem((1:R_H)', n); repelem((1:R_N)', n); repelem((1:R_Hy)', n)];
    Value_col = [H(:); N(:); Hy(:)];
    RepLong = table(string(Node_col), string(Cond_col), Repl_col, Value_col, ...
        'VariableNames', {'Node','Condition','Replicate','Value'});
end

function [SS_all, SS_mu, SS_sd] = run_replicates_clamped3(HypoVal,NLVal,HLVal,clamp_idx,NumOfNodes,tspan,num_iterations,solver_handle,ode_opts,BASE_VARS)
    SS_all = nan(num_iterations, NumOfNodes);
    poolActive = ~isempty(gcp('nocreate'));

    clamp_vals = [HypoVal; NLVal; HLVal];
    clamp_mask = false(NumOfNodes,1);
    clamp_mask(clamp_idx) = true;

    if poolActive
        parfor r = 1:num_iterations
            if ~isempty(BASE_VARS); BASE_VARS.Value; end
            opts  = ode_opts;
            Xinit = rand(NumOfNodes,1);
            Xinit(clamp_idx) = clamp_vals;

            [~, Xout] = solver_handle(@(t,x) ODESysFun_Clamped_safe3(t,x,clamp_idx,clamp_vals,clamp_mask), tspan, Xinit, opts);
            SS_all(r,:) = Xout(end,:);
        end
    else
        opts = ode_opts;
        for r = 1:num_iterations
            Xinit = rand(NumOfNodes,1);
            Xinit(clamp_idx) = clamp_vals;

            [~, Xout] = solver_handle(@(t,x) ODESysFun_Clamped_safe3(t,x,clamp_idx,clamp_vals,clamp_mask), tspan, Xinit, opts);
            SS_all(r,:) = Xout(end,:);
        end
    end

    SS_mu = mean(SS_all,1,'omitnan');
    SS_sd = std(SS_all,0,1,'omitnan');
end

function dxdt = ODESysFun_Clamped_safe3(~, x, clamp_idx, clamp_vals, clamp_mask)
    x = apply_clamp3(x, clamp_idx, clamp_vals);
    dxdt = ODESysFun(0, x);
    dxdt(clamp_mask) = 0;
end

function y = apply_clamp3(x, clamp_idx, clamp_vals)
    y = x;
    y(clamp_idx) = clamp_vals;
end

function ensure_dir(d)
  if exist(d,'dir')~=7
    [ok,msg] = mkdir(d);
    if ~ok, error('Failed to create directory: %s%s', d, msg); end
  end
end

function savefig_safe(h, filenameNoExt)
    try
        ensure_dir(fileparts(filenameNoExt));

        exportgraphics(h, [filenameNoExt '.png'], ...
            'Resolution', 600, ...
            'BackgroundColor', 'white');

        try
            exportgraphics(h, [filenameNoExt '.pdf'], ...
                'ContentType', 'vector', ...
                'BackgroundColor', 'white');
        catch
        end

    catch ME
        warning('exportgraphics failed (%s). Falling back to print.', ME.message);
        try
            print(h, [filenameNoExt '.png'], '-dpng', '-r600');
            try
                print(h, [filenameNoExt '.pdf'], '-dpdf', '-painters');
            catch
            end
        catch ME2
            warning('Could not save figure to "%s": %s', filenameNoExt, ME2.message);
        end
    end
end

function s = sanitize_name(s)
  s = regexprep(s,'[^A-Za-z0-9_]','_');
end

function y = safe_divide(a,b)
  y = nan(size(a));
  mask = abs(b) > eps;
  y(mask) = a(mask)./b(mask);
end

function plot_grouped_bars_with_CI(xlabels,Y,CIlo,CIhi,condOrder,condColors,titleStr,savepath,ylims,isObserved,maxNodesPerFig)
    % Mean-activation bars with asymmetric 95% percentile-bootstrap CIs.
    % CI values are computed once per node/condition from replicate outputs.
    if nargin < 9 || isempty(ylims), ylims = []; end
    if nargin < 10 || isempty(isObserved), isObserved = false; end
    if nargin < 11 || isempty(maxNodesPerFig)
        if isObserved, maxNodesPerFig = 8; else, maxNodesPerFig = 40; end
    end

    assert(isequal(size(Y),size(CIlo),size(CIhi)), ...
        'Y, CIlo and CIhi must have matching dimensions');
    assert(size(Y,2) == numel(condOrder), ...
        'Columns of Y must match condOrder');

    N = size(Y,1);
    nPages = ceil(N / maxNodesPerFig);
    for pg = 1:nPages
        idx1 = (pg-1)*maxNodesPerFig + 1;
        idx2 = min(pg*maxNodesPerFig,N);
        idx = idx1:idx2;
        xlabels_pg = xlabels(idx);
        Ypg  = Y(idx,:);
        Lpg  = CIlo(idx,:);
        Hpg  = CIhi(idx,:);
        Npg  = numel(idx);
        FS = get_pub_style();

        if Npg <= 10, f = create_pub_figure('quarter');
        else, f = create_pub_figure('quarterwide'); end
        ax = axes('Parent',f); hold(ax,'on');
        x = 1:Npg;
        bh = bar(ax,x,Ypg,'grouped','BarWidth',0.86);
        for k = 1:numel(bh)
            bh(k).FaceColor = condColors(k,:);
            bh(k).EdgeColor = 'none';
            bh(k).LineWidth = 0.5;
        end

        % Use the individual bar centers (not group centers) for each CI.
        % Mask missing intervals without introducing invented error sizes.
        for k = 1:numel(bh)
            valid = isfinite(Ypg(:,k)) & isfinite(Lpg(:,k)) & isfinite(Hpg(:,k));
            if ~any(valid), continue; end
            lowerLength = max(0,Ypg(valid,k)-Lpg(valid,k));
            upperLength = max(0,Hpg(valid,k)-Ypg(valid,k));
            errorbar(ax,bh(k).XEndPoints(valid),Ypg(valid,k), ...
                lowerLength,upperLength, ...
                'LineStyle','none','Color',[0.10 0.10 0.10], ...
                'LineWidth',0.75,'CapSize',2.5,'HandleVisibility','off');
        end

        ylabel(ax,'Mean activation (a.u.)', ...
            'FontSize',FS.label,'FontWeight','bold');
        if nPages > 1
            title_here = sprintf('%s (%d/%d)',titleStr,pg,nPages);
            savepath_here = sprintf('%s_part_%02d',savepath,pg);
        else
            title_here = titleStr;
            savepath_here = savepath;
        end
        title(ax,title_here,'Interpreter','none', ...
            'FontSize',FS.title,'FontWeight','bold');

        if ~isempty(ylims), ylim(ax,ylims); end
        xlim(ax,[0.35,Npg+0.65]);
        labels = make_readable_labels(xlabels_pg,10,true);
        if Npg <= 8
            rot=0; kind='bar';
        elseif Npg <= 14
            rot=45; kind='bar_rot45';
        else
            rot=60; kind='bar_rot60';
        end
        set(ax,'XTick',1:Npg,'XTickLabel',labels,'XTickLabelRotation',rot);
        ax.TickLabelInterpreter = 'none';
        legend(ax,bh,condOrder,'Location','northoutside', ...
            'Orientation','horizontal','Box','off', ...
            'FontSize',FS.legend,'FontWeight','bold');
        apply_pub_axes(ax,FS,kind);
        savefig_safe(f,savepath_here);
        close(f);
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

function [ciLow, ciHigh] = bootstrap_diff_CI(repA, repB, nBoot, alpha)
    [~, nNodes] = size(repA);
    ciLow  = nan(1,nNodes);
    ciHigh = nan(1,nNodes);
    for j = 1:nNodes
        vA = repA(:,j); vA = vA(~isnan(vA));
        vB = repB(:,j); vB = vB(~isnan(vB));
        if isempty(vA) || isempty(vB)
            ciLow(j) = NaN; ciHigh(j) = NaN; continue;
        end
        diffs = nan(nBoot,1);
        nA = numel(vA); nB = numel(vB);
        for b = 1:nBoot
            idxA = randi(nA, nA, 1);
            idxB = randi(nB, nB, 1);
            diffs(b) = mean(vA(idxA)) - mean(vB(idxB));
        end
        diffs = sort(diffs);
        lo = max(1, floor((alpha/2) * nBoot));
        hi = min(nBoot, ceil((1-alpha/2) * nBoot));
        ciLow(j)  = diffs(lo);
        ciHigh(j) = diffs(hi);
    end
end

function [ciLo,ciHi] = bootstrap_mean_CI(replicates,nBoot,alpha)
    % Percentile bootstrap CI for condition-specific mean node activation.
    % One column per node; replicate observations are sampled with replacement.
    [~,nNodes] = size(replicates);
    ciLo = nan(1,nNodes); ciHi = nan(1,nNodes);
    for j = 1:nNodes
        x = replicates(:,j); x = x(isfinite(x));
        if isempty(x), continue; end
        if numel(x) == 1
            ciLo(j)=x; ciHi(j)=x; continue;
        end
        n = numel(x);
        % Vectorized resampling: n replicates x nBoot bootstrap draws.
        draws = randi(n,n,nBoot);
        bootMeans = mean(x(draws),1);
        bounds = prctile(bootMeans,100*[alpha/2,1-alpha/2]);
        ciLo(j)=bounds(1); ciHi(j)=bounds(2);
    end
end

function Stats = baseline_nonparametric_stats(HypoRep,NLRep,HLRep,NodeNames,clampIdx)
    % Distributional comparisons across independent initializations.
    % Mann-Whitney-Wilcoxon rank-sum: unpaired, no normality assumption.
    % KW: three-condition omnibus test. p/q are NOT experimental p-values.
    % Fixed input nodes are not tested; their values are imposed by design.
    nNodes = numel(NodeNames);
    pKW = nan(nNodes,1);
    pHN = nan(nNodes,1);
    pHH = nan(nNodes,1);
    pNH = nan(nNodes,1);
    med = nan(nNodes,3);
    q1 = nan(nNodes,3);
    q3 = nan(nNodes,3);
    sampleCount = zeros(nNodes,3);

    for j=1:nNodes
        h = HypoRep(:,j); h = h(isfinite(h));
        n = NLRep(:,j); n = n(isfinite(n));
        l = HLRep(:,j); l = l(isfinite(l));
        vals = {h,n,l};
        for c=1:3
            v=vals{c};
            sampleCount(j,c)=numel(v);
            if isempty(v), continue; end
            med(j,c)=median(v);
            qt=prctile(v,[25,75]);
            q1(j,c)=qt(1); q3(j,c)=qt(2);
        end
        if ismember(j,clampIdx) || isempty(h)||isempty(n)||isempty(l)
            continue;
        end
        if all([h;n;l]==h(1))
            pKW(j)=1;
        else
            pKW(j)=kruskalwallis([h;n;l], ...
                [ones(numel(h),1);2*ones(numel(n),1);3*ones(numel(l),1)],'off');
        end
        pHN(j)=ranksum_safe(h,n);
        pHH(j)=ranksum_safe(h,l);
        pNH(j)=ranksum_safe(n,l);
    end
    Stats=table(string(NodeNames(:)), ...
        sampleCount(:,1),sampleCount(:,2),sampleCount(:,3), ...
        med(:,1),q1(:,1),q3(:,1), ...
        med(:,2),q1(:,2),q3(:,2), ...
        med(:,3),q1(:,3),q3(:,3), ...
        pKW,bh_fdr(pKW),pHN,bh_fdr(pHN), ...
        pHH,bh_fdr(pHH),pNH,bh_fdr(pNH), ...
        'VariableNames',{'Node','N_Hypo','N_NL','N_HL', ...
        'Hypo_Median','Hypo_Q1','Hypo_Q3', ...
        'NL_Median','NL_Q1','NL_Q3', ...
        'HL_Median','HL_Q1','HL_Q3', ...
        'p_KW_3group','q_KW_3group','p_Hypo_NL','q_Hypo_NL', ...
        'p_Hypo_HL','q_Hypo_HL','p_NL_HL','q_NL_HL'});
end

function p=ranksum_safe(a,b)
    if isempty(a)||isempty(b), p=NaN; return; end
    if all(a==a(1)) && all(b==b(1)) && a(1)==b(1)
        p=1; return;
    end
    p=ranksum(a,b);
end

function q=bh_fdr(p)
    % Benjamini-Hochberg adjusted p-values (q-values).
    q=nan(size(p));
    valid=find(isfinite(p) & p>=0 & p<=1);
    if isempty(valid), return; end
    [ps,ord]=sort(p(valid));
    m=numel(ps);
    raw=ps(:).*m./(1:m)';
    qs=flipud(cummin(flipud(raw)));
    qs=min(1,qs);
    q(valid(ord))=qs;
end


function plot_falsification_AB_v12(F, B, out_dir)
% Integrated unchanged v12 visual design, with current-run tables as inputs.
    must_have(F,{'Node','Class','Mean_Normal','Mean_Hyper','Delta_NminusH','CI_low','CI_high'});
    must_have(B,{'Node','NL_Mean','NL_CI_Low','NL_CI_High', ...
                 'HL_Mean','HL_CI_Low','HL_CI_High'});
    classes = lower(strtrim(string(F.Class)));
    F = F(classes=="anabolic" | classes=="catabolic",:);
    classes = lower(strtrim(string(F.Class)));
    nodes = strtrim(string(F.Node));
    nl = nums(F.Mean_Normal); hl = nums(F.Mean_Hyper);
    delta = nums(F.Delta_NminusH);
    ciLo = nums(F.CI_low); ciHi = nums(F.CI_high);
    n = numel(nodes);
    assert(n == 39,'The falsification output must have exactly 39 rows; found %d.',n);
    assert(numel(unique(nodes)) == n,'Benchmark names are not unique.');
    assert(sum(classes=="anabolic")==6 && sum(classes=="catabolic")==33, ...
        'Benchmark class totals must be exactly 6 and 33.');
    [delta,ord] = sort(delta,'descend');
    nodes=nodes(ord); classes=classes(ord);
    nl=nl(ord); hl=hl(ord); ciLo=ciLo(ord); ciHi=ciHi(ord);
    [matched,idx] = ismember(nodes,strtrim(string(B.Node)));
    assert(all(matched),'Some benchmark nodes are missing from baseline CIs.');
    nlLo=nums(B.NL_CI_Low(idx)); nlHi=nums(B.NL_CI_High(idx));
    hlLo=nums(B.HL_CI_Low(idx)); hlHi=nums(B.HL_CI_High(idx));
    bNL=nums(B.NL_Mean(idx)); bHL=nums(B.HL_Mean(idx));
    assert(all(isfinite([nl;hl;delta;ciLo;ciHi;nlLo;nlHi;hlLo;hlHi])), ...
        'One or more benchmark means or confidence intervals are not finite.');
    assert(all(ciLo<=ciHi & nlLo<=nlHi & hlLo<=hlHi),'Reversed CI bounds.');
    assert(all(abs(delta-(nl-hl))<1e-5),'Falsification differences mismatch means.');
    assert(all(abs(nl-bNL)<1e-5 & abs(hl-bHL)<1e-5), ...
        'Falsification and baseline tables do not match.');

    % One coherent color scheme.
    GREEN = [0.00 0.69 0.14]; RED = [0.92 0.075 0.095];
    DARKGREEN = [0.00 0.34 0.06]; DARKRED = [0.58 0.00 0.04];
    GRAY=[0.78 0.80 0.82]; DARK=[0.14 0.14 0.14];
    EDGE=[0.45 0.45 0.45];

    out = fullfile(out_dir,'Supplementary_Falsification_AB');
    if exist(out,'dir')~=7, mkdir(out); end

    %% FIGURE A: ONE AXIS, TWO PARALLEL BAR COLUMNS (NL and HL)
    % This avoids the two-axes geometry that repeatedly hid the NL rows.
    figA = figure('Visible','off','Color','white','Units','pixels', ...
        'Position',[50 50 1200 900],'Renderer','painters');
    axA = axes('Parent',figA,'Units','normalized','Position',[.17 .11 .81 .85]);
    hold(axA,'on');
    xHL = 1.24;  % identical 0-1 scale, with a modest gutter
    ymax = n+0.5;
    xMax = xHL+1.04;

    % Subtle full-row stripes, then fine grid lines behind bars.
    for i=1:2:n
        patch(axA,[0 xMax xMax 0],[i-.5 i-.5 i+.5 i+.5], ...
            [0.979 0.979 0.979],'EdgeColor','none','HandleVisibility','off');
    end
    for t=0:.2:1
        for offset=[0 xHL]
            line(axA,[offset+t offset+t],[.5 ymax], ...
                'Color',[.90 .90 .90],'LineWidth',.42,'HandleVisibility','off');
        end
    end
    for i=1:n
        bar_at(axA,0,nl(i),i,GREEN);
        bar_at(axA,xHL,hl(i),i,RED);
        ci_at(axA,nlLo(i),nlHi(i),i,DARKGREEN);
        ci_at(axA,xHL+hlLo(i),xHL+hlHi(i),i,DARKRED);
    end

    % Place headings INSIDE the axis; no title or overlapping annotations.
    text(axA,.50,-.43,'NL','HorizontalAlignment','center', ...
        'FontName','Arial','FontWeight','bold','FontSize',15,'Color',DARK);
    text(axA,xHL+.50,-.43,'HL','HorizontalAlignment','center', ...
        'FontName','Arial','FontWeight','bold','FontSize',15,'Color',DARK);

    xx=[0:.2:1, xHL+(0:.2:1)];
    lab=repmat({'0','0.2','0.4','0.6','0.8','1'},1,2);
    set(axA,'YDir','reverse','YLim',[-.90 ymax],'XLim',[-.02 xMax], ...
        'YTick',1:n,'YTickLabel',cellstr(nodes), ...
        'XTick',xx,'XTickLabel',lab,'TickLabelInterpreter','none', ...
        'FontName','Arial','FontWeight','bold','FontSize',10.5, ...
        'TickDir','out','TickLength',[.003 .003], ...
        'Box','off','Layer','top','XColor',DARK,'YColor',DARK);
    % Enlarged and BOLD x-axis numbers; keep benchmark labels compact.
    axA.XAxis.FontSize=15;
    axA.XAxis.FontWeight='bold';
    axA.YAxis.FontSize=10.5;
    axA.YAxis.FontWeight='bold';
    xlabel(axA,'Mean steady-state activation (0–1)', ...
        'FontName','Arial','FontSize',17,'FontWeight','bold','Color',DARK);
    drawnow;
    assert(numel(axA.YTick)==n,'Figure A does not display all benchmark row ticks.');
    export_all(figA,axA,fullfile(out,'Figure_A_NL_HL_Parallel_v12'));
    close(figA);

    %% FIGURE B: ONE AXIS, 39 BENCHMARK ROWS AND 2 DIAGONAL COLUMNS
    expDir=ones(n,1); expDir(classes=="catabolic")=-1;
    simDir=zeros(n,1); simDir(ciLo>0)=1; simDir(ciHi<0)=-1;
    support=sum(simDir==expDir);
    opposite=sum(simDir==-expDir);
    uncertain=sum(simDir==0);

    figB = figure('Visible','off','Color','white','Units','pixels', ...
        'Position',[50 50 850 1120],'Renderer','painters');
    axB = axes('Parent',figB,'Units','normalized','Position',[.23 .045 .75 .93]);
    hold(axB,'on');
    bounds=[.56 1.51; 1.69 2.64];  % white gutter between NL and HL
    centers=mean(bounds,2);

    
for i=1:n
    yt=i-.5; yb=i+.5;
    for col=1:2
        code=3-col;
        simC=color_for(simDir(i),code,GREEN,RED,GRAY);
        expC=color_for(expDir(i),code,GREEN,RED,GRAY);
        xl=bounds(col,1); xr=bounds(col,2);

        patch(axB,[xl xr xl],[yt yt yb],simC, ...
            'EdgeColor','none','HandleVisibility','off');
        patch(axB,[xr xr xl],[yt yb yb],expC, ...
            'EdgeColor','none','HandleVisibility','off');

        % Diagonal dividing line: 2 times thicker
        line(axB,[xr xl],[yt yb], ...
            'Color',EDGE,'LineWidth',0.76, ...
            'HandleVisibility','off');
    end
end

for c=1:2
    xl=bounds(c,1); xr=bounds(c,2);

    % Vertical rectangle borders: 3 times thicker
    line(axB,[xl xl],[.5 n+.5], ...
        'Color',EDGE,'LineWidth',1.65);

    line(axB,[xr xr],[.5 n+.5], ...
        'Color',EDGE,'LineWidth',1.65);

    % Horizontal rectangle borders: 3 times thicker
    for j=0:n
        line(axB,[xl xr],[j+.5 j+.5], ...
            'Color',EDGE,'LineWidth',0.90);
    end
end


    % Legend lives within this very same exported axis.
    % Scale positions to keep the summary visible in the 39-row panel.
    keyScale = n/49;
    draw_key(axB,3.04,12*keyScale,GREEN,'Higher under loading',DARK);
    draw_key(axB,3.04,17*keyScale,RED,'Lower under loading',DARK);
    draw_key(axB,3.04,22*keyScale,GRAY,'Simulated CI includes zero',DARK);
    text(axB,3.04,28*keyScale,'Upper-left: simulation', ...
        'FontName','Arial','FontSize',9.0,'Color',DARK);
    text(axB,3.04,31*keyScale,'Lower-right: experiment', ...
        'FontName','Arial','FontSize',9.0,'Color',DARK);
    text(axB,3.04,40*keyScale,sprintf('Supported: %d',support), ...
        'FontName','Arial','FontSize',9.6,'FontWeight','bold','Color',DARK);
    text(axB,3.04,43*keyScale,sprintf('Contradicted: %d',opposite), ...
        'FontName','Arial','FontSize',9.6,'FontWeight','bold','Color',DARK);
    text(axB,3.04,46*keyScale,sprintf('Inconclusive: %d',uncertain), ...
        'FontName','Arial','FontSize',9.6,'FontWeight','bold','Color',DARK);

    set(axB,'YDir','reverse','YLim',[.5 n+.5], ...
        'XLim',[.45 5.8], ...
        'YTick',1:n,'YTickLabel',cellstr(nodes), ...
        'XTick',centers,'XTickLabel',{'NL','HL'}, ...
        'XAxisLocation','top','TickLabelInterpreter','none', ...
        'FontName','Arial','FontSize',9.6,'FontWeight','bold', ...
        'TickLength',[0 0],'Box','off','Layer','top', ...
        'XColor',DARK,'YColor',DARK);
    xlabel(axB,''); ylabel(axB,'');
    drawnow;
    assert(numel(axB.YTick)==n,'Figure B does not display all benchmark row ticks.');
    export_all(figB,axB,fullfile(out,'Figure_B_Vertical_Signs_v12'));
    close(figB);

    fprintf('\nSaved two SEPARATE figures with all %d benchmark nodes.\n',n);
    fprintf('A: Figure_A_NL_HL_Parallel_v12  (Preview.png, Print600dpi.png, pdf, fig)\n');
    fprintf('B: Figure_B_Vertical_Signs_v12 (Preview.png, Print600dpi.png, pdf, fig)\n');
    fprintf('Location: %s\n',out);
    fprintf('Supported %d | Contradicted %d | Inconclusive %d\n', ...
        support,opposite,uncertain);
    fprintf('Upload the *_Preview.png files when checking the layout.\n');
end

function bar_at(ax,offset,v,y,color)
    w=.33;
    patch(ax,[offset offset+v offset+v offset], ...
        [y-w y-w y+w y+w],color,'EdgeColor','none','HandleVisibility','off');
end
function ci_at(ax,lo,hi,y,color)
    line(ax,[lo hi],[y y],'Color',color,'LineWidth',1.1,'HandleVisibility','off');
    line(ax,[lo lo],[y-.17 y+.17],'Color',color,'LineWidth',.95,'HandleVisibility','off');
    line(ax,[hi hi],[y-.17 y+.17],'Color',color,'LineWidth',.95,'HandleVisibility','off');
end
function rgb=color_for(d,code,g,r,neutral)
    if d==0
        rgb=neutral;
    elseif (d==1 && code==2) || (d==-1 && code==1)
        rgb=g;
    else
        rgb=r;
    end
end
function symbol_at(ax,x,y,d,code)
    positive=(d==1 && code==2) || (d==-1 && code==1);
    if positive
        sign='+';
        col=[0.05 0.12 0.05];
    else
        sign='-';
        col=[1 1 1];
    end
    text(ax,x,y,sign,'FontName','Arial','FontSize',9.5, ...
        'FontWeight','bold','HorizontalAlignment','center', ...
        'VerticalAlignment','middle','Color',col, ...
        'Interpreter','none','Clipping','on','HandleVisibility','off');
end
function draw_key(ax,x,y,col,label,dark)
    patch(ax,[x x+.22 x+.22 x],[y-.55 y-.55 y+.45 y+.45], ...
        col,'EdgeColor',[.50 .50 .50],'LineWidth',.5,'HandleVisibility','off');
    text(ax,x+.31,y,label,'FontName','Arial','FontSize',8.7, ...
        'Interpreter','none','VerticalAlignment','middle','Color',dark);
end
function export_all(fig,ax,base)
    drawnow;
    % Preview is intentionally < 2048 px, avoiding cropping in viewers.
    exportgraphics(ax,[base '_Preview.png'], ...
        'Resolution',140,'BackgroundColor','white');
    % This full-resolution publication image can be much taller than a preview.
    exportgraphics(ax,[base '_Print600dpi.png'], ...
        'Resolution',600,'BackgroundColor','white');
    exportgraphics(ax,[base '.pdf'], ...
        'ContentType','vector','BackgroundColor','white');
    savefig(fig,[base '.fig']);
end
function must_have(T,names)
    absent=setdiff(names,T.Properties.VariableNames);
    if ~isempty(absent)
        error('Saved file is missing columns: %s',strjoin(absent,', '));
    end
end

function v=nums(x)
    if isnumeric(x) || islogical(x)
        v=double(x(:));
    else
        v=str2double(string(x(:)));
    end
end


function validate_benchmark_registry(anab, cat, NodeNames)
    expected = [anab(:);cat(:)];
    assert(numel(anab)==6 && numel(cat)==33 && numel(expected)==39, ...
        'Revised registry must contain 6 NL-higher and 33 HL-higher markers.');
    normalized = normalize_for_match(NodeNames);
    errors = {};
    matched = zeros(numel(expected),1);
    for k=1:numel(expected)
        tag = expected{k};
        key = normalize_for_match(fix_common_synonyms(tag));
        hits = find(strcmpi(key{1},normalized));
        if numel(hits)~=1
            errors{end+1} = sprintf('%s (%d matching network nodes)',tag,numel(hits)); %#ok<AGROW>
        else
            matched(k)=hits;
        end
    end
    if ~isempty(errors)
        error('Revised benchmark registry cannot be mapped unambiguously: %s. Update MT_PRIMARY4_1.xlsx or the benchmark names before running.', ...
            strjoin(errors,', '));
    end
    assert(numel(unique(matched))==39,'Benchmark registry maps more than one marker to the same network node.');
end
