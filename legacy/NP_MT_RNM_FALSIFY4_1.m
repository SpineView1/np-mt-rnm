%% MTN_AllInOne_EXCLUSIVE_BASELINES_ONLY_v4_FALSIFY_3INPUT.m
% Corrected version:
% 1) Uses the same exclusive 3-input chronic loading scheme as the main baseline code
% 2) Grouped bar plots sized for quarter-A4 PowerPoint composition
% 3) Large node sets automatically split across multiple figures
% 4) All nodes plotted correctly in all-node grouped bar plots
% 5) Robust clamping logic works even when a clamped value is 0

clc; close all; clear;

%% ------------------------- USER SETTINGS -------------------------------
netfile   = 'MT_PRIMARY4_1.xlsx';
out_xlsx  = 'Mechanotransduction_Results.xlsx';
fig_dir   = fullfile(pwd,'figures_FALSIF4_1');  ensure_dir(fig_dir);
out_dir   = fullfile(pwd,'outputs_FALSIF4_1');  ensure_dir(out_dir);

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
         'CAT','GPX1','SOD1','SOD2','HO-1','PHD2','VHL','Rheb','NutD','MitD'}};
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
    {'Rho GTPases, cytoskeletal & Hippo regulators', ...
        {'RhoA-M','RhoA-E','RAC1-M','RAC1-E','CDC42','ROCK-M','ROCK-E','PAK1','PKN1','FAK-M','FAK-E','MST1/2','LATS1/2'}};
};

tspan_baseline = [0 100];
rng(1);

% Three-input exclusive regime definitions
Hypo_hypo = 0.20; NL_hypo = 0.01; HL_hypo = 0.01;
Hypo_norm = 0.01; NL_norm = 0.80; HL_norm = 0.01;
Hypo_hype = 0.01; NL_hype = 0.10; HL_hype = 0.80;

YLIM_ACTIVATION = [0 1];

USE_PAR      = true;
PAR_WORKERS  = [];

ANABOLIC_LIST  = {...
    'SOX9','ACAN','COL2A1','TIMP3','Bcl2','AQP1','AQP5',...
    'HO-1','GPX1','SOD2','CAT','HSP27','HSP70','NRF2',...
    'HIF-1α','SIRT1','AMPK'};
CATABOLIC_LIST = {...
    'COL1A1','COL10A1','IL1β','IL6','IL8','TNF','MMP3','MMP13',...
    'ADAMTS4/5','CCL2','CXCL1','CXCL3',...
    'COX-2','PGE2','BAX','CASP3','CASP9',...
    'Wnt3a','β-catenin','NF-κB','MRTF-A','ROS',...
    'YAP/TAZ','TonEBP','p38','JNK','HSF1','TGFβ'};
FALS_TOL   = 0.02;
FALS_ALPHA = 0.05;
FALS_NBOOT = 10000;

%% ------------------------- LOAD NETWORK --------------------------------
[Mact, Minh, NodeNames, NumOfNodes] = CreateMatrices_new(netfile);
NodeLabels = ensure_unique_labels(NodeNames);

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
        Pass_NH(i) = (ciLow_NH(i) > FALS_TOL);
    elseif Classes(i)=="catabolic"
        Pass_NH(i) = (ciHigh_NH(i) < -FALS_TOL);
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
    FALS_TOL, FALS_ALPHA, FALS_NBOOT, ...
    'VariableNames',{'TotalRules','Passed','Failed', ...
    'Anabolic_N','Anabolic_Passed','Anabolic_PassRate', ...
    'Catabolic_N','Catabolic_Passed','Catabolic_PassRate', ...
    'Tolerance','Alpha','NBoot'});

writetable(FalsTbl_NH_Export, fullfile(out_dir,'Falsification_Normal_vs_Hyper.csv'));
try, writetable(FalsTbl_NH_Export, out_xlsx, 'Sheet','Falsification_NvH'); catch, end
writetable(FalsSummary_NH, fullfile(out_dir,'Falsification_Summary_NvH.csv'));
try, writetable(FalsSummary_NH, out_xlsx, 'Sheet','Falsification_Summary_NvH'); catch, end

DiffNL = MeanN - MeanL;
[ciLow_NL, ciHigh_NL] = bootstrap_diff_CI(SS_norm_all, SS_hypo_all, FALS_NBOOT, FALS_ALPHA);

Pass_NL = false(NumOfNodes,1);
for i = 1:NumOfNodes
    if Classes(i)=="anabolic"
        Pass_NL(i) = (ciLow_NL(i) > FALS_TOL);
    elseif Classes(i)=="catabolic"
        Pass_NL(i) = (ciHigh_NL(i) < -FALS_TOL);
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
    FALS_TOL, FALS_ALPHA, FALS_NBOOT, ...
    'VariableNames',{'TotalRules','Passed','Failed', ...
    'Anabolic_N','Anabolic_Passed','Anabolic_PassRate', ...
    'Catabolic_N','Catabolic_Passed','Catabolic_PassRate', ...
    'Tolerance','Alpha','NBoot'});

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

%% Forest plot: Normal vs Hyper
testedIdx_NH = find(maskTest);
[~, order_NH] = sort(DiffNH(testedIdx_NH), 'descend');
ordIdx_NH = testedIdx_NH(order_NH);

if ~isempty(ordIdx_NH)
    f = create_pub_figure('forest');
    ax = axes('Parent',f); hold(ax,'on');
    y = 1:numel(ordIdx_NH);

    for k = 1:numel(ordIdx_NH)
        i = ordIdx_NH(k);
        line(ax,[ciLow_NH(i) ciHigh_NH(i)], [y(k) y(k)], 'LineWidth', 1.8, 'Color', [0.55 0.55 0.55]);
    end
    for k = 1:numel(ordIdx_NH)
        i = ordIdx_NH(k);
        col = COL_ANAB;
        if Classes(i)=="catabolic", col = COL_CATA; end
        mk = 'o';
        if ~Pass_NH(i), mk = 's'; end
        plot(ax, DiffNH(i), y(k), mk, 'MarkerFaceColor', col, 'MarkerEdgeColor','k', ...
            'MarkerSize', 6.5, 'LineWidth',0.9);
    end

    xline(ax,0,'k-','LineWidth',1.1);
    xline(ax,FALS_TOL,'--','Color',COL_ANAB,'LineWidth',1.0);
    xline(ax,-FALS_TOL,'--','Color',COL_CATA,'LineWidth',1.0);

    labels = make_readable_labels(NodeNames(ordIdx_NH), 28, true);
    set(ax,'YTick', y, 'YTickLabel', labels, 'YDir','reverse');
    xlabel(ax,'Mean difference (Normal - Hyper)','FontSize',FS.label,'FontWeight','bold');
    ylabel(ax,'Tested nodes','FontSize',FS.label,'FontWeight','bold');
    title(ax,sprintf('Falsification forest plot (Normal - Hyper, %.0f%% CI)',(1-FALS_ALPHA)*100), ...
        'FontSize',FS.title,'FontWeight','bold');
    apply_pub_axes(ax,FS,'forest');
    ax.TickLabelInterpreter = 'none';
    savefig_safe(f, fullfile(out_dir,'falsification_forest_NormalMinusHyper'));
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

%% Replicate distributions: top nodes N vs H
aOrd_NH = ordIdx_NH(Classes(ordIdx_NH)=="anabolic");
cOrd_NH = ordIdx_NH(Classes(ordIdx_NH)=="catabolic");
[~, aSort_NH] = sort(abs(DiffNH(aOrd_NH)), 'descend');
aPick_NH = aOrd_NH(aSort_NH(1:min(3,end)));
[~, cSort_NH] = sort(abs(DiffNH(cOrd_NH)), 'descend');
cPick_NH = cOrd_NH(cSort_NH(1:min(3,end)));
sel_NH = [aPick_NH; cPick_NH];

for s = 1:numel(sel_NH)
    i = sel_NH(s);
    f = create_pub_figure('quarter');
    ax = axes('Parent',f); hold(ax,'on');
    boxplot(ax,[SS_norm_all(:,i), SS_hype_all(:,i)], 'Labels', {'Normal','Hyper'}, ...
        'Symbol','', 'Widths',0.55, 'Whisker',1.5);
    set(findobj(ax,'Tag','Box'),'LineWidth',1.2);
    set(findobj(ax,'Tag','Median'),'LineWidth',1.3,'Color','k');
    title(ax,sprintf('%s | %s | Pass=%d | \\Delta=%.3f', NodeNames{i}, Classes(i), Pass_NH(i), DiffNH(i)), ...
        'Interpreter','none','FontSize',FS.title-1,'FontWeight','bold');
    ylabel(ax,'Steady-state activation','FontSize',FS.label,'FontWeight','bold');
    if ~isempty(YLIM_ACTIVATION), ylim(ax,YLIM_ACTIVATION); end
    apply_pub_axes(ax,FS,'default');
    savefig_safe(f, fullfile(out_dir,sprintf('falsification_replicates_NvH_%s', sanitize_name(NodeNames{i}))));
    close(f);
end

%% Summary pass-rate bars
f = create_pub_figure('quarter');
ax = axes('Parent',f);
cats = categorical({'Anabolic','Catabolic','Overall'});
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
legend(ax,{'Anabolic (pass)','Anabolic (fail)','Catabolic (pass)','Catabolic (fail)','Neutral'}, ...
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

        f = create_pub_figure('quarterwide');
        ax = axes('Parent',f); hold(ax,'on');

        % Data for this page
        vals = DiffNH(ordIdx(idxp));

        bh = bar(ax, vals, 0.92, 'FaceColor','flat','EdgeColor','none');

        % Colors
        C = zeros(numel(idxp),3);
        C(Classes(ordIdx(idxp))=="anabolic",:) = repmat(COL_ANAB, sum(Classes(ordIdx(idxp))=="anabolic"),1);
        C(Classes(ordIdx(idxp))=="catabolic",:) = repmat(COL_CATA, sum(Classes(ordIdx(idxp))=="catabolic"),1);
        bh.CData = C;

        % Reference lines
        yline(ax,0,'k-','LineWidth',1.0);
        yline(ax,FALS_TOL,'--','Color',COL_ANAB,'LineWidth',0.9);
        yline(ax,-FALS_TOL,'--','Color',COL_CATA,'LineWidth',0.9);

        % Labels (ALL labels per page)
        labels = make_readable_labels(NodeNames(ordIdx(idxp)), 14, false);

        set(ax, ...
            'XTick',1:numel(idxp), ...
            'XTickLabel',labels, ...
            'XTickLabelRotation',60);

        ylabel(ax,'\Delta (Normal - Hyper)','FontSize',FS.label,'FontWeight','bold');

        if nPages > 1
            title_here = sprintf('Directional effects (Normal - Hyper) (%d/%d)', pg, nPages);
            savepath_here = fullfile(out_dir, sprintf('directional_effects_bar_NvH_part_%02d', pg));
        else
            title_here = sprintf('Directional effects by node (Normal - Hyper, %.0f%% CI)',(1-FALS_ALPHA)*100);
            savepath_here = fullfile(out_dir,'directional_effects_bar_NvH');
        end

        title(ax, title_here, 'FontSize',FS.title,'FontWeight','bold');

        apply_pub_axes(ax,FS,'bar_rot60');
        ax.TickLabelInterpreter = 'none';

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
legend(ax,{'Anabolic','Catabolic','Neutral'}, 'Location','southoutside', ...
    'Orientation','horizontal','Box','off','FontSize',FS.legend,'FontWeight','bold');
apply_pub_axes(ax,FS,'default');
savefig_safe(f, fullfile(out_dir,'stability_map_CV_NvH_enhanced'));
close(f);

%% -------- PART C: Grouped bars for observed ----------------------------
if ~isempty(observed)
  dataObs = [SS_hypo_mu(observed)', SS_norm_mu(observed)', SS_hype_mu(observed)'];
  plot_grouped_bars_no_errors(...
      observed_names, dataObs, {'Hypo','Normal','Hyper'}, ...
      [COL_HYPO; COL_NORM; COL_HYPER], ...
      'Observed nodes across exclusive load states', ...
      fullfile(fig_dir,'baseline_observed_grouped'), ...
      YLIM_ACTIVATION, true, 8);
end

%% --------- PART D: Category plots with grouped bars --------------------
X_first  = SS_hypo_mu(:);
X_second = SS_norm_mu(:);
X_third  = SS_hype_mu(:);

for i = 1:length(group_categories)
  category_title = group_categories{i}{1};
  cat_nodes = map_custom_names(group_categories{i}{2}, NodeNames);
  if isempty(cat_nodes)
      warning('Category "%s" empty or unmatched; skipping.', category_title);
      continue;
  end
  Y = [X_first(cat_nodes), X_second(cat_nodes), X_third(cat_nodes)];
  labels = NodeLabels(cat_nodes);
  plot_grouped_bars_no_errors(...
      labels, Y, {'Hypo','Normal','Hyper'}, [COL_HYPO; COL_NORM; COL_HYPER], ...
      category_title, fullfile(fig_dir, sprintf('group_%s', sanitize_name(category_title))), ...
      YLIM_ACTIVATION, false, 50);
end

%% --------- PART F: All-nodes comparisons with grouped bars ------------
plot_grouped_bars_no_errors(NodeLabels, [SS_hypo_mu(:), SS_norm_mu(:)], {'Hypo','Normal'}, ...
    [COL_HYPO; COL_NORM], 'All nodes: Hypo vs Normal', ...
    fullfile(fig_dir,'allnodes_Hypo_vs_Normal'), YLIM_ACTIVATION, false, 150);

plot_grouped_bars_no_errors(NodeLabels, [SS_norm_mu(:), SS_hype_mu(:)], {'Normal','Hyper'}, ...
    [COL_NORM; COL_HYPER], 'All nodes: Normal vs Hyper', ...
    fullfile(fig_dir,'allnodes_Normal_vs_Hyper'), YLIM_ACTIVATION, false, 150);

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
            pos = [1 1 13.5 13.5];
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
            ax.Position = [0.28 0.10 0.67 0.84];
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

function plot_grouped_bars_no_errors(xlabels, Y, condOrder, condColors, titleStr, savepath, ylims, isObserved, maxNodesPerFig)
    if nargin < 7 || isempty(ylims), ylims = []; end
    if nargin < 8 || isempty(isObserved), isObserved = false; end
    if nargin < 9 || isempty(maxNodesPerFig)
        if isObserved
            maxNodesPerFig = 8;
        else
            maxNodesPerFig = 50;
        end
    end

    assert(size(Y,2) == numel(condOrder), 'Columns of Y must match condOrder');

    N = size(Y,1);
    nPages = ceil(N / maxNodesPerFig);

    for pg = 1:nPages
        idx1 = (pg-1)*maxNodesPerFig + 1;
        idx2 = min(pg*maxNodesPerFig, N);
        idx  = idx1:idx2;

        xlabels_pg = xlabels(idx);
        Y_pg       = Y(idx,:);
        Npg        = numel(idx);

        FS = get_pub_style();

        if Npg <= 10
            f = create_pub_figure('quarter');
        else
            f = create_pub_figure('quarterwide');
        end

        ax = axes('Parent',f);
        hold(ax,'on');

        x = 1:Npg;
        bh = bar(ax, x, Y_pg, 'grouped', 'BarWidth', 0.86);

        for k = 1:numel(bh)
            bh(k).FaceColor = condColors(k,:);
            bh(k).EdgeColor = 'none';
            bh(k).LineWidth = 0.5;
        end

        ylabel(ax,'Activation (a.u.)','FontSize',FS.label,'FontWeight','bold');

        if nPages > 1
            title_here = sprintf('%s (%d/%d)', titleStr, pg, nPages);
            savepath_here = sprintf('%s_part_%02d', savepath, pg);
        else
            title_here = titleStr;
            savepath_here = savepath;
        end

        title(ax, title_here, 'Interpreter','none', 'FontSize',FS.title, 'FontWeight','bold');

        if ~isempty(ylims)
            ylim(ax, ylims);
        end

        xlim(ax, [0.35, Npg + 0.65]);

        labels = make_readable_labels(xlabels_pg, 10, true);

        if Npg <= 8
            rot = 0;
            kind = 'bar';
        elseif Npg <= 14
            rot = 45;
            kind = 'bar_rot45';
        else
            rot = 60;
            kind = 'bar_rot60';
        end

        set(ax, ...
            'XTick', 1:Npg, ...
            'XTickLabel', labels, ...
            'XTickLabelRotation', rot);

        ax.TickLabelInterpreter = 'none';

        legend(ax, condOrder, ...
            'Location','northoutside', ...
            'Orientation','horizontal', ...
            'Box','off', ...
            'FontSize',FS.legend, ...
            'FontWeight','bold');

        apply_pub_axes(ax, FS, kind);

        savefig_safe(f, savepath_here);
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