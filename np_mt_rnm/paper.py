"""Manuscript metadata, verbatim from the manuscript sources.

Source: NP_MT1.tex and NP_MT_SUP1.tex (manuscript version of 2026-11-26).
Text is copied word for word; only LaTeX markup is translated (``\\%`` -> ``%``,
``$\\Delta_i$`` -> ``Δ<sub>i</sub>``, ``\\kappa`` -> ``κ``). Used for the SBML model
notes, the webapp and the README so that all of them carry the manuscript's own
words. If the manuscript text changes, update this file.
"""
from __future__ import annotations

TITLE = (
    "A systems-level view of nucleus pulposus mechanotransduction and its "
    "regulatory control"
)

AUTHORS: tuple[str, ...] = (
    "Zerihun G. Workineh",
    "Francis Kiptengwer Chemorion",
    "Jérôme Noailly",
)
CORRESPONDING_AUTHOR = "Zerihun G. Workineh"
CORRESPONDING_EMAIL = "zerihungetahun.workineh@upf.edu"

AFFILIATION = (
    "Barcelona Centre for New Medical Technologies (BCN MedTech), Department of "
    "Engineering, Universitat Pompeu Fabra, Barcelona, 08018, Spain"
)

KEYWORDS: tuple[str, ...] = (
    "Nucleus pulposus",
    "Mechanotransduction",
    "Regulatory network modeling",
    "Mechanical loading",
    "Intervertebral disc degeneration",
    "Systems biology",
    "Network control",
)

ABSTRACT = (
    'Intervertebral disc degeneration (IDD), a major contributor to chronic low back '
    'pain, is strongly influenced by mechanical loading-dependent regulation of '
    'nucleus pulposus (NP) cell phenotype. Previous studies have characterized '
    'individual mechanosensors and specific signalling pathways involved in NP '
    'mechanotransduction. However, the regulatory interactions through which these '
    'components integrate mechanical cues into coordinated cellular responses remain '
    'incompletely understood. Here, we present an NP-specific, signed, directed '
    'regulatory network model (RNM) comprising 146 molecular nodes and 353 '
    'literature-supported interactions. The network connects mechanosensors with '
    'signalling cascades, metabolic and redox regulators, transcription factors, '
    'extracellular matrix (ECM) effectors, inflammatory mediators, and cell-fate '
    'pathways. Using the SQUADS semi-quantitative formalism, which accommodates '
    'limited kinetic parameter availability, we simulated steady-state responses to '
    'hypo-loading (Hypo), normal loading (NL), and hyper-loading (HL). We also '
    'analyzed network topology and screened 48 single- and dual-node perturbations '
    'under sustained HL. NL maintained high activity of the ECM-associated regulators '
    'ACAN, COL2A1, SOX9, and TIMP3, supported by AKT1, FAK, and antioxidant '
    'regulation. HL promoted inflammatory and catabolic regulation through mediators '
    'such as NF-κB, AP-1, MMP13, and ADAMTS4/5, accompanied by increased ROS, BAX, '
    'and CASP3 activity. Hypo reduced anabolic activity while retaining calcium, ROS, '
    'AMPK, and apoptotic signalling, suggesting a possible adhesion-dependent stress '
    'mechanism. Model falsification against 39 literature-derived molecular '
    'directionality benchmarks yielded 37 concordant predictions (94.9%), providing '
    'indirect, direction-level support for the RNM. HSF1 and CCN2 showed '
    'contradictory responses. Simulated regulatory recovery was module-specific. SOX9 '
    'activation achieved 82.8% ECM-profile recovery, increasing to 85.9% with ROS '
    'suppression. SOX9 activation combined with NF-κB suppression achieved 57.0% '
    'inflammatory/proteolytic-module recovery. FOXO activation with ROS suppression '
    'achieved 68.9% oxidative-stress/proteostasis-module recovery. NRF2 activation '
    'with PIEZO1 suppression achieved 37.5% cell-survival-module recovery. These '
    'predictions identify candidate molecular targets and intervention combinations '
    'for anabolic reinforcement, inflammatory suppression, redox protection, and '
    'mechanosensory modulation. The RNM provides a mechanistic framework for '
    'computational prioritization of candidate therapeutic targets and intervention '
    'combinations. Its predictions can guide subsequent preclinical validation in '
    'NP-cell and disc-explant models under controlled mechanical loading. '
    'Experimental validation, incorporation of time-dependent loading, and '
    'calibration of regulatory parameters are necessary to establish the therapeutic '
    "relevance of these predictions and improve the model's translational "
    'applicability.'
)

# Figure number and caption (HTML) for every computational figure, keyed by the
# manuscript's file name. Main-text figures 1-9 (Figures 1, 8 and 9 are
# schematics); supplementary figures S1-S3 in NP_MT_SUP1.tex order (S1 is the
# curation-workflow schematic).
_RESCUE_SET1 = (
    "Top 20 node-clamping strategies ranked by HL-to-NL profile recovery in four "
    "modules: (A) ECM anabolism and phenotype markers, (B) growth factors, (C) "
    "transcription factors, and (D) cytokines, chemokines, proteases, and related "
    "mediators. Each percentage is the reduction in the module-specific Euclidean "
    "distance to the NL reference relative to the unperturbed HL-to-NL distance."
)
_RESCUE_SET2 = (
    "Top 20 node-clamping strategies ranked by HL-to-NL profile recovery in four "
    "modules: (A) metabolic and related nodes, (B) ion channels and related nodes, "
    "(C) oxidative-stress defense and proteostasis, and (D) cell survival, "
    "apoptosis, mitophagy, and DNA-damage responses. Recovery percentages are "
    "defined relative to the original module-specific HL-to-NL Euclidean distance."
)
_FALS = (
    "RNM falsification results for 39 experimental benchmarks. (A) Mean "
    "steady-state node activation under NL and HL, with standard deviations (SD) "
    "across 100 replicate simulations per condition. (B) Directional comparison of "
    "RNM predictions and literature-derived responses. (C) Forest plot of "
    "RNM-calculated loading effects (Δ<sub>i</sub>=x̄<sub>i</sub><sup>NL</sup>−"
    "x̄<sub>i</sub><sup>HL</sup>) with 95% bootstrap confidence intervals (CIs)."
)
_NODE_RESOLVED_TAIL = (
    " Blue circles and orange squares identify the first- and second-ranked "
    "interventions within each module, respectively. Points indicate mean paired "
    "intervention-minus-HL activation changes; horizontal error bars show SD across "
    "100 matched simulations."
)
FIGURES: dict[str, tuple[str, str]] = {
    "TOPO_STATS": (
        "Figure 2",
        "Topological metrics of the NP mechanotransduction network. (A) nodes ranked "
        "by signed out-degree, including the number of inhibitory (red) and "
        "activating (blue) outgoing edges, (B) nodes ranked by betweenness "
        "centrality, and (C) nodes ranked by harmonic closeness centrality.",
    ),
    "BSL_ECM_GF": (
        "Figure 3",
        "Baseline mean activation under Hypo (blue), NL (green), and HL (red). Bars "
        "represent means across 100 replicate simulations per condition; error bars "
        "indicate ±1 SD. (A) Extracellular-matrix and phenotype markers; (B) growth "
        "factors; (C) transcription factors; and (D) cytokines, chemokines, "
        "proteases, and related mediators. Activations are bounded between 0 and 1; "
        "symmetric SD whiskers can extend beyond these bounds.",
    ),
    "BSL_MET_OXI": (
        "Figure 4",
        "Baseline mean activation under Hypo (blue), NL (green), and HL (red). Bars "
        "represent means across 100 replicate simulations per condition; error bars "
        "indicate ±1 SD. (A) Metabolic and related signaling nodes; (B) ion channels "
        "and calcium-associated nodes; (C) oxidative-stress defense and "
        "proteostasis; and (D) cell survival, apoptosis, and mitochondrial-stress "
        "markers. Activations are bounded between 0 and 1; symmetric SD whiskers can "
        "extend beyond these bounds.",
    ),
    "NP_MT_FALS_RECT": ("Figure 5", _FALS),
    "NP_MT_FALS": ("Figure 5", _FALS),
    "Top20_ECM_TrueRescue": ("Figure 6", _RESCUE_SET1),
    "Top20_GF_TrueRescue": ("Figure 6", _RESCUE_SET1),
    "Top20_TF_TrueRescue": ("Figure 6", _RESCUE_SET1),
    "Top20_CYT_TrueRescue": ("Figure 6", _RESCUE_SET1),
    "Top20_MET_TrueRescue": ("Figure 7", _RESCUE_SET2),
    "Top20_ION_TrueRescue": ("Figure 7", _RESCUE_SET2),
    "Top20_OXI_TrueRescue": ("Figure 7", _RESCUE_SET2),
    "Top20_CS_TrueRescue": ("Figure 7", _RESCUE_SET2),
    "NodeResolved_ECM_GF_TF_CYT": (
        "Figure S2",
        "Node-resolved responses for the two highest-ranked interventions in four "
        "functional modules: (A) ECM anabolism and phenotype markers, (B) growth "
        "factors, (C) transcription factors, and (D) cytokines, chemokines, "
        "proteases, and others." + _NODE_RESOLVED_TAIL,
    ),
    "NodeResolved_MET_ION_SURV_OXI": (
        "Figure S3",
        "Node-resolved responses for the two highest-ranked interventions in four "
        "functional modules: (A) metabolic and related regulators, (B) ion channels "
        "and related processes, (C) cell survival, apoptosis, mitophagy, and DNA "
        "damage, and (D) oxidative-stress defense and proteostasis."
        + _NODE_RESOLVED_TAIL,
    ),
}

# Supplementary Results, "Node-resolved intervention responses".
NODE_RESOLVED_INTRO = (
    "Figs. S2 and S3 show node-level responses for eight selected modules. Within "
    "each module, the two interventions with the highest module-level rescue "
    "scores were selected from the 48 screened strategies. Points show mean "
    "steady-state activation differences between each intervention and its paired "
    "HL replicate; error bars show the standard deviation of these differences "
    "across 100 simulations."
)

# Results, "Model falsification using independent NP mechanobiology experiments".
FALSIFICATION_SUMMARY = (
    "Of the 39 experimental benchmarks, 37 (94.9%) showed RNM-predicted effects "
    "concordant with their literature-derived directionality. Agreement was "
    "observed for all six benchmarks with positive literature-derived "
    "directionality (\"+\"; 6/6, 100%) and for 31 of the 33 benchmarks with "
    "negative directionality (\"−\"; 31/33, 93.9%). The remaining two benchmarks, "
    "HSF1 and CCN2, were contradictory (5.1%), and none was indeterminate."
)

# Results, "Hyper-to-Normal rescue screening via node clamping" (first paragraph).
RESCUE_SUMMARY = (
    "The 48 single- and dual-node clamping strategies produced different degrees "
    "of HL-to-NL profile recovery across the eight functional modules. The highest "
    "module-specific scores ranged from 10.4% for ion channels and related nodes "
    "to 85.9% for ECM anabolism and phenotype markers. Here, recovery denotes the "
    "percentage reduction in the Euclidean distance from the intervention profile "
    "to the NL reference, normalized to the original HL-to-NL distance within each "
    "module."
)
