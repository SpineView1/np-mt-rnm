"""Manuscript metadata, verbatim from the manuscript sources.

Source: NP_MT1.tex and NP_MT_SUP1.tex (manuscript version of 2026-10-06).
Text is copied word for word; only LaTeX markup is translated (``\\%`` -> ``%``,
``$\\Delta$`` -> ``Δ``, ``\\textbf{A}`` -> ``A``). Used for the SBML model notes,
the webapp and the README so that all of them carry the manuscript's own words.
If the manuscript text changes, update this file.
"""
from __future__ import annotations

TITLE = (
    "A systems-level network model reveals mechanical regulation of nucleus "
    "pulposus cell states."
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

# The source writes "(~95\%)"; "~" is shown as written.
ABSTRACT = (
    "Intervertebral disc degeneration (IDD) is a leading cause of chronic low back "
    "pain and is strongly influenced by mechanical loading–dependent regulation of "
    "nucleus pulposus (NP) cell phenotype. Although individual mechanosensors and "
    "signaling pathways have been characterized, the systems-level principles "
    "governing how NP cells integrate mechanical cues into coordinated regulatory "
    "states remain unclear. Here, we present a systems-level regulatory network "
    "model of NP mechanotransduction comprising 147 molecular nodes and 356 "
    "experimentally supported interactions spanning mechanosensory inputs, "
    "signaling cascades, metabolic and redox regulators, transcription factors, "
    "extracellular matrix (ECM) effectors, inflammatory mediators, and cell-fate "
    "modules. Mechanical environments are represented as hypo-, normal-, and "
    "hyper-loading inputs that initiate distinct signaling programs. "
    "Semi-quantitative regulatory network simulations reveal three stable regimes: "
    "a hypo-loading state characterized by reduced ECM-anabolic activity, impaired "
    "adhesion-mediated survival signaling, and metabolic stress with features "
    "consistent with an anoikis-like phenotype; a normal-loading state associated "
    "with coordinated ECM maintenance, metabolic balance, and redox stability; and "
    "a hyper-loading state dominated by inflammatory amplification, oxidative "
    "stress, matrix degradation, and apoptosis. Falsification tests against "
    "independent experimental data demonstrate high directional concordance "
    "(~95%), supporting the biological plausibility of the network. Systematic "
    "perturbation analysis further identifies a distributed control structure in "
    "which mechanotransductive, redox, and transcriptional regulators jointly "
    "determine state stability, with combined attenuation of mechanically driven "
    "and stress-amplifying pathways together with activation of anabolic or "
    "cytoprotective programs most effectively restoring normal-like states. These "
    "findings provide a systems-level framework for understanding load-dependent "
    "NP cell regulation and for guiding multi-target therapeutic strategies in "
    "intervertebral disc mechanobiology."
)

# Figure number and caption (HTML) for every computational figure, keyed by the
# manuscript's file name. Main-text figures are numbered 1-9; supplementary
# figures S1-S17 in NP_MT_SUP1.tex order.
_NODE_RESOLVED = "Node-resolved perturbation responses within the {}."
FIGURES: dict[str, tuple[str, str]] = {
    "TOPO_STATS": (
        "Figure 3",
        "Topological metrics of the NP mechanotransduction network.(A) nodes ranked "
        "by signed out-degree,(B) nodes ranked by betweenness centrality, and (C) "
        "nodes ranked by harmonic closeness centrality.",
    ),
    "BSL_ECM_GF": (
        "Figure 4",
        "Baseline steady-state responses under Hypo (blue), Normal (green), and Hyper "
        "(red) loading. (<b>A</b>) ECM anabolism and phenotype markers; (<b>B</b>) "
        "Growth factors; (<b>C</b>) Transcription factors; (<b>D</b>) Cytokines, "
        "chemokines, proteases, and related mediators.",
    ),
    "BSL_MET_OXI": (
        "Figure 5",
        "Baseline steady-state responses under Hypo (blue), Normal (green), and Hyper "
        "(red) loading. (<b>A</b>) Metabolic and signaling intermediates; (<b>B</b>) "
        "Ion channels and calcium-dependent pathways; (<b>C</b>) Oxidative-stress "
        "defense and proteostatic responses; (<b>D</b>) Cell survival, apoptosis, and "
        "mitochondrial stress modules.",
    ),
    "ECM_transition": (
        "Figure 6",
        "Transition heatmaps of ECM anabolism and phenotype markers along constrained "
        "loading paths. (A) Hypo-to-Normal transition. (B) Normal-to-Hyper "
        "transition. Values represent mean steady-state activation levels at each "
        "path step.",
    ),
    "Representative_transition_paths": (
        "Figure 7",
        "Representative node trajectories along constrained loading paths. (A) "
        "Hypo-to-Normal transition. (B) Normal-to-Hyper transition. Selected nodes "
        "span ECM phenotype markers, transcriptional regulators, inflammatory and "
        "oxidative-stress mediators, growth factors, and cell survival/apoptosis "
        "regulators. Values represent mean steady-state activation levels at each "
        "path step.",
    ),
    "NP_MT_FALS": (
        "Figure 8",
        "Directional agreement and uncertainty in the Normal-to-Hyper falsification "
        "benchmark. (A) Directional effects for benchmarked nodes, expressed as "
        "Δ=x̄<sub>Normal</sub>−x̄<sub>Hyper</sub>. Dashed horizontal lines indicate "
        "the tolerance thresholds (±F<sub>TOL</sub>=0.02). (B) Bootstrap-derived 95% "
        "confidence intervals for Δ based on N<sub>boot</sub>=10,000 resamples.",
    ),
    "ECM_rescue": (
        "Figure 9",
        "Top 20 Hyper-to-Normal rescue interventions for ECM anabolism and phenotype "
        "markers. Interventions are ranked by the percentage reduction in Euclidean "
        "distance between the ensemble-mean perturbed module profile and the "
        "Normal-loading reference profile relative to the original Hyper-to-Normal "
        "distance. Larger positive values indicate greater restoration toward the "
        "Normal state.",
    ),
    "TF_transition": ("Figure S2", "Transition heatmaps of transcription-factor activity."),
    "GF_transition": ("Figure S3", "Transition heatmaps of growth-factor signaling."),
    "CYT_transition": (
        "Figure S4",
        "Transition heatmaps of cytokines, chemokines, proteases, and related "
        "inflammatory mediators.",
    ),
    "OXI_transition": (
        "Figure S5",
        "Transition heatmaps of oxidative-stress defense and proteostasis.",
    ),
    "CSF_transition": (
        "Figure S6",
        "Transition heatmaps of cell survival, apoptosis, mitophagy, and DNA-damage "
        "responses.",
    ),
    "GF_rescue": ("Figure S7", "Top 20 Hyper-to-Normal rescue interventions for the growth-factor module."),
    "TF_rescue": ("Figure S8", "Top 20 Hyper-to-Normal rescue interventions for the transcription-factor module."),
    "CYT_rescue": (
        "Figure S9",
        "Top 20 Hyper-to-Normal rescue interventions for the cytokine, chemokine, "
        "protease, and related inflammatory module.",
    ),
    "OXI_rescue": (
        "Figure S10",
        "Top 20 Hyper-to-Normal rescue interventions for the oxidative-stress defense "
        "and proteostasis module.",
    ),
    "CSF_rescue": (
        "Figure S11",
        "Top 20 Hyper-to-Normal rescue interventions for the cell-survival, "
        "apoptosis, mitophagy, and DNA-damage module.",
    ),
    "ECM_rescue1": ("Figure S12", _NODE_RESOLVED.format("ECM anabolism and phenotype-marker module")),
    "GF_rescue1": ("Figure S13", _NODE_RESOLVED.format("growth-factor module")),
    "TF_rescue1": ("Figure S14", _NODE_RESOLVED.format("transcription-factor module")),
    "CYT_rescue1": (
        "Figure S15",
        _NODE_RESOLVED.format("cytokine, chemokine, protease, and related inflammatory module"),
    ),
    "OXI_rescue1": ("Figure S16", _NODE_RESOLVED.format("oxidative-stress defense and proteostasis module")),
    "CSF_rescue1": (
        "Figure S17",
        _NODE_RESOLVED.format("cell-survival, apoptosis, mitophagy, and DNA-damage module"),
    ),
}

# Supplementary Section S6, introducing Figures S12-S17.
NODE_RESOLVED_INTRO = (
    "Node-resolved perturbation responses are provided for all six functional "
    "modules. For these figures, Δ denotes the paired steady-state change in each "
    "node relative to its corresponding Hyper-loading replicate. Panel A shows the "
    "mean node-level response across all screened perturbations. Panel B shows the "
    "10 highest-ranked perturbations according to the mean absolute node-level "
    "change within the indicated module. Panels C and D show the node-level "
    "response profiles for the highest-ranked dual- and single-node perturbations, "
    "respectively. Mean absolute change measures the magnitude of the perturbation "
    "response and is separate from the Hyper-to-Normal rescue metric used above."
)
