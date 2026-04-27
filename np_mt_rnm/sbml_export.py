"""SBML Level 3 Version 2 export for the NP-MT regulatory network.

Adapted from RNMUpdate/rnm/sbml_export.py. The mathematical framework is
identical (Mendoza ODE = SQUADS with h=10, gamma=1, alpha=beta=1), only the
Network class and data file differ. Each node becomes one species, each
species gets a RateRule encoding the full ODE, and stimuli (Hypo, NL, HL)
are exposed as boundary species so a downstream loader can clamp them.

The generated SBML can be loaded by libroadrunner / tellurium / COPASI / any
SBML-Level-3 simulator.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

try:
    import libsbml
except ImportError:
    libsbml = None

from np_mt_rnm.network import Network


def _sanitize_id(name: str) -> str:
    """Convert protein names to valid SBML SId identifiers.

    Greek letters and special characters are spelled out; SIds must start
    with a letter or underscore and contain only [A-Za-z0-9_].
    """
    s = name
    replacements = {
        " ": "_", "-": "_", "/": "_", "(": "", ")": "", ",": "_", ".": "_",
        "α": "alpha", "β": "beta", "γ": "gamma",
        "κ": "kappa", "π": "pi",
        "+": "plus",
    }
    for old, new in replacements.items():
        s = s.replace(old, new)
    if s and s[0].isdigit():
        s = "_" + s
    return s


def _build_omega_formula(activator_sids: list[str], inhibitor_sids: list[str]) -> str:
    """Build the infix omega_n formula (Mendoza Eq. 4).

    Three regulator-presence cases. Uniform alpha=beta=1 weights, so the
    sum of weights equals the regulator count.
    """
    has_act = len(activator_sids) > 0
    has_inh = len(inhibitor_sids) > 0

    if has_act:
        n_act = len(activator_sids)
        sum_alpha_x = " + ".join(activator_sids)
        act_term = f"(({n_act} + 1) / {n_act}) * (({sum_alpha_x}) / (1 + ({sum_alpha_x})))"
    if has_inh:
        n_inh = len(inhibitor_sids)
        sum_beta_x = " + ".join(inhibitor_sids)
        inh_term = f"(({n_inh} + 1) / {n_inh}) * (({sum_beta_x}) / (1 + ({sum_beta_x})))"

    if has_act and not has_inh:
        return act_term
    if has_inh and not has_act:
        return f"(1 - {inh_term})"
    if has_act and has_inh:
        return f"({act_term}) * (1 - {inh_term})"
    return "0"


def _build_dxdt_formula(node_sid: str, omega_formula: str) -> str:
    """SQUADS / Mendoza state equation:

        dx/dt = sigmoid(omega; h) - gamma * x

    with h, gamma referenced as global parameters in the SBML model.
    """
    w = omega_formula
    return (
        f"((-exp(0.5 * h) + exp(-1 * h * (({w}) - 0.5))) / "
        f"((1 - exp(0.5 * h)) * (1 + exp(-1 * h * (({w}) - 0.5))))) - "
        f"gamma * {node_sid}"
    )


def export_sbml(
    network: Network,
    filepath: str | Path,
    initial_values: np.ndarray | None = None,
    h: float = 10.0,
    gamma: float = 1.0,
    model_id: str = "NP_MT_RNM",
    model_name: str = "NP Mechanotransduction Regulatory Network Model",
) -> bool:
    """Export the np-mt-rnm network as an SBML Level 3 Version 2 model.

    Parameters
    ----------
    network : Network
        Loaded np-mt-rnm Network (147 nodes).
    filepath : str | Path
        Output SBML file path.
    initial_values : ndarray | None
        Per-node initial concentrations. If None, all species default to 0.5
        (mid-range under the SQUADS sigmoid).
    h, gamma : float
        ODE parameters. Defaults match the paper (h=10, gamma=1).
    model_id, model_name : str
        SBML model identifier and human-readable name.

    Returns
    -------
    bool
        True on successful write.
    """
    if libsbml is None:
        raise ImportError(
            "python-libsbml is required for SBML export. "
            "Install with: pip install python-libsbml"
        )

    n = len(network.node_names)
    if initial_values is None:
        initial_values = np.full(n, 0.5)
    if initial_values.shape != (n,):
        raise ValueError(
            f"initial_values must have shape ({n},), got {initial_values.shape}"
        )

    document = libsbml.SBMLDocument(3, 2)
    model = document.createModel()
    model.setId(model_id)
    model.setName(model_name)

    # Time + dimensionless units.
    time_ud = model.createUnitDefinition()
    time_ud.setId("time_unit")
    u = time_ud.createUnit()
    u.setKind(libsbml.UNIT_KIND_SECOND)
    u.setExponent(1); u.setScale(0); u.setMultiplier(1)
    model.setTimeUnits("time_unit")

    dim_ud = model.createUnitDefinition()
    dim_ud.setId("arbitrary_units")
    u = dim_ud.createUnit()
    u.setKind(libsbml.UNIT_KIND_DIMENSIONLESS)
    u.setExponent(1); u.setScale(0); u.setMultiplier(1)

    compartment = model.createCompartment()
    compartment.setId("cell")
    compartment.setName("NP cell")
    compartment.setConstant(True)
    compartment.setSize(1.0)
    compartment.setSpatialDimensions(3)

    p_h = model.createParameter()
    p_h.setId("h"); p_h.setName("Sigmoid gain"); p_h.setValue(h); p_h.setConstant(True)
    p_gamma = model.createParameter()
    p_gamma.setId("gamma"); p_gamma.setName("Linear decay constant")
    p_gamma.setValue(gamma); p_gamma.setConstant(True)

    # SId sanitisation. Force uniqueness.
    sanitized = [_sanitize_id(name) for name in network.node_names]
    seen: dict[str, int] = {}
    for idx, sid in enumerate(sanitized):
        if sid in seen:
            sanitized[idx] = f"{sid}_{idx}"
        seen[sanitized[idx]] = idx

    # Boundary nodes = pure inputs with no activators AND no inhibitors.
    # These are the mechanical loading stimuli (Hypo, NL, HL in this network).
    is_boundary = np.array(
        [(network.mact[i, :].sum() == 0 and network.minh[i, :].sum() == 0) for i in range(n)]
    )

    # Species: all nodes. Boundary nodes get BoundaryCondition=true so a loader
    # can clamp them via setInitialConcentration / r["X"] = v.
    for i, (name, sid) in enumerate(zip(network.node_names, sanitized)):
        sp = model.createSpecies()
        sp.setId(sid)
        sp.setName(name)
        sp.setCompartment("cell")
        sp.setInitialConcentration(float(initial_values[i]))
        sp.setHasOnlySubstanceUnits(False)
        sp.setSubstanceUnits("arbitrary_units")
        sp.setBoundaryCondition(bool(is_boundary[i]))
        sp.setConstant(False)

    # Rate rules: every non-boundary species gets the Mendoza ODE.
    for i, (name, sid) in enumerate(zip(network.node_names, sanitized)):
        if is_boundary[i]:
            continue
        act_idx = np.where(network.mact[i, :] > 0)[0]
        inh_idx = np.where(network.minh[i, :] > 0)[0]
        act_sids = [sanitized[j] for j in act_idx]
        inh_sids = [sanitized[j] for j in inh_idx]

        omega = _build_omega_formula(act_sids, inh_sids)
        formula = _build_dxdt_formula(sid, omega)

        rule = model.createRateRule()
        rule.setVariable(sid)
        ast = libsbml.parseL3Formula(formula)
        if ast is None:
            print(f"WARN: failed to parse formula for {sid}: {formula}")
            continue
        rule.setMath(ast)

    notes = (
        '<body xmlns="http://www.w3.org/1999/xhtml">'
        f"<h1>{model_name}</h1>"
        f"<p>{n} nodes, {int(network.mact.sum())} activation edges, "
        f"{int(network.minh.sum())} inhibition edges.</p>"
        "<p>Mathematical framework: Mendoza/SQUADS ODE "
        "(Mendoza &amp; Xenarios 2006). Each non-stimulus species rate rule: "
        "dx/dt = sigmoid(omega; h) - gamma*x. "
        f"Parameters: h={h}, gamma={gamma}, alpha=beta=1.</p>"
        f"<p>Boundary (clampable input) species: "
        f"{', '.join(name for name, b in zip(network.node_names, is_boundary) if b) or '(none)'}.</p>"
        "</body>"
    )
    model.setNotes(notes)

    # Validate consistency.
    n_consistency = document.checkConsistency()
    if n_consistency > 0:
        n_real = 0
        for k in range(document.getNumErrors()):
            err = document.getError(k)
            if err.getSeverity() >= libsbml.LIBSBML_SEV_ERROR:
                print(f"SBML error: {err.getMessage()}")
                n_real += 1
        if n_real > 0:
            print(f"SBML validation: {n_real} fatal error(s)")

    success = libsbml.writeSBMLToFile(document, str(filepath))
    if success:
        print(f"SBML model written to {filepath}")
    else:
        print(f"Failed to write SBML to {filepath}")
    return bool(success)
