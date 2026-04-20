"""Topology metrics on the signed regulatory network.

Reproduces paper Figure 3: signed out-degree, betweenness centrality,
harmonic closeness centrality (all computed on the directed graph).
"""
from __future__ import annotations

from dataclasses import dataclass

import networkx as nx
import numpy as np

from np_mt_rnm.network import Network


@dataclass(frozen=True)
class TopologyMetrics:
    node_names: list[str]
    signed_out_degree: dict[str, int]
    betweenness: dict[str, float]
    harmonic_closeness: dict[str, float]
    out_activating: dict[str, int]
    out_inhibiting: dict[str, int]


def compute_topology(net: Network) -> TopologyMetrics:
    n = len(net.node_names)
    G = nx.DiGraph()
    G.add_nodes_from(net.node_names)
    # Edge i ← j means j regulates i. We add the edge j → i (regulator → target)
    # so out-degree counts outgoing regulatory influence.
    for i in range(n):
        for j in range(n):
            if net.mact[i, j]:
                G.add_edge(net.node_names[j], net.node_names[i], sign=+1)
            if net.minh[i, j]:
                G.add_edge(net.node_names[j], net.node_names[i], sign=-1)

    out_act = {node: 0 for node in G.nodes()}
    out_inh = {node: 0 for node in G.nodes()}
    for src, dst, data in G.edges(data=True):
        if data["sign"] > 0:
            out_act[src] += 1
        else:
            out_inh[src] += 1
    signed_out = {node: out_act[node] + out_inh[node] for node in G.nodes()}

    bet = nx.betweenness_centrality(G, normalized=True)
    # Harmonic closeness (out-closeness on the directed graph): how easily
    # a node reaches all others. Paper Fig 3C ranks broadcasters (ROS, NF-κB)
    # highest, so we compute on G.reverse() which makes harmonic_centrality
    # measure outgoing reach from the source node.
    hclose = nx.harmonic_centrality(G.reverse(copy=False))
    if n > 1:
        hclose = {k: v / (n - 1) for k, v in hclose.items()}

    return TopologyMetrics(
        node_names=list(net.node_names),
        signed_out_degree=signed_out,
        betweenness=bet,
        harmonic_closeness=hclose,
        out_activating=out_act,
        out_inhibiting=out_inh,
    )
