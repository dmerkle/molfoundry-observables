"""Selected intramolecular [2+2] channel with an active context NoEdge.

The ordinary-covalent A -> P transformation is 1,5-hexadiene to
bicyclo[2.1.1]hexane. This crossed photocycloaddition has experimental precedent:
Srinivasan and Sonntag, JACS 1967, 89, 407-410, DOI 10.1021/ja00978a042;
https://research.ibm.com/publications/bicyclo-211-hexane-preparation-and-photochlorination

This is ONE selected net channel, not a complete or quantitatively calibrated
photochemical network. Published chemistry also gives other isomers and side
products. No elementary-mechanism or exact intermediate-elimination claim is
made. Illumination/sensitization is absorbed into the illustrative rate.

B is a valence-valid context challenge, 1,2-dimethylenecyclobutane. The rule's
declared context guard excludes it; removing that guard produces D, a
valence-valid C6H8 tricyclic graph. Neither B's inertness nor its ablated D
channel is claimed to be experimentally established. RDKit acceptance proves
graph valence consistency, not thermodynamic stability or synthetic access.

Arbitrary mixtures of A, B, P and D define a finite-species, unbounded-population
language. Its NoEdge guard is active, but can be compiled positively on this
finite language; no separation from FRAGGER is implied.
"""
from __future__ import annotations
from mfobs.core.graph import edge_rows, make_graph, node_labels
from mfobs.core.rules import make_rule
import molfoundry as mf

from fractions import Fraction
from functools import lru_cache


from mfobs.models.nonadjacent_chemistry import PRODUCT_NO_EDGES, CONTEXT_NO_EDGES
from mfobs.core.gml import disjoint_union
from mfobs.core.molecules import complete_with_h


PER_MATCH = Fraction(1, 2)
RATE_FACTORS = {'cyclize': PER_MATCH}
RATE_NAMES = {'cyclize': 'k'}


def reaction_sides() -> tuple[mf.Graph, mf.Graph]:
    """Six-carbon positive pattern includes the entire two-carbon tether."""
    nodes = {u: 'C' for u in range(6)}
    tether = [(0, 4, '-'), (4, 5, '-'), (5, 2, '-')]
    left = make_graph(nodes, [(0, 1, '='), (2, 3, '='), *tether])
    right = make_graph(nodes, [(0, 1, '-'), (2, 3, '-'),
                              (0, 3, '-'), (1, 2, '-'), *tether])
    return left, right


@lru_cache(None)
def guarded_rule() -> mf.Rule:
    return make_rule('tethered-isolated-alkene-crossed-photocycloaddition',
                           *reaction_sides(), PRODUCT_NO_EDGES + CONTEXT_NO_EDGES)


@lru_cache(None)
def ablated_rule() -> mf.Rule:
    return make_rule('tethered-context-guard-ablated-photocycloaddition',
                           *reaction_sides(), PRODUCT_NO_EDGES)


def rules() -> dict[str, mf.Rule]:
    return {'cyclize': guarded_rule()}


@lru_cache(None)
def precursor() -> mf.Graph:
    """A: C6H10, 1,5-hexadiene, two injective positive rule maps."""
    return complete_with_h(reaction_sides()[0])


@lru_cache(None)
def blocked() -> mf.Graph:
    """B: C6H8, with the existing unchanged context edge 0-2."""
    left, _ = reaction_sides()
    return complete_with_h(make_graph(node_labels(left), [*edge_rows(left), (0, 2, '-')]))


@lru_cache(None)
def product() -> mf.Graph:
    """P: C6H10, crossed bicyclo[2.1.1]hexane (not bicyclo[2.2.0])."""
    return complete_with_h(reaction_sides()[1])


@lru_cache(None)
def blocked_product() -> mf.Graph:
    """D: the C6H8 ablation successor, included so ablated tests stay closed."""
    _, right = reaction_sides()
    return complete_with_h(make_graph(node_labels(right), [*edge_rows(right), (0, 2, '-')]))


def language_hosts() -> dict[str, mf.Graph]:
    return {'A': precursor(), 'B': blocked(), 'P': product(), 'D': blocked_product()}


def host(populations: dict[str, int]) -> mf.Graph:
    allowed = language_hosts()
    if set(populations) - set(allowed):
        raise ValueError('Unknown molecular species in finite-component language')
    if any(not isinstance(n, int) or isinstance(n, bool) or n < 0 for n in populations.values()):
        raise ValueError('Populations must be nonnegative integers')
    return disjoint_union([graph for name, graph in allowed.items()
                           for _ in range(populations.get(name, 0))])


def seed_graph() -> mf.Graph:
    """One ordinary C=C graph observable, with four embeddings in A and B."""
    return make_graph({0: 'C', 1: 'C'}, [(0, 1, '=')])
