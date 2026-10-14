"""Generated typed DPO projection for the Brusselator--ROP model.

The atom-resolved molecular rules remain the scientific input.  Their complete
component inventories generate small population-token rules; the polymer part
reuses the terminal-word projection.  No observable jump vector or published
ODE is supplied to this module.
"""
from __future__ import annotations
from mfobs.core.graph import node_labels
from mfobs.core.rules import make_rule, rule_graph
from mfobs.core.conditions import engine_conditions, engine_no_edges

from mfobs.core.composition_jump_model import GraphReactionClass

from mfobs.core.gml import add_rule_context
from mfobs.core.jump_closure_worklist import RuleExpression
from mfobs.core.polynomial_jump_closure import Polynomial, PolynomialRing
from mfobs.core.recursive_words import TerminalGrowthGrammar
from mfobs.models.aba.projection import TerminalWordGraphModel
from mfobs.models.population_projection import (
    PopulationChannel,
    TypedTokenGraphModel,
    molecular_component_inventory,
)
from mfobs.models import brusselator_growth_chemistry as chemistry


CORE = (
    chemistry.SOURCE,
    chemistry.AUTOCATALYSIS,
    chemistry.CONVERSION,
    chemistry.SINK,
)


class BrusselatorGrowthGraphModel:
    """Compose generated population and terminal-word structural rules."""

    def __init__(self, *, input_definition_sha256: str = "", state_space=None):
        from mfobs.models.brusselator_growth_chemistry import state_family as brusselator
        expected = brusselator()
        self.omega = expected if state_space is None else state_space.require(expected)
        domain = self.omega.definition
        tags = tuple(TypedTokenGraphModel.tag(name) for name in domain["molecules"])
        grammar = TerminalGrowthGrammar(
            "brusselator-growth",
            domain["polymer_alphabet"],
            (),
            False,
            input_definition_sha256,
        )
        self.word = TerminalWordGraphModel(grammar, extra_tags=tags)
        self.population = TypedTokenGraphModel(
            self.word,
            dynamic_coordinates=dict(domain["dynamic_coordinates"]),
            persistent_components=domain["persistent_components"],
            discarded_components=domain["discarded_components"],
        )
        self.input_definition_sha256 = input_definition_sha256

    def molecular_catalogue(self) -> dict[str, object]:
        return {name: chemistry.species(name) for name in self.omega.definition["molecules"]}

    def population_channels(self) -> tuple[PopulationChannel, ...]:
        catalogue = self.molecular_catalogue()
        rules = chemistry.core_rules()
        return tuple(
            PopulationChannel(
                identifier=name,
                channel=chemistry.RATE_CHANNELS[name],
                left=molecular_component_inventory(rule_graph(rules[name], "left"), catalogue),
                right=molecular_component_inventory(rule_graph(rules[name], "right"), catalogue),
            )
            for name in CORE
        )

    def reaction_classes(self) -> tuple[GraphReactionClass, ...]:
        result = [
            self.population.reaction_class(channel)
            for channel in self.population_channels()
        ]
        word_classes = {
            value.metadata["kind"]: value
            for value in self.word.reaction_classes()
        }
        feed = word_classes["feed"]
        result.append(GraphReactionClass(
            "monomer-inlet",
            chemistry.RATE_CHANNELS[chemistry.MONOMER_INLET],
            feed.rule,
            chemistry.MONOMER_INLET,
            "the unique persistent monomer-feed marker",
            {
                **dict(feed.metadata),
                "kind": "terminal-word-feed",
            },
        ))
        append = word_classes["append"]
        growth = add_rule_context(
            append.rule,
            self.population.token_graph("X"),
        )
        result.append(GraphReactionClass(
            chemistry.GROWTH,
            chemistry.RATE_CHANNELS[chemistry.GROWTH],
            growth,
            chemistry.GROWTH,
            "one X token, one free lactone token, and the unique chain end",
            {
                **dict(append.metadata),
                "kind": "terminal-word-append-with-population-catalyst",
                "catalyst_coordinate": "X",
            },
        ))
        return tuple(result)

    def observable_expressions(self) -> dict[str, RuleExpression]:
        word = self.word.observable_expressions()
        return {
            "X": self.population.observable_expression("X"),
            "Y": self.population.observable_expression("Y"),
            "M": word["M_B"],
            "N": word["N_B"],
        }

    def reduce(self, rule):
        """Exact count of a computed input, without constructing events."""
        import sympy as sp
        from mfobs.core.identity import certified_input_identity
        from mfobs.core.gml import component_graphs, union
        from mfobs.core.pattern_support import Image, UnsupportedReduction
        from mfobs.models.aba.counting import reduce_input

        identity = certified_input_identity(rule, self.dangling_certificate())
        graph = rule_graph(identity, "left")
        no_edges = engine_no_edges(identity)
        adjacency = engine_conditions(identity)
        populations, words = [], []
        for component in component_graphs(graph):
            target = populations if any(label.startswith("Population:") for label in node_labels(component).values()) else words
            # Population tokens are isolated Position--kind--tag components on Ω.
            if target is populations and component.numVertices > 2:
                return Image(sp.Integer(0), {})
            target.append(component)
        pop_graph = union(*populations)
        if any(set(pair) & set(node_labels(pop_graph)) for pair in no_edges) or any(
                c.center in node_labels(pop_graph) for c in adjacency):
            raise UnsupportedReduction("conditional population token count")
        ring = PolynomialRing(("X", "Y", "M", "N"))
        count = self.population.activity_polynomial(ring, None,
            make_rule("population-count", pop_graph, pop_graph))
        expression = sum(sp.Rational(t["coefficient"]) *
            sp.prod(sp.Symbol(n)**p for n, p in t["monomial"].items()) for t in ring.record(count)["terms"])
        if words:
            graph = union(*words)
            reduced = reduce_input(make_rule("word-count", graph, graph,
                no_edges, adjacency))
            expression *= reduced.expression.subs({sp.Symbol("M_B"):sp.Symbol("M"),
                                                   sp.Symbol("N_B"):sp.Symbol("N")})
        return Image(sp.expand(expression), {})

    def compose(self, observable, source):
        from mfobs.core.engine import commutator
        from mfobs.core.composition import rule_labels
        from mfobs.core.gml import rule_side
        from mfobs.core.typed_incidence import TypedIncidencePredicate, typed_incidence_constraints
        from mfobs.models.aba.counting import structural_violation
        spec = self.specification()
        predicate = TypedIncidencePredicate(spec)
        return commutator(observable, source.rule, include_empty=False,
            constraints=typed_incidence_constraints(spec, *rule_labels(observable, source.rule)),
            rule_rejection=lambda r: structural_violation(r) or structural_violation(r,"right") or
                                     predicate.violation(rule_side(r,"left")) or
                                     predicate.violation(rule_side(r,"right")))

    def activity_polynomial(
        self,
        ring: PolynomialRing,
        reaction: GraphReactionClass,
        activity_rule,
    ) -> Polynomial:
        kind = reaction.metadata["kind"]
        if kind == "finite-component-population":
            return self.population.activity_polynomial(
                ring, reaction, activity_rule
            )
        if kind == "terminal-word-feed":
            return ring.constant(1)
        if kind == "terminal-word-append-with-population-catalyst":
            return ring.multiply(
                ring.coordinate("M"),
                ring.coordinate(reaction.metadata["catalyst_coordinate"]),
            )
        raise ValueError(f"unsupported generated activity class {kind!r}")

    def specification(self) -> dict:
        return self.word.specification()

    def dangling_certificate(self):
        """Token positions by the token invariant; all else by the typed theory."""
        from mfobs.core.typed_incidence import dangling_certificate

        return self.population.dangling_certificate(
            dangling_certificate(self.specification())
        )

    def event_partition_certificate(self) -> dict:
        return {
            "status": "VERIFIED_ALL_STATES",
            "method": (
                "complete-component population projection plus the exact "
                "terminal-word congruence"
            ),
            "argument": (
                "Each complete X/Y molecule maps to one typed population token; "
                "persistent A/B components map to unit boundary tokens and D/E "
                "products are explicitly discarded by the open-reactor quotient. "
                "The caprolactone and unique chain end map bijectively to the "
                "terminal-word feed/append classes. Token injectivity gives exact "
                "falling factorials."
            ),
            "finite_enumeration_used": False,
        }

    def projection_record(self) -> dict:
        channels = self.population_channels()
        molecular_factors = {
            name: str(value)
            for name, value in chemistry.rate_factors().items()
        }
        return {
            "schema_version": 1,
            "kind": "brusselator-populations-plus-terminal-word-DPO-projection",
            "population": self.population.projection_record(channels),
            "word": self.word.projection_record(),
            "atom_level_per_match_normalisation": molecular_factors,
            "all_state_match_bijection": (
                "Automorphism-normalized complete molecular matches correspond "
                "one-to-one to injective token matches; the polymer insertion "
                "correspondence is the terminal-word match bijection."
            ),
            "supplied_delta_table": False,
        }
