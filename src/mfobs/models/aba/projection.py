"""Directed terminal-word graph encoding for ROP projection comparisons."""
from __future__ import annotations
from mfobs.core.graph import make_graph, edge_rows, node_labels
from mfobs.core.rules import make_rule
from mfobs.core.engine import mf
from dataclasses import dataclass
from typing import Sequence
from mfobs.core.composition_jump_model import GraphReactionClass
from mfobs.core.conditions import Adjacency

from mfobs.core.jump_closure_worklist import RuleExpression
from mfobs.core.polynomial_jump_closure import Polynomial, PolynomialRing
from mfobs.core.recursive_words import TerminalGrowthGrammar, WordState


POSITION = "Position"


END = "End"


ROOT_TAG = "Root"


NEXT = "Next"


TERMINAL = "Terminal"


SOURCE_BOND = "src"


TARGET_BOND = "tgt"


KIND_BOND = "kind"


BOUND_BOND = "bound"


GATE_SOURCE = "GateSource"


GATE_TARGET = "GateTarget"


def symbol_tag(symbol: str) -> str:
    return f"Symbol:{symbol}"


def feed_tag(symbol: str) -> str:
    return f"Feed:{symbol}"


@dataclass
class TypedGraphBuilder:
    """Small public builder for the shared typed-incidence graph schema."""

    nodes: dict[int, str]
    edges: list[tuple[int, int, str]]
    next_id: int = 0

    @classmethod
    def empty(cls) -> "TypedGraphBuilder":
        return cls({}, [])

    def node(self, label: str) -> int:
        result = self.next_id
        self.next_id += 1
        self.nodes[result] = label
        return result

    def edge(self, source: int, target: int, label: str) -> None:
        self.edges.append((source, target, label))

    def relation(self, source: int, target: int, label: str) -> None:
        self.edge(source, target, label)

    def directed(self, source: int, target: int, edge_type: str) -> int:
        edge = self.node(edge_type)
        self.edge(source, edge, SOURCE_BOND)
        self.edge(edge, target, TARGET_BOND)
        return edge

    def graph(self) -> mf.Graph:
        return make_graph(self.nodes, self.edges)


@dataclass(frozen=True)
class _TerminalFragment:
    graph: mf.Graph
    last: int
    end: int
    terminal_edge: int


class TerminalWordGraphModel:
    """Generate structural rules and observables from one word grammar."""

    def __init__(
        self,
        grammar: TerminalGrowthGrammar,
        *,
        extra_tags: tuple[str, ...] = (),
    ):
        self.grammar = grammar
        self.extra_tags = tuple(extra_tags)
        if len(self.extra_tags) != len(set(self.extra_tags)):
            raise ValueError("extra tags must be unique")

    @property
    def tags(self) -> tuple[str, ...]:
        base = (ROOT_TAG,) + tuple(
            symbol_tag(value) for value in self.grammar.alphabet
        )
        if set(base) & set(self.extra_tags):
            raise ValueError("extra tags overlap terminal-word tags")
        return base + self.extra_tags

    def specification(self) -> dict:
        objects = [POSITION, END, *self.tags]
        objects.extend(feed_tag(value) for value in self.grammar.alphabet)
        relation_bonds = {
            KIND_BOND: [[POSITION, tag] for tag in self.tags],
        }
        relation_bounds = [
            {
                "center": POSITION,
                "bond": KIND_BOND,
                "neighbourLabels": list(self.tags),
                "maximum": 1,
                "name": "one-position-kind",
            },
            *(
                {
                    "center": tag,
                    "bond": KIND_BOND,
                    "neighbourLabels": [POSITION],
                    "maximum": 1,
                    "name": f"one-owner-{tag}",
                }
                for tag in self.tags
            ),
        ]
        if self.grammar.gated:
            objects.extend((GATE_SOURCE, GATE_TARGET))
            relation_bonds[BOUND_BOND] = [[GATE_SOURCE, GATE_TARGET]]
            relation_bounds.extend((
                {
                    "center": GATE_SOURCE,
                    "bond": BOUND_BOND,
                    "neighbourLabels": [GATE_TARGET],
                    "maximum": 1,
                    "name": "one-gate-bond-at-source",
                },
                {
                    "center": GATE_TARGET,
                    "bond": BOUND_BOND,
                    "neighbourLabels": [GATE_SOURCE],
                    "maximum": 1,
                    "name": "one-gate-bond-at-target",
                },
            ))
        return {
            "name": f"terminal-word-incidence-{self.grammar.name}",
            "objectLabels": objects,
            "edgeNodeLabels": [NEXT, TERMINAL],
            "sourceBond": SOURCE_BOND,
            "targetBond": TARGET_BOND,
            "atMostOneOutgoing": True,
            "atMostOneIncomingTypes": [NEXT, TERMINAL],
            "forbidSuccessorAfterTypes": [TERMINAL],
            "simpleDirected": True,
            "simpleUnderlying": True,
            "acyclicDirected": True,
            "edgeEndpointTypes": {
                NEXT: {"source": [POSITION], "target": [POSITION]},
                TERMINAL: {"source": [POSITION], "target": [END]},
            },
            "relationBonds": relation_bonds,
            "relationIncidentBounds": relation_bounds,
        }

    @staticmethod
    def tagged_position(builder: TypedGraphBuilder, tag: str) -> int:
        """Add one position with exactly one schema-level kind tag."""

        position = builder.node(POSITION)
        marker = builder.node(tag)
        builder.relation(position, marker, KIND_BOND)
        return position

    _tagged_position = tagged_position

    def _terminal_fragment(self, suffix: tuple[str, ...]) -> _TerminalFragment:
        depth = self.grammar.suffix_depth
        if suffix not in self.grammar.suffix_classes:
            raise ValueError("suffix is not a class of the grammar congruence")
        builder = TypedGraphBuilder.empty()
        if len(suffix) < depth:
            last = self._tagged_position(builder, ROOT_TAG)
        elif not suffix:
            # The unique depth-zero class covers every terminal position.
            last = builder.node(POSITION)
        else:
            last = None
        for symbol in suffix:
            current = self._tagged_position(builder, symbol_tag(symbol))
            if last is not None:
                builder.directed(last, current, NEXT)
            last = current
        if last is None:
            raise AssertionError("a terminal fragment needs a final position")
        end = builder.node(END)
        terminal = builder.directed(last, end, TERMINAL)
        return _TerminalFragment(builder.graph(), last, end, terminal)

    @staticmethod
    def copy_graph(graph: mf.Graph) -> TypedGraphBuilder:
        return TypedGraphBuilder(
            dict(node_labels(graph)),
            list(edge_rows(graph)),
            max(node_labels(graph), default=-1) + 1,
        )

    _copy = copy_graph

    def _gate_context(self, builder: TypedGraphBuilder) -> tuple[int, int]:
        source = builder.node(GATE_SOURCE)
        target = builder.node(GATE_TARGET)
        return source, target

    def feed_rule(self, symbol: str) -> mf.Rule:
        left_builder = TypedGraphBuilder.empty()
        left_builder.node(feed_tag(symbol))
        left = left_builder.graph()
        right_builder = self._copy(left)
        self._tagged_position(right_builder, symbol_tag(symbol))
        return make_rule(f"feed_{symbol}", left, right_builder.graph())

    def insert_rule(self, symbol: str, suffix: tuple[str, ...]) -> mf.Rule:
        fragment = self._terminal_fragment(suffix)
        left_builder = self._copy(fragment.graph)
        free = self._tagged_position(left_builder, symbol_tag(symbol))
        no_edges = ()
        if self.grammar.gated:
            gate_source, gate_target = self._gate_context(left_builder)
            no_edges = ((gate_source, gate_target),)
        left = left_builder.graph()

        right_nodes = dict(node_labels(left))
        right_nodes.pop(fragment.terminal_edge)
        right_edges = [
            edge for edge in edge_rows(left) if fragment.terminal_edge not in edge[:2]
        ]
        right_builder = TypedGraphBuilder(
            right_nodes,
            right_edges,
            max(node_labels(left), default=-1) + 1,
        )
        right_builder.directed(fragment.last, free, NEXT)
        right_builder.directed(free, fragment.end, TERMINAL)
        free_condition = Adjacency(
            free,
            (NEXT, TERMINAL),
            (SOURCE_BOND, TARGET_BOND),
            "=",
            0,
        )
        return make_rule(
            f"insert_{symbol}_{''.join(suffix) or 'empty'}",
            left,
            right_builder.graph(),
            no_edges,
            (free_condition,),
        )

    def donor_rules(self) -> tuple[mf.Rule, mf.Rule]:
        builder = TypedGraphBuilder.empty()
        source, target = self._gate_context(builder)
        opened = builder.graph()
        closed_builder = self._copy(opened)
        closed_builder.relation(source, target, BOUND_BOND)
        closed = closed_builder.graph()
        return (
            make_rule("donor_bind", opened, closed, ((source, target),)),
            make_rule("donor_unbind", closed, opened),
        )

    def reaction_classes(self) -> tuple[GraphReactionClass, ...]:
        """Generate the finite rule refinement as native MolFoundry rules."""

        result = []
        for symbol in self.grammar.alphabet:
            result.append(GraphReactionClass(
                f"feed-{symbol}",
                f"beta_{symbol}",
                self.feed_rule(symbol),
                f"feed_{symbol}",
                "the unique persistent feed marker",
                {"kind": "feed", "symbol": symbol},
            ))
        for suffix in self.grammar.suffix_classes:
            for symbol in self.grammar.alphabet:
                result.append(GraphReactionClass(
                    f"append-{symbol}-from-{''.join(suffix) or 'empty'}",
                    f"k_{symbol}",
                    self.insert_rule(symbol, suffix),
                    f"insert_{symbol}",
                    (
                        f"suffix={''.join(suffix) or 'empty'}"
                        + (" and NoEdge(GateSource,GateTarget)" if self.grammar.gated else "")
                    ),
                    {"kind": "append", "symbol": symbol, "suffix": list(suffix)},
                ))
        if self.grammar.gated:
            bind, unbind = self.donor_rules()
            result.extend((
                GraphReactionClass(
                    "donor-bind", "a", bind, "donor_bind",
                    "NoEdge(GateSource,GateTarget)", {"kind": "donor-bind"},
                ),
                GraphReactionClass(
                    "donor-unbind", "b", unbind, "donor_unbind",
                    "the gate relation edge is present", {"kind": "donor-unbind"},
                ),
            ))
        return tuple(result)

    def projection_record(self) -> dict:
        """Machine-readable record of the generated structural projection."""

        classes = self.reaction_classes()
        source_to_classes = {}
        for reaction in classes:
            source_to_classes.setdefault(reaction.source_rule, []).append(
                reaction.identifier
            )
        return {
            "schema_version": 1,
            "kind": "typed-terminal-word-DPO-projection",
            "source_to_structural_classes": source_to_classes,
            "structural_class_count": len(classes),
            "typed_incidence_specification": self.specification(),
            "construction": (
                "Source insertions are refined by the finite suffix congruence "
                "induced by the selected motif length. No jump vector or "
                "published equation is supplied to the projection."
            ),
        }

    @staticmethod
    def _identity(name: str, graph: mf.Graph, no_edges=(), adjacency=()):
        return make_rule(
            name, graph, graph, tuple(no_edges), tuple(adjacency)
        )

    def _symbol_observable(self, symbol: str, *, free: bool):
        builder = TypedGraphBuilder.empty()
        position = self._tagged_position(builder, symbol_tag(symbol))
        condition = Adjacency(
            position,
            (NEXT, TERMINAL),
            (SOURCE_BOND, TARGET_BOND),
            "=" if free else ">",
            0,
        )
        return self._identity(
            f"M_{symbol}" if free else f"N_{symbol}",
            builder.graph(),
            adjacency=(condition,),
        )

    def _word_fragment(self, word: tuple[str, ...], *, terminal: bool) -> mf.Graph:
        if not word:
            raise ValueError("word observables must be nonempty")
        builder = TypedGraphBuilder.empty()
        last = None
        for symbol in word:
            current = self._tagged_position(builder, symbol_tag(symbol))
            if last is not None:
                builder.directed(last, current, NEXT)
            last = current
        if terminal:
            end = builder.node(END)
            builder.directed(last, end, TERMINAL)
        return builder.graph()

    def observable_rules(self) -> dict[str, mf.Rule]:
        """Generate the selected counting observables as native MolFoundry rules."""

        result = {}
        for observable in self.grammar.observables():
            if observable.kind == "symbol-count":
                rule = self._symbol_observable(observable.symbol, free=False)
            elif observable.kind == "pool-count":
                rule = self._symbol_observable(observable.symbol, free=True)
            elif observable.kind == "factor-count":
                graph = self._word_fragment(observable.word, terminal=False)
                rule = self._identity(observable.name, graph)
            elif observable.kind == "suffix-class":
                graph = self._word_fragment(observable.word, terminal=True)
                rule = self._identity(observable.name, graph)
            elif observable.kind == "gate-open":
                builder = TypedGraphBuilder.empty()
                source, target = self._gate_context(builder)
                rule = self._identity(
                    observable.name,
                    builder.graph(),
                    no_edges=((source, target),),
                )
            else:
                raise AssertionError(f"unsupported word observable {observable.kind!r}")
            result[observable.name] = rule
        return result

    def observable_expressions(self) -> dict[str, RuleExpression]:
        return {
            name: RuleExpression.from_rules((rule,))
            for name, rule in self.observable_rules().items()
        }

    def activity_polynomial(
        self,
        ring: PolynomialRing,
        reaction: GraphReactionClass,
        _activity_rule,
    ) -> Polynomial:
        """Exact all-word count of a generated structural activity identity."""

        metadata = reaction.metadata
        kind = metadata["kind"]
        if kind == "feed":
            return ring.constant(1)
        if kind == "append":
            suffix = tuple(metadata["suffix"])
            propensity = ring.multiply(
                ring.coordinate(f"M_{metadata['symbol']}"),
                self.grammar.suffix_indicator(ring, suffix),
            )
            if self.grammar.gated:
                propensity = ring.multiply(propensity, ring.coordinate("X"))
            return propensity
        if kind == "donor-bind":
            return ring.coordinate("X")
        if kind == "donor-unbind":
            return ring.subtract(ring.constant(1), ring.coordinate("X"))
        raise AssertionError(f"unsupported reaction class {kind!r}")

    def event_partition_certificate(self) -> dict:
        return {
            "status": "VERIFIED_ALL_STATES",
            "status_scope": (
                "model-level certificate required by the pJC constructor; its "
                "all-state coverage rests on the stated terminal-word argument, "
                "not on a generic machine-checked proof of arbitrary graph languages"
            ),
            "method": (
                "typed terminal-path rule refinement by the finite suffix "
                "congruence induced mechanically by selected motif length"
            ),
            "alphabet": list(self.grammar.alphabet),
            "motifs": ["".join(value) for value in self.grammar.motifs],
            "suffix_depth": self.grammar.suffix_depth,
            "suffix_classes": [
                "".join(value) or "empty" for value in self.grammar.suffix_classes
            ],
            "argument": (
                "Every valid state has one terminal incidence edge. For depth k, "
                "its terminal path has exactly one class: the whole word when "
                "shorter than k, otherwise its length-k suffix. The generated "
                "insert-rule refinements are therefore disjoint and exhaustive; "
                "the free-position adjacency condition and optional NoEdge guard "
                "preserve source-rule match multiplicity."
            ),
            "enumeration_used": False,
            "finite_enumeration_used": False,
            "source_rules_refined_not_replaced": True,
        }

    def refinement_consistency_check(
        self,
        reactions: Sequence[GraphReactionClass] | None = None,
    ) -> dict:
        """Check the finite generated class inventory, not host semantics.

        This mechanically checks that the construction contains exactly the
        generated source-rule classes.  It deliberately does not establish the
        terminal-word induction argument, molecular/structural correspondence,
        or a property of an arbitrary graph language.
        """

        reactions = tuple(
            self.reaction_classes() if reactions is None else reactions
        )
        expected = {
            **{
                f"feed_{symbol}": (f"feed-{symbol}",)
                for symbol in self.grammar.alphabet
            },
            **{
                f"insert_{symbol}": tuple(
                    f"append-{symbol}-from-{''.join(suffix) or 'empty'}"
                    for suffix in self.grammar.suffix_classes
                )
                for symbol in self.grammar.alphabet
            },
        }
        if self.grammar.gated:
            expected.update({
                "donor_bind": ("donor-bind",),
                "donor_unbind": ("donor-unbind",),
            })
        actual: dict[str, tuple[str, ...]] = {}
        for reaction in reactions:
            actual.setdefault(reaction.source_rule, ())
            actual[reaction.source_rule] += (reaction.identifier,)
        failures = []
        for source_rule in sorted(set(expected) | set(actual)):
            wanted = tuple(sorted(expected.get(source_rule, ())))
            found = tuple(sorted(actual.get(source_rule, ())))
            if found != wanted:
                failures.append({
                    "source_rule": source_rule,
                    "expected_classes": list(wanted),
                    "actual_classes": list(found),
                })
        return {
            "status": "PASS" if not failures else "FAIL",
            "scope": (
                "finite class-identifier and source-rule inventory check of the "
                "generated terminal-word refinement only; it does not inspect "
                "rule guards, establish semantic disjointness, or prove an "
                "all-state graph-language property"
            ),
            "alphabet": list(self.grammar.alphabet),
            "motifs": ["".join(word) for word in self.grammar.motifs],
            "maximum_motif_length": max(
                (len(word) for word in self.grammar.motifs), default=1
            ),
            "suffix_depth": self.grammar.suffix_depth,
            "expected_source_rules": {
                source: list(classes) for source, classes in expected.items()
            },
            "actual_source_rules": {
                source: list(classes) for source, classes in actual.items()
            },
            "rule_predicates_checked": False,
            "semantic_disjointness_checked": False,
            "failures": failures,
        }

    def state_graph(self, state: WordState) -> mf.Graph:
        if len(state.pools) != len(self.grammar.alphabet):
            raise ValueError("pool dimension differs from the grammar alphabet")
        builder = TypedGraphBuilder.empty()
        last = self._tagged_position(builder, ROOT_TAG)
        for symbol in state.word:
            current = self._tagged_position(builder, symbol_tag(symbol))
            builder.directed(last, current, NEXT)
            last = current
        end = builder.node(END)
        builder.directed(last, end, TERMINAL)
        for symbol, count in zip(self.grammar.alphabet, state.pools):
            for _ in range(count):
                self._tagged_position(builder, symbol_tag(symbol))
            builder.node(feed_tag(symbol))
        if self.grammar.gated:
            source, target = self._gate_context(builder)
            if not state.gate_open:
                builder.relation(source, target, BOUND_BOND)
        return builder.graph()
