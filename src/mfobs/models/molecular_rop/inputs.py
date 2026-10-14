"""Molecular ROP: constructor-defined Ω, lactones, native rules and selected motifs."""
from __future__ import annotations
from mfobs.core.graph import edge_rows, make_graph, node_labels, node_rows
from mfobs.core.rules import make_rule
import molfoundry as mf
from collections import Counter
from dataclasses import dataclass
from functools import lru_cache
from hashlib import sha256
import json
from typing import Mapping

from mfobs.core.gml import graph_gml, rule_gml, union, shift, disjoint_union

from mfobs.core.recursive_words import TerminalGrowthGrammar, WordState
from mfobs.core.molecules import complete_with_h
from mfobs.state_space import StateSpace, alphabet_symbols, covalent_theory
from mfobs.core.conditions import identity_rule
from mfobs.core.pattern_support import Pattern
from mfobs.models.aba.inputs import SourceRule


def state_family(alphabet=("A", "B"), gated=False):
    """Declare Ω for this model, independently of its selected queries."""
    alphabet = alphabet_symbols(alphabet)
    if set(alphabet) - {"A", "B"}:
        raise ValueError("molecular constructors support only lactones A and B")
    if type(gated) is not bool:
        raise ValueError("gated must be a boolean")
    return StateSpace("molecular-paths", {
        "alphabet": alphabet, "gated": gated,
        "path": "one finite initiator-rooted polyester chain with a terminal OH",
        "repeats": {a: {"spacer_carbons": {"A": 4, "B": 5}[a]} for a in alphabet},
        "pools": "arbitrary finite intact lactones and one isolated feed marker per symbol",
        "gate": "one palladium coordination complex, with or without donor coordination" if gated else None,
    }, covalent_theory(coordination=gated, feed_symbols=alphabet))


SPACER_CARBONS = {a: repeat["spacer_carbons"]
                  for a, repeat in state_family().definition["repeats"].items()}


MONOMER_NAMES = {
    "A": "delta-valerolactone",
    "B": "epsilon-caprolactone",
}


@dataclass(frozen=True)
class RopConfiguration:
    """One declarative choice of the existing A/B ROP analysis family.

    The configuration controls only the analysis query and the optional
    NoEdge gate.  It deliberately does not introduce a new reaction DSL: the
    atom-resolved A/B chemistry below remains the sole molecular model.
    """

    identifier: str
    title: str
    alphabet: tuple[str, ...]
    gated: bool
    motifs: tuple[tuple[str, ...], ...]
    purpose: str

    def __post_init__(self) -> None:
        if not self.identifier:
            raise ValueError("ROP configuration identifier must be nonempty")
        if not self.alphabet:
            raise ValueError("ROP configuration alphabet must be nonempty")
        if len(self.alphabet) != len(set(self.alphabet)):
            raise ValueError("ROP configuration alphabet symbols must be unique")
        unknown_symbols = set(self.alphabet) - set(SPACER_CARBONS)
        if unknown_symbols:
            raise ValueError(
                "ROP configuration alphabet may contain only chemical repeat "
                f"symbols A and B, not {sorted(unknown_symbols)!r}"
            )
        if not isinstance(self.gated, bool):
            raise ValueError("ROP configuration gate must be a Boolean")
        for motif in self.motifs:
            if len(motif) < 2:
                raise ValueError(
                    "ROP motif queries must have length at least two; "
                    f"got {''.join(motif)!r}"
                )
            outside = set(motif) - set(self.alphabet)
            if outside:
                raise ValueError(
                    f"ROP motif {''.join(motif)!r} uses symbols outside the "
                    f"selected alphabet {self.alphabet!r}: {sorted(outside)!r}"
                )
        if len(self.motifs) != len(set(self.motifs)):
            raise ValueError("ROP motif queries must be unique")

    def record(self) -> dict:
        """Canonical, non-derived part of the analysis configuration."""

        return {
            "identifier": self.identifier,
            "title": self.title,
            "alphabet": list(self.alphabet),
            "motifs": ["".join(motif) for motif in self.motifs],
            "gate": self.gated,
            "purpose": self.purpose,
        }


VARIANTS = {
    "baseline": RopConfiguration(
        "baseline", "ROP baseline", ("B",), False, (),
        "single-monomer first-moment reference",
    ),
    "noedge": RopConfiguration(
        "noedge", "ROP with active NoEdge gate", ("B",), True, (),
        "isolates conditional-DPO gate semantics",
    ),
    "copolymer": RopConfiguration(
        "copolymer", "Sequence-resolved copolymer ROP", ("A", "B"), False,
        tuple(tuple(word) for word in ("AA", "AB", "BA", "BB", "ABA")),
        "chemical recursive-combinatorial bridge to the Rémy analysis",
    ),
    "combined": RopConfiguration(
        "combined", "Gated sequence-resolved copolymer ROP", ("A", "B"), True,
        tuple(tuple(word) for word in ("AA", "AB", "BA", "BB", "ABA")),
        "full model combining sequence statistics and an active NoEdge gate",
    ),
}


def variant(name: str) -> RopConfiguration:
    try:
        return VARIANTS[name]
    except KeyError as exc:
        raise ValueError(f"unknown ROP-family variant {name!r}") from exc


def configured_variant(
    *,
    motifs: tuple[str | tuple[str, ...], ...] | list[str | tuple[str, ...]],
    gated: bool,
    alphabet: tuple[str, ...] = ("A", "B"),
    title: str | None = None,
) -> RopConfiguration:
    """Build a reproducible custom query over the fixed literal A/B chemistry.

    ``motifs`` is intentionally the only sequence-language input.  The default
    analysis discovers its support from unrefined DPO generator images, without
    supplied events, jumps or a target moment basis.
    """

    normalized_alphabet = tuple(alphabet)
    normalized_motifs = tuple(
        tuple(motif) if isinstance(motif, str) else tuple(motif)
        for motif in motifs
    )
    words = "-".join("".join(motif) for motif in normalized_motifs) or "units"
    identifier = f"custom-{''.join(normalized_alphabet)}-{words}-{'gated' if gated else 'ungated'}"
    return RopConfiguration(
        identifier,
        title or "Custom sequence-resolved ROP",
        normalized_alphabet,
        gated,
        normalized_motifs,
        "user-configured sequence queries over the fixed literal A/B chemistry",
    )


def configuration(value: str | RopConfiguration) -> RopConfiguration:
    """Resolve a named preset or validate an explicit configuration."""

    if isinstance(value, str):
        return variant(value)
    if isinstance(value, RopConfiguration):
        return value
    raise TypeError(
        "ROP analysis selection must be a named preset or RopConfiguration"
    )


def effective_configuration(
    value: str | RopConfiguration,
    *,
    linear_limit: int,
) -> dict:
    """Return the complete analysis input, including the worklist budget."""

    if not isinstance(linear_limit, int) or isinstance(linear_limit, bool) or linear_limit < 1:
        raise ValueError("linear dimension limit must be a positive integer")
    selected = configuration(value)
    record = {
        "schema_version": 1,
        "kind": "ROP-family-analysis-configuration",
        "preset": value if isinstance(value, str) else None,
        "model": selected.record(),
        "linear_dimension_limit": linear_limit,
    }
    record["sha256"] = sha256(
        json.dumps(record, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    return record


@lru_cache(None)
def lactone_core(symbol: str) -> mf.Graph:
    """Heavy-atom cyclic ester distinguished by its methylene spacer length."""

    spacer = SPACER_CARBONS[symbol]
    nodes = {0: "C", 1: "O", 2: "O"}
    methylenes = list(range(3, 3 + spacer))
    nodes.update({node: "C" for node in methylenes})
    edges = [(0, 1, "="), (0, 2, "-"), (2, methylenes[0], "-")]
    edges.extend((left, right, "-") for left, right in zip(methylenes, methylenes[1:]))
    edges.append((methylenes[-1], 0, "-"))
    return make_graph(nodes, edges)


@lru_cache(None)
def monomer(symbol: str) -> mf.Graph:
    return complete_with_h(lactone_core(symbol))


def chain_heavy(word: tuple[str, ...]) -> tuple[mf.Graph, int]:
    """BnO-[A/B]_w-H and the terminal oxygen identifier."""

    if set(word) - set(SPACER_CARBONS):
        raise ValueError("unknown repeat-unit symbol")
    nodes = {index: "C" for index in range(7)}
    nodes[7] = "O"
    edges = [
        (index, (index + 1) % 6, "=" if index % 2 == 0 else "-")
        for index in range(6)
    ]
    edges.extend([(0, 6, "-"), (6, 7, "-")])
    terminal, nxt = 7, 8
    for symbol in word:
        spacer = SPACER_CARBONS[symbol]
        carbonyl, oxo = nxt, nxt + 1
        methylenes = list(range(nxt + 2, nxt + 2 + spacer))
        outgoing = methylenes[-1] + 1
        nodes.update({carbonyl: "C", oxo: "O", outgoing: "O"})
        nodes.update({node: "C" for node in methylenes})
        edges.extend([
            (terminal, carbonyl, "-"),
            (carbonyl, oxo, "="),
            (carbonyl, methylenes[0], "-"),
        ])
        edges.extend((left, right, "-") for left, right in zip(methylenes, methylenes[1:]))
        edges.append((methylenes[-1], outgoing, "-"))
        terminal, nxt = outgoing, outgoing + 1
    return make_graph(nodes, edges), terminal


@lru_cache(None)
def chain(word: tuple[str, ...] = ()) -> mf.Graph:
    return complete_with_h(chain_heavy(tuple(word))[0])


def repeat_fragment(word: tuple[str, ...], *, terminal: bool = False) -> mf.Graph:
    """Connected heavy-atom sequence pattern; optionally end-anchored by OH."""

    if not word:
        raise ValueError("a repeat fragment needs at least one unit")
    nodes = {0: "O"}
    edges = []
    incoming, nxt = 0, 1
    for symbol in word:
        spacer = SPACER_CARBONS[symbol]
        carbonyl, oxo = nxt, nxt + 1
        methylenes = list(range(nxt + 2, nxt + 2 + spacer))
        outgoing = methylenes[-1] + 1
        nodes.update({carbonyl: "C", oxo: "O", outgoing: "O"})
        nodes.update({node: "C" for node in methylenes})
        edges.extend([
            (incoming, carbonyl, "-"),
            (carbonyl, oxo, "="),
            (carbonyl, methylenes[0], "-"),
        ])
        edges.extend((left, right, "-") for left, right in zip(methylenes, methylenes[1:]))
        edges.append((methylenes[-1], outgoing, "-"))
        incoming, nxt = outgoing, outgoing + 1
    if terminal:
        nodes[nxt] = "H"
        edges.append((incoming, nxt, "-"))
    return make_graph(nodes, edges)


@lru_cache(None)
def catalyst(closed: bool = False) -> mf.Graph:
    """Minimal asymmetric hemilabile Pd fragment used by the software model."""

    nodes = {0: "Pd", 1: "P", 2: "C", 3: "C", 4: "N", 5: "O"}
    edges = [(0, 1, "coord"), (1, 2, "-"), (2, 3, "-"), (3, 4, "-"), (0, 5, "coord")]
    if closed:
        edges.append((0, 4, "coord"))
    return make_graph(nodes, edges)


def reservoir(symbol: str) -> mf.Graph:
    return make_graph({0: f"Feed-{symbol}"})


def propagation_sides(symbol: str) -> tuple[mf.Graph, mf.Graph]:
    """Atom-preserving net alcoholysis of one declared lactone."""

    core = lactone_core(symbol)
    offset = max(node_labels(core)) + 1
    left = union(core, make_graph({offset: "O", offset + 1: "H"}, [(offset, offset + 1, "-")]))
    edges = set(edge_rows(left))
    edges.remove((0, 2, "-"))
    edges.remove((offset, offset + 1, "-"))
    edges.update({(0, offset, "-"), (2, offset + 1, "-")})
    return left, make_graph(node_labels(left), edges)


def feed_rule(symbol: str) -> mf.Rule:
    left = reservoir(symbol)
    right = union(left, shift(monomer(symbol), 1))
    return make_rule(f"feed_{symbol}", left, right)


def insertion_rule(symbol: str, *, gated: bool) -> mf.Rule:
    left, right = propagation_sides(symbol)
    no_edges = ()
    if gated:
        offset = max(node_labels(left)) + 1
        context = shift(catalyst(False), offset)
        left, right = union(left, context), union(right, context)
        no_edges = ((offset, offset + 4),)
    return make_rule(f"insert_{symbol}", left, right, no_edges)


def donor_rules() -> tuple[mf.Rule, mf.Rule]:
    open_catalyst, closed_catalyst = catalyst(False), catalyst(True)
    return (
        make_rule("donor_bind", open_catalyst, closed_catalyst, ((0, 4),)),
        make_rule("donor_unbind", closed_catalyst, open_catalyst),
    )


def rules(selection: str | RopConfiguration) -> dict[str, mf.Rule]:
    selected = configuration(selection)
    result = {f"feed_{symbol}": feed_rule(symbol) for symbol in selected.alphabet}
    result.update({
        f"insert_{symbol}": insertion_rule(symbol, gated=selected.gated)
        for symbol in selected.alphabet
    })
    if selected.gated:
        bind, unbind = donor_rules()
        result.update({bind.name: bind, unbind.name: unbind})
    return result


def state_graph(selection: str | RopConfiguration, state: WordState) -> mf.Graph:
    selected = configuration(selection)
    if len(state.pools) != len(selected.alphabet):
        raise ValueError("state pool dimension differs from the selected variant")
    components = [reservoir(symbol) for symbol in selected.alphabet]
    components.append(chain(state.word))
    for symbol, count in zip(selected.alphabet, state.pools):
        components.extend(monomer(symbol) for _ in range(count))
    if selected.gated:
        components.append(catalyst(not state.gate_open))
    return disjoint_union(components)


def _observable_rule(name: str, graph: mf.Graph, no_edges=()) -> mf.Rule:
    return make_rule(name, graph, graph, tuple(no_edges))


def observable_rules(grammar: TerminalGrowthGrammar) -> dict[str, mf.Rule]:
    result = {}
    for observable in grammar.observables():
        if observable.kind == "symbol-count":
            graph = repeat_fragment((observable.symbol,))
            result[observable.name] = _observable_rule(observable.name, graph)
        elif observable.kind == "factor-count":
            result[observable.name] = _observable_rule(
                observable.name, repeat_fragment(observable.word)
            )
        elif observable.kind == "pool-count":
            result[observable.name] = _observable_rule(
                observable.name, lactone_core(observable.symbol)
            )
        elif observable.kind == "gate-open":
            result[observable.name] = _observable_rule(
                observable.name, catalyst(False), ((0, 4),)
            )
        elif observable.kind == "suffix-class":
            result[observable.name] = _observable_rule(
                observable.name, repeat_fragment(observable.word, terminal=True)
            )
        else:
            raise AssertionError(f"unsupported observable kind {observable.kind!r}")
    return result


def exported_gml(selection: str | RopConfiguration) -> dict[str, Mapping[str, str]]:
    selected = configuration(selection)
    language = TerminalGrowthGrammar(
        selected.identifier, selected.alphabet, selected.motifs, selected.gated, "",
    )
    return {
        "rules": {identifier: rule_gml(rule) for identifier, rule in rules(selected).items()},
        "observables": {
            identifier: rule_gml(rule)
            for identifier, rule in observable_rules(language).items()
        },
        "example_states": {
            "empty": graph_gml(state_graph(selected, WordState(
                (), (1,) * len(language.alphabet), True
            ))),
            "AB": graph_gml(state_graph(selected, WordState(
                tuple(symbol for symbol in ("A", "B") if symbol in language.alphabet),
                (1,) * len(language.alphabet), True,
            ))),
        },
    }


def _identity(name, graph, no_edges=()):
    return identity_rule(graph, name, tuple(no_edges))


def source_rules(selection, *, state_space=None):
    selected = configuration(selection)
    expected = state_family(selected.alphabet, selected.gated)
    omega = expected if state_space is None else state_space.require(expected)
    sources = []
    for name, rule in rules(selected).items():
        kind = "feed" if name.startswith("feed_") else "append" if name.startswith("insert_") else "bind" if name == "donor_bind" else "unbind"
        symbol = name.rsplit("_", 1)[1] if kind in ("feed", "append") else None
        rate = ("beta_" if kind == "feed" else "k_") + symbol if symbol else "a" if kind == "bind" else "b"
        omega.validate_rule(rule)
        sources.append(SourceRule(rule=rule, rate=rate, identifier=name, kind=kind, symbol=symbol))
    return tuple(sources)


def query_pattern(word):
    name = ("N_" if len(word) == 1 else "C_") + word
    return Pattern(name, _identity(name, repeat_fragment(tuple(word))), "motif", {"word": word, "representation": "atom-resolved"})


def forbidden_specification(selection=None, *, maximum_candidates=2000000):
    selected = configuration(selection) if selection is not None else None
    omega = state_family(selected.alphabet, selected.gated) if selected else state_family(gated=True)
    specification = omega.chemical_constraints()
    if selected is not None:
        # Every constructor atom has one of these radius-one port profiles.
        # One-unit chains plus the empty initiator cover root/interior/terminal
        # atoms; the shared ester oxygen has the same C,C ports at all lengths.
        hosts = [chain(())] + [chain((a,)) for a in selected.alphabet]
        hosts += [monomer(a) for a in selected.alphabet]
        hosts += [reservoir(a) for a in selected.alphabet]
        if selected.gated: hosts += [catalyst(False),catalyst(True)]
        profiles = {}
        for host in hosts:
            omega.validate_graph(host)
            labels = node_labels(host)
            for atom,label in node_rows(host):
                profile = dict(Counter(bond+"|"+labels[v if u == atom else u] for u,v,bond in edge_rows(host) if atom in (u,v)))
                if profile not in profiles.setdefault(label,[]): profiles[label].append(profile)
        specification["allowedIncidentProfiles"] = profiles
    specification["maximumCandidateOverlaps"] = maximum_candidates
    return specification
