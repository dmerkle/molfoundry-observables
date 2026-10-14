"""Direct molecular applications and comparison with the typed-path result."""
from __future__ import annotations
import molfoundry as mf

import mfobs.models.aba.events as terminal_events
from collections import Counter
from itertools import product

import sympy as sp

from mfobs.core.gml import rule_side
from mfobs.core.mixtures import fast_count
from mfobs.models.molecular_rop import inputs as rop
from mfobs.models.molecular_rop.events import polynomial_jump_model

def validate_molecular_hosts(selection, support, sources, *, jump_model, maximum_length=3, extra_words=()):
    """Atom-graph matching/application and full joint-jump checks.

    No word successor formula is used to construct observed transitions.
    MolFoundry enumerates the rule maps, applies each one, and counts the
    selected AND newly discovered patterns.
    """
    from mfobs.core.conditions import matches
    from mfobs.core.graph import apply_match
    from mfobs.core.recursive_words import WordState
    from mfobs.checks.terminal_paths import coordinate_values
    selected = rop.configuration(selection)
    if support.status != "POLYNOMIAL_SUPPORT_CLOSED":
        return {"status": "NOT_EXECUTED", "reason": support.status}
    atom_rules = rop.rules(selected)
    patterns = {}
    for name, pattern in support.patterns.items():
        patterns[name] = (rule_side(pattern.rule,"left"), tuple((c.source, c.target) for c in pattern.rule.matchConstraintsNoEdge))
    names = tuple(patterns)
    counts_cache = {}
    def counts(host):
        if host not in counts_cache:
            counts_cache[host] = {n: fast_count(g,host,nac) for n,(g,nac) in patterns.items()}
        return counts_cache[host]
    states = rule_checks = image_checks = distributions = target_checks = occurrence_checks = 0
    failures = []
    words = {tuple(w) for length in range(maximum_length+1) for w in product(selected.alphabet,repeat=length)}
    words.update(tuple(w) for w in extra_words)
    for word in sorted(words,key=lambda w:(len(w),w)):
        if set(word)-set(selected.alphabet):
            raise ValueError("validation word outside the selected alphabet")
        for pools in product((0,1,2),repeat=len(selected.alphabet)):
            for opened in ((False,True) if selected.gated else (True,)):
                host = rop.state_graph(selected,WordState(word,pools,opened))
                before = counts(host)
                word_counts = coordinate_values(support.patterns,"".join(word),dict(zip(selected.alphabet,pools)),opened)
                if before != word_counts:
                    failures.append({"kind":"coordinate-correspondence","word":"".join(word),"pools":pools,"gate_open":opened,"molecular":before,"word_counts":word_counts})
                occurrence_checks += len(names)
                substitutions = {sp.Symbol(n):v for n,v in before.items()}
                for source in sources:
                    rule = atom_rules[source.identifier]
                    observed = Counter(); number = 0; changes = Counter()
                    for mapping in matches(rule, host):
                        target = apply_match(rule,host,mapping)
                        if target is None:
                            continue
                        number += 1
                        after = counts(target)
                        delta = tuple(after[n]-before[n] for n in names)
                        if any(delta): observed[delta] += 1
                        for n,change in zip(names,delta): changes[n] += change
                        target_checks += 1
                    expected_activity = support.activities[source.identifier].subs(substitutions)
                    if number != expected_activity:
                        failures.append({"kind":"source-activity","source":source.identifier,"actual":number,"expected":str(expected_activity),"word":"".join(word),"pools":pools,"gate_open":opened})
                    for n in names:
                        expected_image = support.images[n][source.identifier].subs(substitutions)
                        if changes[n] != expected_image:
                            failures.append({"kind":"generator-image","source":source.identifier,"coordinate":n,"actual":changes[n],"expected":str(expected_image),"word":"".join(word),"pools":pools,"gate_open":opened})
                        image_checks += 1
                    predicted = Counter()
                    for event in jump_model.events:
                        if event.rule == source.identifier:
                            rate = jump_model.ring.evaluate(event.propensity,before)
                            if rate < 0:
                                raise AssertionError("negative event activity on a molecular host")
                            if rate: predicted[event.delta] += rate
                    if observed != predicted:
                        failures.append({"kind":"joint-jump-distribution","source":source.identifier,"actual":repr(observed),"expected":repr(predicted),"word":"".join(word),"pools":pools,"gate_open":opened})
                    rule_checks += 1; distributions += 1
                states += 1
    return {"status":"PASS" if not failures else "FAIL", "states":states,
        "maximum_exhaustive_word_length":maximum_length,"additional_words":["".join(w) for w in sorted(words) if len(w)>maximum_length],
        "pool_values":[0,1,2],"observable_match_checks":occurrence_checks,"rule_match_checks":rule_checks,
        "scalar_generator_identities":image_checks,"channel_joint_jump_distributions":distributions,
        "individual_molecular_successors":target_checks,"failures":failures,
        "scope":"bounded molecular-host validation after discovery; not the all-chain counting argument",
        "implementations":"MolFoundry rule-map enumeration, DPO application and observable counting"}


def compare_typed_reference(support, sources, motifs, *, gated=False):
    """Run the existing typed derivation only after chemical discovery is complete."""
    from mfobs.models.aba import counting as terminal
    if support.status != "POLYNOMIAL_SUPPORT_CLOSED":
        return {"status":"NOT_EXECUTED","reason":support.status}
    reference, ref_sources = terminal.derive(motifs,alphabet=tuple(r.symbol for r in sources if r.kind=="feed"),gated=gated)
    if reference.status != "POLYNOMIAL_SUPPORT_CLOSED":
        return {"status":"UNKNOWN_REFERENCE","reason":reference.status}
    meaning = lambda p:(p.kind,p.metadata.get("word",""))
    reference_names = {meaning(p):n for n,p in reference.patterns.items()}
    correspondence = {n:reference_names[meaning(p)] for n,p in support.patterns.items() if meaning(p) in reference_names}
    if len(correspondence) != len(support.patterns) or set(correspondence.values()) != set(reference.patterns):
        return {"status":"UNKNOWN_COORDINATE_CORRESPONDENCE","mapping":correspondence}
    subs = {sp.Symbol(n):sp.Symbol(other) for n,other in correspondence.items()}
    checked = 0
    for source in sources:
        if sp.expand(support.activities[source.identifier].xreplace(subs)-reference.activities[source.identifier]) != 0:
            raise AssertionError("molecular and typed activities disagree")
        for name,other in correspondence.items():
            if sp.expand(support.images[name][source.identifier].xreplace(subs)-reference.images[other][source.identifier]) != 0:
                raise AssertionError("molecular and typed generator images disagree")
            checked += 1
    molecular_model = polynomial_jump_model(support,sources,name="molecular-validation",input_definition_sha256="validation")
    reference_model = terminal_events.polynomial_jump_model(reference,ref_sources,name="typed-validation",input_definition_sha256="validation")
    order = tuple(reference.patterns)
    def joint_table(model, rename):
        table = {}
        for event in model.events:
            changes = {rename[n]:d for n,d in zip(model.ring.variables,event.delta)}
            key = (event.rule, tuple(changes[n] for n in order))
            polynomial = sp.Add(*[sp.Rational(t["coefficient"])*sp.prod(sp.Symbol(rename[n])**k for n,k in t["monomial"].items())
                for t in model.ring.record(event.propensity)["terms"]])
            table[key] = sp.expand(table.get(key,0)+polynomial)
        return table
    actual = joint_table(molecular_model,correspondence)
    expected = joint_table(reference_model,{n:n for n in order})
    if any(sp.expand(actual.get(k,0)-expected.get(k,0)) != 0 for k in set(actual)|set(expected)):
        raise AssertionError("molecular and typed joint event rates disagree")
    return {"status":"PASS","coordinate_correspondence":correspondence,"rate_channel_generator_identities":checked,
        "joint_jump_rate_identities":len(set(actual)|set(expected)),
        "scope":"exact symbolic activities, coordinate images and full joint event rates under the graph-count correspondence",
        "reference_role":"validation only; typed rules were not chemical composition operands"}


from mfobs.core.conditions import matches, transitions

from mfobs.core.polynomial_jump_closure import (
    generator_action,
)
from mfobs.core.recursive_words import WordState
from mfobs.models.aba.projection import TerminalWordGraphModel
from mfobs.models.molecular_rop.inputs import RopConfiguration, configuration, observable_rules, rules, state_graph
from mfobs.checks.rop_projection import composition_compilation, grammar


def graph_semantics_validation(
    selection: str | RopConfiguration,
    maximum_word_length: int = 3,
) -> dict:
    """Finite regression oracle for the literal graph encoding, never a proof input."""

    selected = configuration(selection)
    language = grammar(selected)
    observed = observable_rules(language)
    checked_states = observable_checks = rule_checks = 0
    failures = []
    words = [()]
    for length in range(1, maximum_word_length + 1):
        from itertools import product
        words.extend(product(selected.alphabet, repeat=length))
    for word in words:
        for pools in ((0,) * len(selected.alphabet), (1,) * len(selected.alphabet), (2,) * len(selected.alphabet)):
            for gate_open in ((False, True) if selected.gated else (True,)):
                semantic = WordState(tuple(word), tuple(pools), gate_open)
                host = state_graph(selected, semantic)
                expected = language.coordinate_values(semantic)
                for observable_name, rule in observed.items():
                    actual = len(matches(rule, host))
                    if actual != expected[observable_name]:
                        failures.append({
                            "kind": "observable", "state": repr(semantic),
                            "observable": observable_name, "expected": expected[observable_name],
                            "actual": actual,
                        })
                    observable_checks += 1
                for rule_name, rule in rules(selected).items():
                    actual = len(matches(rule, host))
                    if rule_name.startswith("feed_"):
                        target = 1
                    elif rule_name.startswith("insert_"):
                        symbol = rule_name[-1]
                        target = pools[selected.alphabet.index(symbol)]
                        if selected.gated and not gate_open:
                            target = 0
                    elif rule_name == "donor_bind":
                        target = int(gate_open)
                    else:
                        target = int(not gate_open)
                    if actual != target:
                        failures.append({
                            "kind": "rule", "state": repr(semantic),
                            "rule": rule_name, "expected": target, "actual": actual,
                        })
                    rule_checks += 1
                checked_states += 1
    return {
        "status": "PASS" if not failures else "FAIL",
        "scope": "finite graph regression only; the state-language induction is the all-size argument",
        "input_definition_sha256": language.input_definition_sha256,
        "maximum_word_length": maximum_word_length,
        "states": checked_states,
        "observable_match_checks": observable_checks,
        "rule_match_checks": rule_checks,
        "failures": failures,
        "gml_exported": True,
        "finite_enumeration_used_by_analysis": False,
    }


def composition_projection_validation(
    selection: str | RopConfiguration,
    maximum_word_length: int = 2,
) -> dict:
    """Finite independent regression of the structural DPO projection.

    This checks graph matches and graph targets directly.  It validates the
    implementation but is never used to construct the all-state certificate,
    the event partition, a jump vector, or an equation.
    """

    from itertools import product

    selected = configuration(selection)
    language = grammar(selected)
    graph_model = TerminalWordGraphModel(language)
    classes = graph_model.reaction_classes()
    structural_observables = graph_model.observable_rules()

    compiled = composition_compilation(selected).model
    events = {event.identifier: event for event in compiled.events}
    molecular_rules = rules(selected)
    words = [()]
    for length in range(1, maximum_word_length + 1):
        words.extend(product(selected.alphabet, repeat=length))

    states = observable_checks = class_checks = target_checks = source_checks = 0
    failures = []

    def isomorphic(left: mf.Graph, right: mf.Graph) -> bool:
        return (
            left.numVertices == right.numVertices
            and left.numEdges == right.numEdges
            and bool(matches(left, right))
        )

    for word in words:
        for pools in product(range(3), repeat=len(selected.alphabet)):
            for gate_open in ((False, True) if selected.gated else (True,)):
                state = WordState(tuple(word), tuple(pools), gate_open)
                structural_host = graph_model.state_graph(state)
                molecular_host = state_graph(selected, state)
                coordinate_values = language.coordinate_values(state)
                for observable_name, rule in structural_observables.items():
                    actual = len(matches(rule, structural_host))
                    expected = coordinate_values[observable_name]
                    if actual != expected:
                        failures.append({
                            "kind": "structural-observable",
                            "state": repr(state),
                            "observable": observable_name,
                            "expected": expected,
                            "actual": actual,
                        })
                    observable_checks += 1

                source_totals = Counter()
                for reaction in classes:
                    rule = reaction.rule
                    actual = len(matches(rule, structural_host))
                    expected = int(compiled.ring.evaluate(
                        events[reaction.identifier].propensity,
                        coordinate_values,
                    ))
                    if actual != expected:
                        failures.append({
                            "kind": "structural-activity",
                            "state": repr(state),
                            "event": reaction.identifier,
                            "expected": expected,
                            "actual": actual,
                        })
                    source_totals[reaction.source_rule] += actual
                    targets = transitions(rule, structural_host)
                    if sum(targets.values()) != actual:
                        failures.append({
                            "kind": "structural-target-multiplicity",
                            "state": repr(state),
                            "event": reaction.identifier,
                            "expected": actual,
                            "actual": sum(targets.values()),
                        })
                    if actual:
                        metadata = reaction.metadata
                        if metadata["kind"] == "feed":
                            pool_values = list(state.pools)
                            pool_values[selected.alphabet.index(metadata["symbol"])] += 1
                            successor = WordState(state.word, tuple(pool_values), state.gate_open)
                        elif metadata["kind"] == "append":
                            pool_values = list(state.pools)
                            pool_values[selected.alphabet.index(metadata["symbol"])] -= 1
                            successor = WordState(
                                state.word + (metadata["symbol"],),
                                tuple(pool_values),
                                state.gate_open,
                            )
                        else:
                            successor = WordState(
                                state.word, state.pools, not state.gate_open
                            )
                        expected_target = graph_model.state_graph(successor)
                        for target in targets:
                            if not isomorphic(target, expected_target):
                                failures.append({
                                    "kind": "structural-target",
                                    "state": repr(state),
                                    "event": reaction.identifier,
                                })
                            target_checks += 1
                    class_checks += 1

                for source_rule, structural_count in source_totals.items():
                    rule = molecular_rules[source_rule]
                    molecular_count = len(matches(rule, molecular_host))
                    if structural_count != molecular_count:
                        failures.append({
                            "kind": "chemical-projection-multiplicity",
                            "state": repr(state),
                            "source_rule": source_rule,
                            "molecular": molecular_count,
                            "structural": structural_count,
                        })
                    source_checks += 1
                states += 1
    return {
        "status": "PASS" if not failures else "FAIL",
        "scope": (
            "finite independent graph regression only; no enumerated value is "
            "used by composition, event generation, or either closure analysis"
        ),
        "input_definition_sha256": language.input_definition_sha256,
        "maximum_word_length": maximum_word_length,
        "states": states,
        "observable_match_checks": observable_checks,
        "structural_class_match_checks": class_checks,
        "structural_target_checks": target_checks,
        "molecular_to_structural_multiplicity_checks": source_checks,
        "failures": failures,
        "finite_enumeration_used_by_analysis": False,
    }


def polynomial_operator_validation(
    selection: str | RopConfiguration,
    maximum_word_length: int = 4,
) -> dict:
    """Independent finite regression of event propensities, jumps, and moments."""

    selected = configuration(selection)
    language = grammar(selected)
    model = language.polynomial_jump_model()
    words = [()]
    from itertools import product
    for length in range(1, maximum_word_length + 1):
        words.extend(product(selected.alphabet, repeat=length))
    states = 0
    event_checks = first_moment_checks = 0
    failures = []
    for word in words:
        for pools in product(range(3), repeat=len(selected.alphabet)):
            for gate_open in ((False, True) if selected.gated else (True,)):
                state = WordState(tuple(word), tuple(pools), gate_open)
                coordinates = language.coordinate_values(state)
                direct = [0] * len(model.ring.variables)
                source_class = language.classify(state.word)
                for event in model.events:
                    actual_propensity = model.ring.evaluate(event.propensity, coordinates)
                    if event.identifier.startswith("feed-"):
                        expected_propensity = 1
                    elif event.identifier.startswith("append-"):
                        pieces = event.identifier.split("-from-")
                        symbol = pieces[0].split("-")[1]
                        source = pieces[1]
                        expected_class = "".join(source_class) or "empty"
                        expected_propensity = (
                            state.pools[selected.alphabet.index(symbol)]
                            if source == expected_class and (state.gate_open or not selected.gated)
                            else 0
                        )
                    elif event.identifier == "donor-bind":
                        expected_propensity = int(state.gate_open)
                    else:
                        expected_propensity = int(not state.gate_open)
                    if actual_propensity != expected_propensity:
                        failures.append({
                            "kind": "propensity", "state": repr(state),
                            "event": event.identifier, "expected": expected_propensity,
                            "actual": str(actual_propensity),
                        })
                    for index, jump in enumerate(event.delta):
                        direct[index] += expected_propensity * jump
                    event_checks += 1
                for index, coordinate in enumerate(model.ring.variables):
                    seed = model.ring.coordinate(coordinate)
                    generated = model.ring.evaluate(
                        generator_action(model, seed),
                        coordinates,
                    )
                    if generated != direct[index]:
                        failures.append({
                            "kind": "first-moment", "state": repr(state),
                            "observable": coordinate, "expected": direct[index],
                            "actual": str(generated),
                        })
                    first_moment_checks += 1
                states += 1
    return {
        "status": "PASS" if not failures else "FAIL",
        "scope": "finite regression only; construction uses the all-word suffix proof",
        "input_definition_sha256": language.input_definition_sha256,
        "maximum_word_length": maximum_word_length,
        "states": states,
        "event_propensity_checks": event_checks,
        "first_moment_checks": first_moment_checks,
        "failures": failures,
        "finite_enumeration_used_by_analysis": False,
    }


def validate_result(record):
    """Run molecular, word and typed comparisons outside the analysis pipeline."""
    from mfobs.checks.terminal_paths import result_components, validate_words
    from mfobs.models.molecular_rop.events import compile_uniform_events
    selected, support, sources, model = result_components(record, molecular=True)
    if model is None:
        return {key: {"status": "NOT_EXECUTED", "reason": support.status} for key in
                ("word_oracle_comparison", "molecular_graph_validation", "typed_reference_comparison")}
    motifs = tuple("".join(w) for w in selected.motifs)
    return {
        "word_oracle_comparison": validate_words(support, sources, jump_model=model,
                                                  event_compiler=compile_uniform_events),
        "molecular_graph_validation": validate_molecular_hosts(selected, support, sources, jump_model=model,
            extra_words=tuple(w+w[1:] for w in motifs)+tuple(w+w for w in motifs)),
        "typed_reference_comparison": compare_typed_reference(support, sources,
            tuple(record["support_discovery"]["patterns"][n]["word"]
                  for n in record["support_discovery"]["selected_queries"]), gated=selected.gated),
    }
