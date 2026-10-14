#!/usr/bin/env python3
"""Exhaustive independent verification AFTER non-enumerative synthesis.

Checks every assignment on each allowed topology, then demonstrates how three
all-free initial graphs supply different vectors for the SAME computed ODE.
"""
import argparse
import csv
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import sys
from time import perf_counter

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
try:
    from .check_local_slots import parse_pattern
except ImportError:  # Direct script invocation.
    from check_local_slots import parse_pattern
from mfobs.checks import eight_slots as ref
from mfobs.checks import local_slots_reference as graph_ref
from mfobs.inputs import load
from mfobs.models.eight_slots import polynomial_count


def canonical_state(bits, lengths):
    components = []; offset = 0
    for n in lengths:
        word = tuple(bits[offset:offset+n]); reverse = word[::-1]
        canonical = min([word[i:]+word[:i] for i in range(n)] +
                        [reverse[i:]+reverse[:i] for i in range(n)])
        components.append((n, canonical)); offset += n
    return tuple(sorted(components))


def numeric_trajectories(result, families, out, *, write_observables=True):
    import numpy as np
    from scipy.sparse import csr_matrix
    from scipy.sparse.linalg import expm_multiply
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    names = result['basis']; index = {name: i for i, name in enumerate(names)}
    rates = {'basal_bind': 0.08, 'cooperative_bind': 1.2, 'unbind': 0.4}
    numeric = np.zeros((len(names), len(names)))
    for source, channels in result['equations_by_reaction'].items():
        for channel, row in channels.items():
            for target, coefficient in row.items():
                numeric[index[source], index[target]] += rates[channel]*float(Fraction(coefficient))
    parameters = {'basal': rates['basal_bind'], 'cooperative': rates['cooperative_bind'], 'off': rates['unbind']}
    times = np.linspace(0, 30, 151); curves = {}; initial = {}; maximum_error = 0.0
    for tag, family in families.items():
        states = family['states']; counts = family['counts']; edges = family['edges']
        state_index = {bits: i for i, bits in enumerate(states)}
        zero = (0,)*8; vector = np.array([counts[zero][name] for name in names], dtype=float)
        initial[tag] = {'topology': list(ref.TOPOLOGIES[tag]), 'occupancy_bits': list(zero),
                        'observable_initial_values': {name: str(counts[zero][name]) for name in names}}
        ode = expm_multiply(csr_matrix(numeric), vector, start=0, stop=30, num=len(times))
        q = np.zeros((256, 256))
        for bits in states:
            i = state_index[bits]
            for target, _, rate in graph_ref.events(bits, edges, **parameters):
                q[i, state_index[target]] += float(rate); q[i, i] -= float(rate)
        p0 = np.zeros(256); p0[state_index[zero]] = 1
        probabilities = expm_multiply(csr_matrix(q.T), p0, start=0, stop=30, num=len(times))
        evaluations = np.array([[counts[bits][name] for name in names] for bits in states])
        expected = probabilities @ evaluations
        error = float(np.max(np.abs(expected-ode))); maximum_error = max(maximum_error, error)
        assert error < 1e-8, (tag, error)
        assert np.min(ode[:, index['P0']]) > -1e-9 and np.max(ode[:, index['P0']]) < 8+1e-9
        curves[tag] = ode[:, index['P0']]
    with (out/'trajectory.csv').open('w', newline='', encoding='utf-8') as handle:
        writer = csv.writer(handle); writer.writerow(['time', *[f'occupied_{tag}' for tag in families]])
        writer.writerows([float(t), *[float(curves[tag][i]) for tag in families]] for i, t in enumerate(times))
    (out/'initial-conditions.json').write_text(json.dumps(initial, indent=2)+'\n', encoding='utf-8')
    differing = {name: {tag: row['observable_initial_values'][name] for tag, row in initial.items()}
                 for name in names if len({row['observable_initial_values'][name] for row in initial.values()}) > 1}
    assert differing, 'Different topology trajectories need distinguishable initial observable vectors'
    if write_observables:
        last_round = max(row['round'] for row in result['rounds'])
        table = ['# Discovered observables and three initial vectors', '',
            'All entries evaluate the exported exact observable expressions in a completely free initial graph. '
            'Only P0 was supplied; the other coordinates were discovered by composition. '
            '`0` means free, `1` occupied, `*` unspecified. GML files are conditional strict identity rules.', '',
            '| Observable | Discovery depth | Pattern | Initial 8 | Initial 4+4 | Initial 3+5 |',
            '|---|---:|---|---:|---:|---:|']
        for name in names:
            description = result.get('identity_patterns', {}).get(name, {}).get('meaning', 'exact linear combination')
            initial_columns = ' | '.join(str(initial[tag]['observable_initial_values'][name]) for tag in families)
            link = f'observables/{name}.json' if 'observable_expressions' in result else f'identity-patterns/round-{last_round}/{name}.gml'
            table.append(f'| [{name}]({link}) | {result["depths"][name]} | {description} | {initial_columns} |')
        table += ['', '[Computed ODE](ode.txt) · [Exact verification](verification.json) · [Initial vectors](initial-conditions.json)', '',
                  '![One ODE, different topology-dependent initial vectors](topology-trajectories.svg)']
        (out/'observables.md').write_text('\n'.join(table)+'\n', encoding='utf-8')
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.1), constrained_layout=True)
    for tag, curve in curves.items():
        axes[0].plot(times, curve, label=tag)
        if tag != '8':
            axes[1].plot(times, curve-curves['8'], label=f'{tag} minus 8')
    axes[0].set(title='Same ODE, three all-free initial graphs', xlabel='time', ylabel='Expected occupied slots')
    axes[1].set(title='Topology-dependent difference', xlabel='time', ylabel='Difference in expected occupancy')
    for axis in axes:
        axis.grid(alpha=0.2); axis.legend()
    fig.savefig(out/'topology-trajectories.svg')
    fig.savefig(out/'topology-trajectories.png', dpi=150)
    plt.close(fig)
    return {'rates': rates, 'time_interval': [0, 30], 'time_points': len(times),
            'same_matrix_for_all_topologies': True, 'initially_differing_observables': differing,
            'max_ODE_vs_independent_CTMC_error_all_observables': maximum_error,
            'max_occupancy_difference_from_C8': {tag: float(np.max(np.abs(curve-curves['8'])))
                                               for tag, curve in curves.items() if tag != '8'}}


def verify(out: Path, *, result=None):
    from mfobs.provenance import source_hashes, check_runtime
    import sympy as sp
    from sympy.polys.matrices import DomainMatrix
    start = perf_counter(); out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    path = out/'result.json'
    stored_result = result is None
    if stored_result:
        result = json.loads(path.read_text(encoding='utf-8'))
    assert result['closed'] and not result['pending'] and not result['unknown']
    assert result['input_definition_sha256'] == load(ROOT/'examples/eight-slots/model.json').definition_digest()
    check_runtime(result)
    assert result['source_sha256'] == source_hashes(), 'Synthesis sources changed; recompute first'
    for source, expected in result['source_sha256'].items():
        assert sha256((ROOT/source).read_bytes()).hexdigest() == expected, (
            f'Saved synthesis used a different version of {source}. '
            'Regenerate the result (notebook MODE="compute") before verifying.'
        )
    names = result['basis']
    expressions = result.get('observable_expressions', {
        n: [{'coefficient': '1', 'gml': item['gml']}]
        for n, item in result['identity_patterns'].items()})
    parsed = {term['gml']: parse_pattern(term['gml'])
              for terms in expressions.values() for term in terms}
    families = {}; family_reports = {}; count_checks = 0; row_checks = 0; evaluation_rows = set()
    for tag, lengths in ref.TOPOLOGIES.items():
        edges = ref.topology_edges(lengths); states = list(ref.assignments(lengths)); counts = {}
        for bits in states:
            host = graph_ref.typed_slots(bits, edges)
            term_counts = {gml: graph_ref.typed_match_count(graph, host, no_edges)
                           for gml, (graph, no_edges) in parsed.items()}
            counts[bits] = {name: sum(Fraction(t['coefficient']) * term_counts[t['gml']] for t in terms)
                            for name, terms in expressions.items()}
            for name in names:
                polynomial = (result['basis_polynomials'][name] if 'basis_polynomials' in result
                              else result['symbolic_patterns'][name]['polynomial'])
                assert counts[bits][name] == polynomial_count(polynomial, bits, lengths), (tag, bits, name)
                count_checks += 1
            evaluation_rows.add(tuple(counts[bits][name] for name in names))
        channels = {'basal_bind': 'basal', 'cooperative_bind': 'cooperative', 'unbind': 'off'}
        for bits in states:
            for reaction, channel in channels.items():
                parameters = {key: int(key == channel) for key in channels.values()}
                events = list(graph_ref.events(bits, edges, **parameters))
                for name in names:
                    direct = sum(rate*(counts[target][name]-counts[bits][name]) for target, _, rate in events)
                    predicted = sum(Fraction(c)*counts[bits][target]
                                    for target, c in result['equations_by_reaction'][name][reaction].items())
                    assert direct == predicted, (tag, bits, name, reaction, direct, predicted)
                    row_checks += 1
        family_reports[tag] = {'occupancy_assignments_checked': len(states),
            'state_orbits_counted_only_after_synthesis': len({canonical_state(bits, lengths) for bits in states}),
            'patterns_checked': len(names), 'rate_separated_rows_checked': len(states)*len(names)*3}
        families[tag] = {'states': states, 'counts': counts, 'edges': edges}
        print(tag, json.dumps(family_reports[tag]), flush=True)
    # Entries are exact rationals. Use rational elimination rather than the
    # general symbolic rank routine, which is costly for residual coordinates.
    evaluation = DomainMatrix.from_Matrix(sp.Matrix(sorted(evaluation_rows)))
    rank = evaluation.rank()
    dependencies = [{name: str(value) for name, value in zip(names, vector) if value}
                    for vector in evaluation.nullspace().to_Matrix().tolist()]
    trajectories = numeric_trajectories(result, families, out, write_observables=stored_result)
    sources = [Path(__file__), Path(ref.__file__), Path(graph_ref.__file__), ROOT/'scripts/check_local_slots.py']
    report = {'passed': True, 'result_sha256': sha256(path.read_bytes()).hexdigest() if stored_result else None,
        'source_sha256': {p.relative_to(ROOT).as_posix(): sha256(p.read_bytes()).hexdigest() for p in sources},
        'synthesis_used_no_host_or_orbit_evaluations': result['synthesis_state_evaluations'] == result['synthesis_orbit_evaluations'] == 0,
        'formal_polynomial_generators': len(names), 'verified_count_function_rank': rank,
        'verified_additional_count_dependencies': dependencies,
        'families': family_reports, 'independent_raw_pattern_vs_polynomial_checks': count_checks,
        'exact_rate_separated_generator_rows': row_checks, 'trajectories': trajectories,
        'seconds': perf_counter()-start,
        'scope': 'Exhaustive on every occupancy of each compatible eight-slot topology up to graph relabeling. Inert orphan P/L vertices do not affect anchored counts or transitions. These state evaluations occur only after synthesis, never as its reducer.'}
    (out/'verification.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
    print(json.dumps({key: value for key, value in report.items() if key not in ('source_sha256', 'trajectories')}, indent=2), flush=True)

    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, default=ROOT/'out/eight-slots',
                        help='directory produced by run_eight_slots.py')
    args = parser.parse_args()
    verify(args.out.resolve())


if __name__ == '__main__':
    main()
