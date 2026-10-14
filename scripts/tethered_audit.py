#!/usr/bin/env python3
"""Build and audit the balanced intramolecular modification of the NAC test.

The adapter receives molecular rules, one graph seed, and an explicit finite
component language with UNBOUNDED copy numbers. Expected equations below are
post-synthesis checks only; they are never input to discovery.
"""
from __future__ import annotations

from fractions import Fraction
from itertools import product
import json
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
from mfobs.core.graph import edge_rows, node_rows
from mfobs.core.rules import make_rule

from mfobs.core.conditions import transitions
from mfobs.core.gml import graph_gml, rule_gml
from mfobs.core.worklist import run
from mfobs.inputs import load
from mfobs.models import tethered_chemistry as chemistry
from mfobs.outputs import ode_text
from mfobs.core.molecules import formula


def write_json(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, sort_keys=True)+'\n', encoding='utf-8')


def build_inputs(directory):
    directory.mkdir(parents=True, exist_ok=True)
    hosts = chemistry.language_hosts()
    seed = chemistry.seed_graph()
    for name, rule in [('guarded', chemistry.guarded_rule()),
                       ('ablated', chemistry.ablated_rule()),
                       ('alkene', make_rule('oriented-alkene', seed, seed))]:
        (directory/f'{name}.gml').write_text(rule_gml(rule), encoding='utf-8')
    for name, host in hosts.items():
        (directory/f'{name}.gml').write_text(graph_gml(host), encoding='utf-8')
    language = {'kind': 'finite-component-multiset',
        'components': {name: {'nodes': list(node_rows(host)), 'edges': list(edge_rows(host))}
                       for name, host in hosts.items()}}
    write_json(directory/'state-language.json', language)
    write_json(directory/'reduction.json',
        {'method': 'exact-connected-components-and-dynamic-span'})
    paths = {}
    for variant in ['guarded', 'ablated']:
        manifest = {'schema_version': 1, 'id': 'tethered_'+variant,
            'title': 'Atom-balanced tethered cycloaddition: '+variant,
            'adapter': 'mfobs.models.finite_components:Adapter',
            'rules': [{'id': 'cyclize', 'file': variant+'.gml', 'rate': 'k',
                       'per_match': '1/2', 'kind': 'chemical'}],
            'observables': [{'id': 'B', 'file': 'alkene.gml'}],
            'initial_observables': ['B'],
            'composition': {'mode': 'public', 'connected': False, 'include_empty': True},
            'forbidden': None, 'state_language': 'state-language.json',
            'reduction': 'reduction.json', 'normalisation': {},
            'initial_state': {'populations': {name: 2 if name in ['A', 'B'] else 0
                                             for name in hosts}},
            'parameters': None, 'reference': None}
        path = directory/f'model-{variant}.json'
        write_json(path, manifest)
        paths[variant] = path
    return paths


def verify_derivatives(result, rule):
    from mfobs.checks.observable_expressions import parse_expressions, count
    patterns = parse_expressions(result)
    checks = []
    # Small mixtures validate counts and derivatives directly on every atom map.
    # Large arbitrary populations are covered by the connected-component proof,
    # not by manufacturing enormous test graphs or extrapolating these samples.
    populations = [dict(A=a, B=b, P=p, D=0) for a, b, p in product(range(3), range(3), range(2))]
    populations += [dict(A=1, B=2, P=1, D=2), dict(A=3, B=1, P=2, D=1)]
    for population in populations:
        host = chemistry.host(population)
        targets = transitions(rule, host)
        values = {name: count(expression, host)
                  for name, expression in patterns.items()}
        for name, expression in patterns.items():
            formula_value = sum(Fraction(c)*(1 if species == '__ONE__' else population[species])
                                for species, c in result['basis_counting_functions'][name].items())
            assert values[name] == formula_value
            direct = sum(Fraction(raw, 2)*(count(expression, target)-values[name])
                         for target, raw in targets.items())
            row = result['equations_by_reaction'][name]['cyclize']
            predicted = sum(Fraction(c)*values[target] for target, c in row.items())
            assert direct == predicted, (population, name, direct, predicted)
            checks.append({'population': population, 'observable': name,
                           'value': str(values[name]), 'unit_rate_derivative': str(direct)})
    return checks


def verify_physical_equations(result, variant):
    """Independent target comparison AFTER composition-based discovery."""
    import sympy as sp
    from mfobs.outputs import matrix
    names = result['basis']
    variables = {name: sp.Symbol(name) for name in chemistry.language_hosts()}
    forms = [sum(sp.Rational(c)*(1 if s == '__ONE__' else variables[s])
                 for s, c in result['basis_counting_functions'][name].items()) for name in names]
    lhs = matrix(result)*sp.Matrix(forms)
    A, B, k = variables['A'], variables['B'], sp.Symbol('k')
    for form, derivative in zip(forms, lhs):
        independent = -k*A*sp.diff(form, A)
        independent += k*A*sp.diff(form, variables['P'])
        if variant == 'ablated':
            independent += -k*B*sp.diff(form, B)+k*B*sp.diff(form, variables['D'])
        assert sp.expand(derivative-independent) == 0
    expected_dimension = 2 if variant == 'guarded' else 1
    assert len(names) == expected_dimension


def main():
    out, examples = ROOT/'out/tethered-nac', ROOT/'examples/tethered-nac'
    out.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    paths = build_inputs(examples)
    outcomes = {}
    for variant, path in paths.items():
        result = run(load(path), max_expansions=8, progress=True)
        write_json(out/f'{variant}-result.json', result)
        if not result['closed']:
            raise RuntimeError(f'{variant}: UNKNOWN, no closed-ODE claim: {result["unknown"]}')
        rule = chemistry.guarded_rule() if variant == 'guarded' else chemistry.ablated_rule()
        checks = verify_derivatives(result, rule)
        verify_physical_equations(result, variant)
        write_json(out/f'{variant}-verification.json', {'passed': True, 'checks': checks,
            'independent_physical_generator_matches': True})
        (out/f'{variant}-ode.txt').write_text(ode_text(result), encoding='utf-8')
        outcomes[variant] = {'closed': result['closed'], 'status': result['status'],
            'basis': result['basis'], 'basis_counting_functions': result['basis_counting_functions'],
            'matrix': result['matrix'], 'seconds': result['seconds'],
            'raw_compositions': result['raw_compositions'],
            'derivative_checks': len(checks), 'input_definition_sha256': result['input_definition_sha256']}
    report = {'passed': True, 'outcomes': outcomes, 'seconds': time.perf_counter()-started,
        'hosts': {name: {'formula': formula(host), 'atoms': host.numVertices,
                         'smiles': host.smiles}
                  for name, host in chemistry.language_hosts().items()},
        'scope': 'Arbitrary independent finite copy numbers of four supplied finite molecular components; no finite batch cap.',
        'chemical_scope': 'Selected crossed intramolecular [2+2] net channel, not a complete or calibrated photochemical mechanism.',
        'not_a_fragger_separation': True,
        'noedge_active_but_positively_compilable': True,
        'physical_guarded_equations': {'expected_total_alkene_bonds': '-2*k*E[A]', 'expected_reactive_precursors': '-k*E[A]'},
        'physical_ablated_equation': {'expected_total_alkene_bonds': '-k*E[total_alkene_bonds]'}}
    write_json(out/'summary.json', report)
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
