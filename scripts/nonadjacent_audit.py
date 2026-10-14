#!/usr/bin/env python3
"""Reproduce a balanced covalent NAC diagnostic and its closure obstruction.

No chemical model or target closure is fitted here. The canonical synthesis
budget returns UNKNOWN; a separate analytical induction establishes the
absence of uniform finite linear observable closure for the alkene seed.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
from mfobs.chemical_depiction import to_rdkit
from mfobs.core.graph import make_graph, node_rows
from mfobs.core.rules import make_rule

from mfobs.checks.pair_depletion import polynomial_certificate, finite_batch
from mfobs.core.conditions import matches, transitions, parse_observable

from mfobs.core.gml import graph_gml, rule_gml
from mfobs.core.engine import commutator
from mfobs.core.worklist import run
from mfobs.inputs import load
from mfobs.models import nonadjacent_chemistry as chemistry
from mfobs.core.gml import bond_order_sums
from mfobs.core.molecules import formula, smiles


def write_json(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, sort_keys=True)+'\n', encoding='utf-8')


def build_inputs(directory):
    directory.mkdir(parents=True, exist_ok=True)
    seed = make_graph({0: 'C', 1: 'C'}, [(0, 1, '=')])
    for name, rule in [('guarded', chemistry.guarded_rule()),
                       ('ablated', chemistry.ablated_rule()),
                       ('alkene', make_rule('oriented-alkene', seed, seed))]:
        (directory/f'{name}.gml').write_text(rule_gml(rule), encoding='utf-8')
    write_json(directory/'state-language.json', {
        'scope': 'All finite explicit-H C/H graphs with ordinary single/double bonds, C valence 4 and H valence 1.',
        'closed_under_rule': 'Every atom and its bond-order sum is preserved.',
        'not_assumed_by_adapter': 'The canonical adapter does not use a chemical language counter.',
        'inventory': 'Each finite closed initial mixture has a finite reachable space; initial inventories have no uniform bound.'})
    write_json(directory/'reduction.json', {'method': 'canonical jump-span equality and direct edge/NoEdge contradiction elimination',
        'no_target_basis': True, 'no_polynomial_certificate_used': True})
    manifest = {'schema_version': 1, 'id': 'nonadjacent_cycloaddition',
        'title': 'Balanced covalent pair-consuming NAC diagnostic (not validated kinetics)',
        'adapter': 'mfobs.core.identity:Adapter',
        'rules': [{'id': 'cycloaddition', 'file': 'guarded.gml', 'rate': 'k',
                   'per_match': '1/8', 'kind': 'chemical'}],
        'observables': [{'id': 'B', 'file': 'alkene.gml'}],
        'initial_observables': ['B'],
        'composition': {'mode': 'public', 'connected': False, 'include_empty': True},
        'forbidden': None, 'state_language': 'state-language.json',
        'reduction': 'reduction.json',
        'normalisation': {'B': {'scale': '2', 'meaning': 'raw oriented C=C matches divided by 2 count alkene bonds'}},
        'initial_state': {}, 'parameters': None, 'reference': None}
    write_json(directory/'model.json', manifest)
    return seed, load(directory/'model.json')


def audit_hosts(seed, out):
    hosts = [('ethylene_0', chemistry.ethylene_mixture(0)),
             ('ethylene_1', chemistry.ethylene()),
             ('ethylene_2', chemistry.ethylene_mixture(2)),
             ('ethylene_3', chemistry.ethylene_mixture(3)),
             ('ethylene_5', chemistry.ethylene_mixture(5)),
             ('conjugated_diene', chemistry.conjugated_diene()),
             ('separated_diene', chemistry.separated_diene())]
    records, composition_records = [], []
    for kind, rule in [('guarded', chemistry.guarded_rule()), ('ablated', chemistry.ablated_rule())]:
        engine_rule = rule
        engine_seed = make_rule('B', seed, seed)
        result = commutator(engine_seed, engine_rule, include_empty=True)
        terms = []
        for key, coef in result.observable_coefficients.items():
            parsed, adj = parse_observable(result.observable_representatives[key])
            assert parsed.jump_identity and not parsed.generic_conditions
            pattern = make_graph(parsed.labels, parsed.edges)
            terms.append((coef, pattern, parsed.no_edges, adj))
        composition_records.append({'kind': kind, 'left_raw': result.left_raw,
            'right_raw': result.right_raw, 'signed_jump_terms': len(terms),
            'terms_with_no_edge': sum(bool(ne) for _, _, ne, _ in terms)})
        for name, host in hosts:
            targets = transitions(rule, host)
            before = len(matches(seed, host))
            direct = sum(Fraction(raw, 8)*(len(matches(seed, target))-before)
                         for target, raw in targets.items())
            composed = sum(Fraction(coef, 8)*len(matches(pattern, host, ne, adj))
                           for coef, pattern, ne, adj in terms)
            assert direct == composed, (kind, name, direct, composed)
            for target in targets:
                assert node_rows(target) == node_rows(host)
                assert bond_order_sums(target) == bond_order_sums(host)
                mol = to_rdkit(target)
                assert all(a.GetFormalCharge() == 0 and a.GetNumRadicalElectrons() == 0
                           for a in mol.GetAtoms())
            record = {'host': name, 'kind': kind, 'formula': formula(host),
                'smiles': smiles(host),
                'successor_smiles': sorted({smiles(target) for target in targets}),
                'raw_admissible_maps': sum(targets.values()),
                'atom_id_distinct_successors': len(targets),
                'generator_of_raw_B': str(direct), 'composition_agrees': True}
            records.append(record)
            (out/f'{name}.gml').write_text(graph_gml(host), encoding='utf-8')
    return records, composition_records


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--max-expansions', type=int, default=1)
    args = parser.parse_args()
    out, examples = ROOT/'out/nonadjacent-pairs', ROOT/'examples/nonadjacent-pairs'
    out.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    seed, inputs = build_inputs(examples)
    records, compositions = audit_hosts(seed, out)
    worklist_result = run(inputs, max_expansions=args.max_expansions)
    write_json(out/'worklist.json', worklist_result)
    certificate = polynomial_certificate(8)
    batches = [finite_batch(n0) for n0 in [2, 4, 8, 16, 32, 33, 64]]
    report = {'status': 'CHECKED_BALANCED_NAC_DIAGNOSTIC_WITH_ANALYTICAL_NONCLOSURE',
        'checks_passed': True, 'closed': False, 'worklist_status': worklist_result['status'],
        'zero_jump_terms_discarded': len(worklist_result.get('zero_jump_terms', {})),
        'not_a_new_finite_closure': True, 'not_a_validated_chemical_kinetic_model': True,
        'no_fragger_impossibility_claim': True,
        'input_definition_sha256': inputs.definition_digest(),
        'rule_gml_sha256': hashlib.sha256((examples/'guarded.gml').read_bytes()).hexdigest(),
        'elapsed_seconds': time.perf_counter()-started,
        'molecular_checks': records, 'composition_checks': compositions,
        'polynomial_induction_checks': certificate, 'finite_batches': batches,
        'proof_scope': 'No finite-dimensional Q-invariant linear function space containing alkene count, uniform over arbitrary ethylene/cyclobutane inventories, for k>0. This excludes exact autonomous moment closure for arbitrary initial laws, not mean-field approximation or finite batch enumeration.',
        'important_caveat': 'NoEdge is active on conjugated diene hosts, but redundant on the ethylene depletion sublanguage used for the nonclosure proof. Nonclosure is not caused by NACs.',
        'chemistry_caveat': 'NoEdge selects a nonconjugated-pair model channel; it does not assert that conjugated diene photocyclization is chemically impossible.'}
    write_json(out/'audit.json', report)
    print(json.dumps({'status': report['status'], 'checks_passed': True,
                      'worklist_status': worklist_result['status'],
                      'molecular_generator_checks': len(records),
                      'finite_batch_dimensions': [(r['initial_ethylene'], r['seed_krylov_dimension']) for r in batches],
                      'elapsed_seconds': report['elapsed_seconds']}, indent=2))


if __name__ == '__main__':
    main()
