"""Compare fresh MolFoundry output with the exported MØD 1.0.1.100 run.

MØD columns are recorded evidence, never input to MolFoundry composition.
Rule identities are checked again on each call; numeric totals alone do not
establish agreement. The sentinel translation is specific to these products.
"""
from collections import Counter
from hashlib import sha256
import json
from pathlib import Path

from mfobs.core.engine import mf
from mfobs.core.graph import edge_rows, make_graph, node_labels
from mfobs.core.rules import added_edges, make_rule, rule_graph

EVIDENCE = Path(__file__).resolve().parents[2] / 'examples/tcs/reference/mod-1.0.1.100-2026-09-25'


def translate_absence_edges(rule):
    """Translate this example's input/context ':' sentinels to NoEdge."""
    if any(label == ':' for _, _, label in added_edges(rule)):
        raise ValueError('Right-only absence sentinels are outside this translation')
    left, right = rule_graph(rule, 'left'), rule_graph(rule, 'right')
    pairs = {(u, v) for u, v, label in edge_rows(left) if label == ':'}
    plain = lambda side: make_graph(node_labels(side), [e for e in edge_rows(side) if e[2] != ':'])
    return make_rule('translated', plain(left), plain(right), tuple(sorted(pairs)))


def compare_recorded_mod(fresh, *, evidence=EVIDENCE):
    """Check all rule classes and multiplicities against fresh composition."""
    evidence = Path(evidence)
    recorded = json.loads((evidence / 'report.json').read_text(encoding='utf-8'))
    if sha256((evidence / 'calc-original.py').read_bytes()).hexdigest() != recorded['source_sha256']:
        raise ValueError('The recorded MØD source fingerprint changed')
    result = {'version': recorded['version'], 'evidence': str(evidence),
              'source_sha256': recorded['source_sha256'], 'runs': {}}
    for name, variant in [('plain-complete', 'plain'),
                          ('encoded-noedge-complete', 'reaction_and_query_noedge')]:
        record = next(r for r in recorded['records'] if r['name'] == name)
        actual = fresh['variants'][variant]['rows']
        by_key = {r['key_sha256']: r for r in actual}
        counts = {order: Counter() for order in ('after', 'before')}
        rows = []
        for row in record['rows']:
            path = evidence / name / row['gml']
            rule = mf.Rule.fromGMLString(path.read_text(encoding='utf-8'), add=False)
            translated_path = None
            if record['colon_encoding']:
                rule = translate_absence_edges(rule)
                translated_path = evidence / 'translated-noedge' / row['gml']
                saved = mf.Rule.fromGMLString(translated_path.read_text(encoding='utf-8'), add=False)
                if saved.canonKey != rule.canonKey:
                    raise ValueError('Saved translated GML differs from current translation')
            key = sha256(rule.canonKey.encode()).hexdigest()
            if key not in by_key:
                raise ValueError(f'MØD r{row["id"]} has no matching MolFoundry rule')
            for order in counts:
                counts[order][key] += row[order]
            rows.append({**row, 'key_sha256': key, 'plain_class': by_key[key]['plain_class'],
                         'original_gml': str(path),
                         'translated_gml': str(translated_path) if translated_path else None})
        for order in counts:
            expected = Counter({r['key_sha256']: r[order] for r in actual if r[order]})
            if +counts[order] != expected:
                raise ValueError(f'MØD/MolFoundry {variant} {order} multiplicities differ')
        result['runs'][variant] = {'status': 'PASS', 'rows': rows,
                                  'totals': {k: sum(v.values()) for k, v in counts.items()}}
    return result
