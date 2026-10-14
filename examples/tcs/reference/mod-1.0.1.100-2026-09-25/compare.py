"""Check the captured MOD outputs against fresh MolFoundry compositions."""
from collections import Counter
from hashlib import sha256
from pathlib import Path
import json
import shutil

from mfobs.core.engine import mf
from mfobs.core.graph import edge_rows, make_graph, node_labels
from mfobs.core.rules import added_edges, make_rule, rule_graph
from examples.tcs import table1
from examples.tcs.compare_conditions import compare_conditions

ROOT = Path(__file__).resolve().parent
CAPTURE = ROOT / 'out'
REPORT = json.loads((CAPTURE / 'report.json').read_text(encoding='utf-8'))
DEST = ROOT / 'evidence'
DEST.mkdir(exist_ok=True)
for record in REPORT['records']:
    shutil.copytree(CAPTURE / record['name'], DEST / record['name'], dirs_exist_ok=True)
for name in ['calc-original.py', 'report.json']:
    shutil.copy2(CAPTURE / name, DEST / name)
for name in ['audit.py', 'compare.py']:
    shutil.copy2(ROOT / name, DEST / name)

fresh = compare_conditions(DEST / 'molfoundry')
alias_map = table1._alias_map()
checks = []
for mod_name, mf_name in [('plain-complete', 'plain'), ('encoded-noedge-complete', 'reaction_and_query_noedge')]:
    record = next(r for r in REPORT['records'] if r['name'] == mod_name)
    mf_record = fresh['variants'][mf_name]
    mf_by_key = {row['key_sha256']: row for row in mf_record['rows']}
    counters = {order: Counter() for order in ['before', 'after']}
    mappings = []
    for row in record['rows']:
        original = mf.Rule.fromGMLString((DEST / mod_name / row['gml']).read_text(encoding='utf-8'), add=False)
        rule = original
        colon_preconditions = []
        if record['colon_encoding']:
            # Only left/context sentinel edges become input absence conditions.
            # The captured example never creates a right-only sentinel edge.
            assert not [e for e in added_edges(original) if e[2] == ':']
            left, right = rule_graph(original, 'left'), rule_graph(original, 'right')
            colon_preconditions = [(u, v) for u, v, label in edge_rows(left) if label == ':']
            strip = lambda graph: make_graph(node_labels(graph), [e for e in edge_rows(graph) if e[2] != ':'])
            rule = make_rule('translated_' + str(row['id']), strip(left), strip(right),
                             tuple(sorted(set(colon_preconditions))))
            gml = rule.getGMLString()
            target = DEST / 'translated-noedge' / row['gml']
            target.parent.mkdir(exist_ok=True)
            target.write_text(gml, encoding='utf-8')
        key = rule.canonKey
        digest = sha256(key.encode()).hexdigest()
        match = mf_by_key.get(digest)
        mappings.append({**row, 'translated_gml': ('translated-noedge/' + row['gml']) if record['colon_encoding'] else None,
                         'colon_preconditions': colon_preconditions, 'molfoundry_key_sha256': digest,
                         'molfoundry_plain_class': match['plain_class'] if match else None,
                         'molfoundry_gml': 'molfoundry/' + match['gml'] if match else None,
                         'isomorphic_after_explicit_translation': match is not None})
        for order in counters:
            counters[order][digest] += row[order]
    comparisons = {}
    for order in counters:
        expected = Counter({r['key_sha256']: r[order] for r in mf_record['rows'] if r[order]})
        actual = +counters[order]
        comparisons[order] = {'equal': actual == expected, 'mod_total': sum(actual.values()),
                              'molfoundry_total': sum(expected.values()),
                              'extra': dict(actual - expected), 'missing': dict(expected - actual)}
    checks.append({'mod_run': mod_name, 'molfoundry_variant': mf_name, 'comparisons': comparisons, 'mappings': mappings})
summary = {'mod_version': REPORT['version'], 'source_sha256': REPORT['source_sha256'],
           'scope': 'Plain compositions and, for this example, explicit conversion of left/context colon sentinel edges into NoEdge input conditions. Comparison is of resulting conditional rules and per-order multiplicities; not a general correctness theorem for complement encodings.',
           'checks': checks}
(DEST / 'cross-tool-comparison.json').write_text(json.dumps(summary, indent=2) + '\n', encoding='utf-8')
for check in checks:
    print(check['mod_run'], json.dumps(check['comparisons']), 'matched classes', sum(r['isomorphic_after_explicit_translation'] for r in check['mappings']), '/', len(check['mappings']))
    for row in check['mappings']:
        if row['id'] == 978:
            print('r978 =', row['molfoundry_plain_class'])
