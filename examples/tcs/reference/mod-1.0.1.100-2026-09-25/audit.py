"""Replay Daniel's exact script, preserving rule allocation and exporting all runs."""
from mod import *
from collections import Counter
from hashlib import sha256
from pathlib import Path
import json
import shutil

ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / 'calc-original.py'
records = []


def handleResult(r1, r2, res, **kwargs):
    original_handleResult(r1, r2, res, **kwargs)
    encoded = r1 is aldolAdd_F_id_noEdge
    boost = bool(config.rc.useBoostCommonSubgraph)
    name = ('encoded-noedge' if encoded else 'plain') + ('-boost' if boost else '-complete')
    folder = ROOT / name
    folder.mkdir(exist_ok=True)
    inputs = {}
    for label, rule in [('identity', r1), ('reaction', r2)]:
        filename = label + '.gml'
        (folder / filename).write_text(rule.getGMLString(), encoding='utf-8')
        inputs[label] = filename
    by_id = {}
    orders = {}
    for direction, (good, bad) in zip(['before', 'after'], res):
        orders[direction] = {'accepted': len(good), 'rejected': len(bad)}
        counts = Counter(r.id for r in good)
        for r in good:
            row = by_id.setdefault(r.id, {'id': r.id, 'name': str(r), 'before': 0, 'after': 0, 'gml': 'r' + str(r.id) + '.gml'})
            row[direction] = counts[r.id]
            (folder / row['gml']).write_text(r.getGMLString(), encoding='utf-8')
    records.append({'name': name, 'useBoostCommonSubgraph': boost, 'colon_encoding': encoded,
                    'inputs': inputs, 'orders': orders, 'rows': list(by_id.values())})


source = SOURCE.read_bytes()
(ROOT / 'calc-original.py').write_bytes(source)
code = source.decode().replace('def handleResult(', 'def original_handleResult(')
assert code != source.decode()
exec(compile(code, str(SOURCE), 'exec'), globals())

# Load reference rules only after the exact replay, so their IDs cannot affect it.
reference = ROOT.parents[5] / 'evidence/tcs-table-recomputation-2026-09-21/gml/computed-table/rows'
refs = {p.stem: Rule.fromGMLString(p.read_text(encoding='utf-8'), add=False) for p in reference.glob('*.gml')}
assert refs
for record in records:
    if record['colon_encoding']:
        continue
    for row in record['rows']:
        r = Rule.fromGMLString((ROOT / record['name'] / row['gml']).read_text(encoding='utf-8'), add=False)
        row['isomorphic_molfoundry_classes'] = [name for name, ref in refs.items() if r.isomorphism(ref) > 0]
        assert len(row['isomorphic_molfoundry_classes']) == 1, row

report = {'version': str(version()), 'source_path': str(SOURCE), 'source_sha256': sha256(source).hexdigest(),
          'method': 'Original script with result-handler wrapper for export only; original match printing and evaluator order preserved. Plain rule equality checked with MOD Rule.isomorphism after all runs.',
          'records': records}
(ROOT / 'report.json').write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')
print('AUDIT_EXPORT=' + str(ROOT / 'report.json'))
