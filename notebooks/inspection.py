"""Inspectable files and source navigation for the visual examples.

Only presentation and serialization live here. Rules, coefficients, reductions,
closure and event construction are supplied by the existing analysis modules.
"""
from __future__ import annotations
from hashlib import sha256
from html import escape
import json
import os
from pathlib import Path
import shutil
from urllib.parse import quote

from IPython.display import HTML, display
from mfobs.visual_links import editor_link, editor_url
from mfobs.visual_links import MAX_GRAPH_LINK_LENGTH

HERE = Path(__file__).resolve().parent
PACKAGE = HERE.parent
FILES = HERE / 'visual/inspection'
STYLE = '''body{font:16px/1.45 system-ui,sans-serif;color:#243345;max-width:1250px;margin:30px auto;padding:0 24px}a{color:#18568a}h1,h2{line-height:1.2}h2{padding-top:20px}table{border-collapse:collapse;width:100%;font-size:14px;margin:15px 0}th,td{padding:6px 10px;border-bottom:1px solid #dce2e7;text-align:left;vertical-align:top}th{background:#edf1f5}tr:nth-child(even){background:#f8fafb}code,pre{font-family:ui-monospace,monospace}pre{white-space:pre-wrap;font-size:12px}summary{cursor:pointer}small{color:#526174}.line{display:block;white-space:pre}.line:target{background:#fff0bd}.number{display:inline-block;width:4em;color:#718096;text-decoration:none}nav a{margin-right:16px}'''


def href(path, label, *, base=None):
    if Path(path).suffix == '.gml':
        gml = Path(path).read_text(encoding='utf-8')
        if not gml.lstrip().startswith('graph') or len(editor_url(gml)) <= MAX_GRAPH_LINK_LENGTH:
            return editor_link(gml, str(label))
    relative=os.path.relpath(Path(path), Path.cwd() if base is None else Path(base))
    return f'<a href="{escape(quote(relative),quote=True)}">{escape(str(label))}</a>'


def page(path,title,body):
    Path(path).parent.mkdir(parents=True,exist_ok=True)
    Path(path).write_text('<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">'
        '<title>'+escape(title)+'</title><style>'+STYLE+'</style><body>'
        +href(PACKAGE/'notebooks/paper_examples.ipynb','← All examples',base=Path(path).parent)
        +'<h1>'+escape(title)+'</h1>'+body+'</body></html>', encoding='utf-8')


def gml_link(gml,label='View in MolFoundry'):
    return editor_link(gml, label)


def input_rules(items):
    display(HTML('<p>'+' · '.join(gml_link(gml,name) for name,gml in items.items())+'</p>'))


def html_table(headers, rows):
    # Cells are explicitly escaped by callers; only generated links carry HTML.
    return '<table><thead><tr>'+''.join('<th>'+escape(h)+'</th>' for h in headers)+'</tr></thead><tbody>'+''.join(
        '<tr>'+''.join('<td>'+str(value)+'</td>' for value in row)+'</tr>' for row in rows)+'</tbody></table>'


def _directory(example,inputs,run):
    binding={'inputs':inputs,'source':run.metadata['source_sha256'],'versions':run.metadata['versions']}
    directory=FILES/(example+'-'+sha256(json.dumps(binding,sort_keys=True).encode()).hexdigest()[:12])
    directory.mkdir(parents=True,exist_ok=True)
    (directory/'provenance.json').write_text(json.dumps({**binding,'run':run.metadata},indent=2)+'\n', encoding='utf-8')
    return directory


def _json(path,data):
    Path(path).write_text(json.dumps(data,indent=2,ensure_ascii=False)+'\n', encoding='utf-8')


def aldol_tables(comparison,composition_dir,plain_dir,run):
    """Show the published/plain table and export every conditional/raw row."""
    directory=_directory('aldol',comparison['input_sha256'],run)
    shutil.copytree(composition_dir,directory/'conditional',dirs_exist_ok=True)
    shutil.copytree(Path(plain_dir)/'computed-table',directory/'plain',dirs_exist_ok=True)
    report=json.loads((directory/'plain/table.json').read_text(encoding='utf-8'))
    published=json.loads((PACKAGE/'examples/tcs/article-table.json').read_text(encoding='utf-8'))
    _json(directory/'published-table.json',published)
    (directory/'chemistry.json').write_bytes((PACKAGE/'examples/tcs/chemistry.json').read_bytes())
    from examples.tcs.table1 import normalized_published_rows
    merged = normalized_published_rows(published)
    _json(directory/'published-normalized.json', {'rows': merged, 'alias_mapping': report['published_comparison']['alias_mapping']})
    actual={row['rule']:row for row in report['rows']}; old={row['rule']:row for row in merged}
    names=list(old)+[name for name in actual if name not in old]
    def plain_table(base):
        rows=[]
        for name in names:
            a,p=actual.get(name,{}),old.get(name,{})
            label = ' / '.join(dict.fromkeys([name, *p.get('published_names', [])]))
            rule_link=href(directory/'plain'/a['gml'],'MolFoundry GML',base=base) if a else 'not generated'
            rows.append([escape(label),rule_link,p.get('r_plus_after_r_id','—'),
                a.get('r_plus_after_r_id',0),
                p.get('r_id_after_r_plus','—'),a.get('r_id_after_r_plus',0),
                '×' if a.get('violates_C1') else '','×' if a.get('violates_C2') else ''])
        rows.append(['<b>Total</b>','',published_total(published,'r_plus_after_r_id'),
                     report['totals']['computed']['r_plus_after_r_id'],
                     published_total(published,'r_id_after_r_plus'),
                     report['totals']['computed']['r_id_after_r_plus'],'',''])
        return html_table(['Class','Rules','After · TCS','After · MolFoundry',
                          'Before · TCS','Before · MolFoundry','C1','C2'],rows)
    intro = ('<p>Plain composition, with both explicit NoEdge conditions removed. '
             'After = reaction then query; before = query then reaction. '
             'Published aliases rA/r73, rB/r120 and rC/r68 are grouped by isomorphism.</p>'
             '<p>C1 marks two double bonds at one carbon; C2 marks two oxygen neighbours. '
             'These annotations do not filter rows.</p>')
    display(HTML('<details><summary>Plain compositions · comparison with TCS Table 1</summary>'
                 +intro+plain_table(Path.cwd())+'</details>'))
    cond_bodies=[]
    for variant,title in [('reaction_and_query_noedge','With NoEdge · reaction and query'),
                          ('reaction_noedge','With NoEdge · reaction only')]:
        rows=comparison['variants'][variant]['rows']
        def cond_table(base):
            body=[]
            for number,row in enumerate(sorted(rows,key=lambda r:r['key_sha256']),1):
                link=href(directory/'conditional'/row['gml'],str(number),base=base)
                body.append([link,row['after'],row['before'],row['net_coefficient']])
            body.append(['<b>Total</b>',sum(r['after'] for r in rows),sum(r['before'] for r in rows),
                         sum(r['net_coefficient'] for r in rows)])
            return html_table(['Rule (GML)','After','Before','After − before'],body)
        target=directory/(variant+'.html')
        scope=('<p>MolFoundry conditional composition. Each numbered link opens the exact rule in MolFoundry, '
               'including its application conditions. All conditional isomorphism classes are retained, '
               'including canceled classes (coefficient zero). '
               'After = reaction then query; before = query then reaction.</p>')
        page(target,title,scope+cond_table(directory))
        if variant=='reaction_and_query_noedge':
            display(HTML('<h3>'+escape(title)+f' · {len(rows)} classes</h3>'+scope+cond_table(Path.cwd())))
        else:
            display(HTML('<details><summary>'+escape(title)+f' · {len(rows)} classes</summary>'+scope+cond_table(Path.cwd())+'</details>'))
        cond_bodies.append(href(target,title,base=directory))
    raw_rows=[]
    for variant,record in comparison['variants'].items():
        for order,rows in record.get('raw_candidates',{}).items():
            for index,row in enumerate(rows,1):
                raw_rows.append([escape(variant),escape(order),index,'kept' if row['accepted_by_valence_filter'] else 'rejected by valence filter',
                    href(directory/'conditional'/row['gml'],'GML',base=directory),escape(str(row['overlap_first_to_second']))])
    page(directory/'raw-composites.html','Aldol addition: all raw composed candidates',
         '<p>Products successfully returned by conditional composition, before the external balanced/max-valence filter. '
         'Each overlap has a row; isomorphic products may occur more than once. Null/contradictory compositions have no resulting GML rule.</p>'
         +html_table(['Input variant','Chronological order','Overlap #','Filter result','Composed rule','Overlap map'],raw_rows))
    page(directory/'index.html','Aldol addition: inspect the calculation',
         intro+plain_table(directory)
         +'<h2>Conditional tables</h2><p>'+' · '.join(cond_bodies)+'</p><p>'+
         href(directory/'raw-composites.html',f'All {len(raw_rows)} raw candidates, with accepted/rejected flags',base=directory)+'</p>')
    return directory


def published_total(report,key):
    return sum(row[key] for row in report['rows'])
