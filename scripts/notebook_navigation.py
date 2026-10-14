#!/usr/bin/env python3
"""Update section links in the paper notebooks without changing executable cells."""
from pathlib import Path
import json
import re
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
from html import escape

WALKTHROUGH_SECTIONS = (
    ("rules", "1 · Rules & state space"),
    ("inputs", "2 · Inputs & semantics"),
    ("observables", "3 · Observable discovery"),
    ("ode", "4 · ODE & pJC"),
    ("initial", "5 · Initial states"),
    ("dynamics", "6 · Dynamics"),
)


WALKTHROUGH_TABS = WALKTHROUGH_SECTIONS


# The same topics and colours identify existing content on pages and notebooks.
TOPIC_COLORS = ("#e7f3e9", "#fbe9e9", "#e8f0fb", "#f0eafb",
                "#e3f3f0", "#fff0de")


def closure_heading(closure):
    """Name only the results that this presentation actually exposes."""
    names = [name for name, key in (("ODE", "ode"), ("pJC", "pjc")) if closure[key][0]]
    return "4 · " + (" & ".join(names) if names else "Closure status")


def section_navigation(targets=None, *, closure=None, notes=None, tabbed=False):
    """Six topics, with separate links for exact ODE and pJC results.

    Targets are (href, location/reason). A missing href is disabled, with its
    reason visible. Cross-surface links state their destination. Inline styling
    also works in saved notebooks. Availability describes displayed results,
    not a claim that an unresolved closure is mathematically impossible.
    """
    targets = targets or {key: ('#' + key, '') for key, _ in WALKTHROUGH_TABS}
    closure = closure or {key: targets['ode'] for key in ('ode', 'pjc')}
    notes = notes or {}
    items = []
    for (key, label), color in zip(WALKTHROUGH_TABS, TOPIC_COLORS):
        href, note = targets[key]
        available = any(t[0] for t in closure.values()) if key == 'ode' else bool(href)
        style = (f'background:{color if available else "#f2f3f4"};color:#243d48;'
                 'padding:9px 12px;border-radius:4px;display:inline-block;'
                 'font-size:13px;font-weight:600;line-height:1.4;white-space:nowrap')

        def item(text, target, *, result=None):
            url, reason = target
            if not url:
                return (f'<span aria-disabled="true" title="{escape(reason, quote=True)}" '
                        'style="color:#737d83;display:inline-block;vertical-align:top">'
                        f'{escape(text)}<small style="display:block;font-size:10px;font-weight:400">'
                        f'{escape(reason)}</small></span>')
            extra = ''
            # Existing tab controllers select the enclosing panel; the shared
            # section script then scrolls to the requested result within it.
            if tabbed and result and url.startswith('#') and url != '#ode':
                extra = f' data-scroll-target="{escape(url[1:], quote=True)}"'
                url = '#ode'
            suffix = f' ↗ {reason}' if reason in ('Web', 'Notebook') else ''
            detail = reason if reason and not suffix else notes.get(key, '')
            return (f'<a href="{escape(url, quote=True)}"{extra} '
                    'style="color:#243d48;padding:0;display:inline-block;font-size:13px;'
                    'font-weight:600;text-decoration:underline;border:0;background:transparent">'
                    f'{escape(text + suffix)}'
                    + (f'<small style="display:block;font-size:10px;font-weight:400">{escape(detail)}</small>' if detail else '')
                    + '</a>')

        if key == 'ode':
            content = '4 · ' + item('ODE', closure['ode'], result='ode') + ' &amp; ' + item('pJC', closure['pjc'], result='pjc')
            if not available and href:
                content += '<small style="display:block">' + item('Status', targets[key]) + '</small>'
        else:
            content = item(label, (href, note))
        identifier = f' id="tab-{key}"' if tabbed else ''
        items.append(f'<span data-topic="{key}"{identifier} style="{style}">{content}</span>')
    kind = 'tabs' if tabbed else 'example-sections'
    return (f'<nav class="{kind}" aria-label="Example sections" '
            'style="display:flex;align-items:stretch;gap:5px;flex-wrap:wrap;margin:16px auto;'
            'padding-top:8px;padding-bottom:8px;border-bottom:1px solid #cddce1">'
            + ''.join(items) + '</nav>')


# Short notebooks link to their existing cells or to the visual presentation.
SHORT = {
    '01_aldol': {'rules':'651f6501', 'inputs':'72e9e137'},
    '02_aba': {'inputs':'bc0784e4', 'observables':'0f534a89', 'ode':'0f534a89',
               'initial':'0de11a31', 'dynamics':'171c073c'},
    '03_chlorination': {'rules':'nac-rules-query', 'inputs':'nac-rules-query',
                       'observables':'b2c2f115', 'ode':'0c2e5456', 'initial':'430082e7',
                       'dynamics':'d38c9bf8'},
    '04_remy': {'observables':'89b5b8ab', 'ode':'89b5b8ab', 'initial':'8c6523a5',
                'dynamics':'8b509420'},
}
WEB = {
    ('02_aba','rules'): 'rop-family/index.html#aba-rules',
    ('04_remy','rules'): 'remy-trees/index.html#rules',
    ('04_remy','inputs'): 'remy-trees/index.html#inputs',
}
OLD_IDS = ('1-Model','2-Query','3-Symbolic-analysis','4-Equations',
           '5-Initial-state','6-Dynamics')


# Availability here is presentation-specific: external results link to their
# existing presentation rather than adding duplicate event tables to notebooks.
CLOSURE = {
    '01_aldol': {'ode': (None, 'Not analysed'), 'pjc': (None, 'Not analysed')},
    '02_aba': {'ode': (None, 'Unresolved'), 'pjc': ('../../experiments/rop-family/index.html#aba-pjc', 'Web')},
    '03_chlorination': {'ode': (None, 'Not derived'), 'pjc': ('#ode', '')},
    '04_remy': {'ode': ('#ode', ''), 'pjc': ('../../experiments/remy-trees/index.html#events', 'Web')},
    '05_eight_slots': {'ode': ('#ode', ''), 'pjc': ('../../experiments/eight-slots/index.html#pjc', 'Web')},
    '06_dna_walker': {'ode': ('#ode', ''), 'pjc': (None, 'Not analysed')},
    '07_molecular_rop': {'ode': (None, 'ABA: unresolved'), 'pjc': ('#ode', '')},
    '08_query_variation': {'ode': (None, 'Unresolved'), 'pjc': ('#observables', '')},
    '09_composition_comparison': {'ode': (None, 'Not analysed'), 'pjc': (None, 'Not analysed')},
    '10_brusselator': {'ode': (None, 'Unresolved'), 'pjc': ('#ode', '')},
    '11_growing_ring': {'ode': (None, 'Unresolved'), 'pjc': (None, 'Unresolved')},
}
NOTES = {'10_brusselator': {'dynamics': 'SSA & approximation'},
         '11_growing_ring': {'dynamics': 'Approximation'}}


def cell_target(cell):
    identifier = ('section-nav-' if cell['cell_type'] == 'code' else '') + cell['id']
    return '#cell-id=' + identifier


def update(path):
    book = json.loads(path.read_text(encoding='utf-8'))
    before = [c for c in book['cells'] if c['cell_type']=='code']
    book['cells'] = [c for c in book['cells'] if not c['id'].startswith('section-nav-')]
    targets, placements = {}, {}
    closure = dict(CLOSURE[path.stem])
    for i, ((key, label), color) in enumerate(zip(WALKTHROUGH_SECTIONS, TOPIC_COLORS)):
        if key == 'ode':
            local = {k: (url if url and url.startswith('#') else None, note)
                     for k, (url, note) in closure.items()}
            label = closure_heading(local)
            if path.stem == '08_query_variation':
                label = '4 · Generator identities'
        if (path.stem,key) in WEB:
            targets[key] = ('https://stochmech.algochem.techfak.de/experiments/'+WEB[path.stem,key], 'Web')
            continue
        if path.stem == '01_aldol' and key not in SHORT[path.stem]:
            targets[key] = (None, 'Composition audit only')
            continue
        if path.stem == '09_composition_comparison' and key in ('observables','ode','initial','dynamics'):
            targets[key] = (None, 'Composition measurements only')
            continue
        if path.stem == '08_query_variation' and key in ('initial','dynamics'):
            aba = json.loads((path.parent/'02_aba.ipynb').read_text(encoding="utf-8"))
            cell = next(c for c in aba['cells'] if c['id'] == SHORT['02_aba'][key])
            targets[key] = ('02_aba.ipynb'+cell_target(cell), 'Notebook')
            continue
        targets[key] = ('#'+key, '')
        if path.stem in SHORT:
            if key not in SHORT[path.stem]:
                continue
            cell_id = SHORT[path.stem][key]
            cell = next(c for c in book['cells'] if c['id'] == cell_id)
            targets[key] = (cell_target(cell), '')
            placements.setdefault(cell_id,[]).append((key,label,color))
        else:
            cell = next(c for c in book['cells'] if c['cell_type']=='markdown' and
                        (f'id="{OLD_IDS[i]}"' in ''.join(c['source']) or f'id="{key}"' in ''.join(c['source'])))
            targets[key] = (cell_target(cell), '')
            text = ''.join(cell['source'])
            text = re.sub(r'<h2[^>]*>.*?</h2>',
                f'<h2 id="{key}" style="background:{color};padding:12px 16px;border-radius:5px">{label}</h2>',text,count=1)
            cell['source'] = text.splitlines(keepends=True)
    # Remove the old coloured link strip, preserving the surrounding introduction.
    first = book['cells'][0]
    first['source'] = re.sub(r'<a href="#[1-7]-[^\n]+\n?', '', ''.join(first['source'])).splitlines(keepends=True)
    for key, (url, note) in closure.items():
        if url and url.startswith('#'):
            closure[key] = (targets[url[1:]][0], note)
        elif url and url.startswith('../../experiments/'):
            closure[key] = ('https://stochmech.algochem.techfak.de/experiments/'+url.removeprefix('../../experiments/'), note)
    navigation = section_navigation(targets, closure=closure, notes=NOTES.get(path.stem))
    navigation = re.sub(r'(<a href="https://[^"]+")',
                        r'\1 target="_blank" rel="noopener noreferrer"', navigation)
    cells = [first, {'cell_type':'markdown','id':'section-nav-links','metadata':{},
                     'source':[navigation]}]
    for cell in book['cells'][1:]:
        if cell['id'] in placements:
            topics = placements[cell['id']]
            # Existing prose is its own link target; code gets a short label.
            if cell['cell_type']=='code':
                content = '<p>' + ' · '.join(
                    f'<strong style="background:{color};padding:4px 8px;border-radius:4px">{label}</strong>'
                    for _,label,color in topics) + '</p>'
                cells.append({'cell_type':'markdown','id':'section-nav-'+cell['id'],'metadata':{},'source':[content]})
        cells.append(cell)
    book['cells']=cells
    assert before == [c for c in cells if c['cell_type']=='code'], 'Executable cells or outputs changed'
    path.write_text(json.dumps(book,indent=1,ensure_ascii=False)+'\n', encoding='utf-8')


if __name__=='__main__':
    for path in sorted((ROOT/'notebooks/examples').glob('[0-9][0-9]_*.ipynb')):
        update(path)
    print('Updated section navigation; executable cells and saved outputs unchanged.')
