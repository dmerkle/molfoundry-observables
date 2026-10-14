"""Conditional tables stand alone and link the exact native rules."""
from pathlib import Path
from html import unescape
from urllib.parse import parse_qs, urlsplit
import json
import re
import sys


HERE=Path(__file__).resolve().parent
sys.path[:0]=[str(HERE),str(HERE.parent/'src'),str(HERE.parent)]
import inspection


def test_native_noedge_table_visible_with_exact_links(tmp_path,monkeypatch):
    run=type('Run',(),{'metadata':{'source_sha256':{},'versions':{}}})()
    monkeypatch.setattr(inspection,'FILES',tmp_path/'inspection')
    outputs=[]
    monkeypatch.setattr(inspection,'display',lambda value:outputs.append(value.data))
    plain=tmp_path/'plain/computed-table';plain.mkdir(parents=True)
    (plain/'rows').mkdir()
    (plain/'rows/r73.gml').write_text('rule [ context [ node [ id 0 label "C" ] ] ]', encoding='utf-8')
    (plain/'table.json').write_text(json.dumps({'rows':[{'rule':'r73','gml':'rows/r73.gml','r_plus_after_r_id':1,'r_id_after_r_plus':1}],
        'totals':{'computed':{'r_plus_after_r_id':1,'r_id_after_r_plus':1}},
        'published_comparison':{'alias_mapping':{'rA':'r73','rB':'r120','rC':'r68'}}}), encoding='utf-8')
    conditional=tmp_path/'conditional';conditional.mkdir()
    gml='rule [ context [ node [ id 0 label "C" ] node [ id 1 label "C" ] ] constrainNoEdge [ source 0 target 1 ] ]'
    (conditional/'conditioned.gml').write_text(gml, encoding='utf-8')
    second=gml.replace('label "C"','label "O"')
    (conditional/'second.gml').write_text(second, encoding='utf-8')
    row={'plain_class':'r73','key_sha256':'a'*64,'gml':'conditioned.gml','after':1,'before':1,'net_coefficient':0,
         'constraints':[{'source':0,'target':1}]}
    rows=[{**row,'key_sha256':'b'*64,'gml':'second.gml','after':0,'before':2,'net_coefficient':-2},row]
    comparison={'input_sha256':{},'variants':{key:{'rows':rows} for key in ('reaction_noedge','reaction_and_query_noedge')}}
    (conditional/'report.json').write_text(json.dumps(comparison), encoding='utf-8')
    directory=inspection.aldol_tables(comparison,conditional,tmp_path/'plain',run)
    visible=[html for html in outputs if html.startswith('<h3>With NoEdge')]
    folded=[html for html in outputs if html.startswith('<details><summary>With NoEdge')]
    assert len(visible)==len(folded)==1
    assert '<details>' not in visible[0]
    assert 'reaction and query · 2 classes' in visible[0]
    assert 'reaction only · 2 classes' in folded[0]
    cases=[(visible[0],Path.cwd()),(folded[0],Path.cwd())]
    cases += [((directory/(variant+'.html')).read_text(encoding='utf-8'),directory)
              for variant in ('reaction_noedge','reaction_and_query_noedge')]
    for html,base in cases:
        assert re.findall(r'<th>(.*?)</th>',html)==['Rule (GML)','After','Before','After − before']
        for clutter in ['MØD','translated','encoding','Native NoEdge pairs','Underlying plain class',
                        'Conditional class hash','both','r73','(0, 1)','a'*12]:
            assert clutter not in html
        assert 'coefficient zero' in html
        links=re.findall(r'<a href="([^"]+)"[^>]*>([^<]+)</a>',html)
        home = [(url, label) for url, label in links if label == '← All examples']
        if home:
            assert (base / home[0][0]).resolve() == inspection.PACKAGE / 'notebooks/paper_examples.ipynb'
        links = [(url, label) for url, label in links if label != '← All examples']
        assert [label for _,label in links]==['1','2']
        assert [urlsplit(unescape(url)).netloc for url,_ in links]==['molfoundry.de']*2
        assert [parse_qs(urlsplit(unescape(url)).query)['gml'][0] for url,_ in links]==[gml,second]
        assert '</a></td><td>1</td><td>1</td><td>0</td></tr>' in html
        assert '<td><b>Total</b></td><td>1</td><td>3</td><td>-2</td>' in html
    assert (directory/'conditional/conditioned.gml').read_text(encoding='utf-8')==gml
    assert json.loads((directory/'conditional/report.json').read_text(encoding='utf-8'))==comparison
    assert not any('Literal published transcription' in html for html in outputs)
