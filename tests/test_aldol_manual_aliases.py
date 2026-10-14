from hashlib import sha256
import json
import molfoundry as mf
from examples.tcs import table1

def test_verified_manual_aliases_match_actual_overlap_products():
    audit=json.loads((table1.HERE/'article-manual-aliases.json').read_text(encoding='utf-8'))
    reaction=mf.Rule.fromGMLString((table1.HERE/'aldol_plain.gml').read_text(encoding='utf-8'),add=False)
    query=mf.Rule.fromGMLString((table1.HERE/'query_unconditioned.gml').read_text(encoding='utf-8'),add=False)
    products={tuple(sorted(c.overlap)):c.rule.canonKey for c in table1._candidates(reaction,query)}
    aliases=table1._alias_map()
    for row in audit['exact_overlap_evidence']:
        key=products[tuple(sorted(tuple(pair) for pair in row['reaction_to_query_overlap']))]
        assert key==row['canonical_key']
        assert aliases[sha256(key.encode()).hexdigest()]==audit['mapping'][row['article_name']]

def test_alias_merge_leaves_only_extra_count_and_retains_transcription():
    original=table1.PUBLISHED.read_bytes()
    result=table1.compute_article_table()
    assert table1.PUBLISHED.read_bytes()==original
    comparison=result['published_comparison']
    assert comparison['alias_mapping']=={'rA':'r73','rB':'r120','rC':'r68'}
    assert comparison['discrepancies']==[{'rule':'new-f7c8786fa6b1','field':'r_plus_after_r_id','published':0,'computed':1}]
    rows={r['rule']:r for r in comparison['normalized_published_rows']}
    for alias,name in comparison['alias_mapping'].items():
        assert rows[name]['published_names']==[name,alias]
        assert rows[name]['r_plus_after_r_id']==rows[name]['r_id_after_r_plus']==1
    assert result['totals']['computed']=={'r_plus_after_r_id':13,'r_id_after_r_plus':18}
