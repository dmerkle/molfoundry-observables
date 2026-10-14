"""MolFoundry links preserve native graph semantics and compact display labels."""
from html import unescape
import re
from urllib.parse import parse_qs, urlsplit


def test_notebook_graph_file_link_uses_depiction(tmp_path):
    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).parents[1] / 'notebooks'))
    from inspection import href
    path = tmp_path / 'host.gml'
    path.write_text('graph [ node [ id 0 label "Position" ] ]', encoding="utf-8")
    url = unescape(re.search(r'href="([^"]+)"', href(path, 'Host')).group(1))
    assert urlsplit(url).path == '/depiction'
    assert parse_qs(urlsplit(url).query)['gml'] == ['graph [ node [ id 0 label "P" ] ]']


def test_compact_display_retains_nested_conditions_and_adjacency_labels():
    from mfobs.visual_links import compact_rule_gml, editor_url as display_url
    gml = '''rule [ context [ node [ id 0 label "Position" ] ]
        constrainAdj [ id 0 op "=" count 0
            nodeLabels [ label "Next" label "Terminal" ]
            edgeLabels [ label "src" label "tgt" ] ]
        condition [ not [ exists [ node [ id 1 label "Next" ]
            edge [ source 0 target 1 label "src" ]
            not [ exists [ node [ id 2 label "Terminal" ]
                edge [ source 1 target 2 label "tgt" ] ] ] ] ] ] ]'''
    url = urlsplit(display_url(gml))
    query = parse_qs(url.query)
    assert url.path == '/depiction'
    assert query['withEdgeLabels'] == ['on']
    assert query['drawRuleConditions'] == ['on']
    compact = query['gml'][0]
    assert compact == compact_rule_gml(gml)
    restored = compact
    for short, full in [('P', 'Position'), ('N', 'Next'), ('T', 'Terminal')]:
        restored = restored.replace(f'label "{short}"', f'label "{full}"')
    assert restored == gml


def test_short_path_rules_keep_concrete_match_counts(aba_records):
    import json
    from pathlib import Path
    import molfoundry as mf
    from mfobs.visual_links import compact_rule_gml
    from mfobs.models.rop_exploration import expand
    for case in ('aba', 'aba-gated'):
        result = aba_records[case]
        expansion = expand(result, depth=1)
        rules = result['support_discovery']['source_rules']
        for record in rules.values():
            original = mf.Rule.fromGMLString(record['gml'], add=False)
            shown = mf.Rule.fromGMLString(compact_rule_gml(record['gml']), add=False)
            for state in expansion['states'][:3]:
                host = mf.Graph.fromGMLString(state['gml'], add=False)
                short_host = mf.Graph.fromGMLString(compact_rule_gml(state['gml']), add=False)
                assert original.numMatches([host], mf.ApplyMode.UseAny) == shown.numMatches([short_host], mf.ApplyMode.UseAny)
