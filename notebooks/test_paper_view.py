"""Display/numerical interface tests; synthetic fixtures are NOT paper evidence."""
from mfobs.numerics.events import affine_generator
from pathlib import Path

import nbformat
import numpy as np
import pytest
import paper_view as v


def poly(*terms):
    return {"terms": [{"coefficient": str(c), "monomial": powers} for c, powers in terms]}


def fixture_record(nonlinear=False):
    # A deliberately simple synthetic immigration/death model for helper tests.
    return {"linear": {"status": "CLOSED"}, "pjc": {
        "status": "CLOSED", "observable_coordinates": ["X"], "operator": {"terms": [
            {"event_class": "birth", "rate": "b", "delta": {"X": 1}, "propensity_polynomial": poly((1, {}))},
            {"event_class": "death", "rate": "d", "delta": {"X": -1}, "propensity_polynomial": poly((1, {"X": 2 if nonlinear else 1}))},
        ]}}}


def test_affine_generator_uses_exported_terms():
    names, A, b = affine_generator(fixture_record(), {"b": 2., "d": .5})
    assert names == ("X",)
    assert np.array_equal(A, [[-.5]])
    assert np.array_equal(b, [2.])


def test_nonlinear_drift_is_not_a_closed_expectation_ode():
    with pytest.raises(ValueError, match="nonlinear"):
        affine_generator(fixture_record(True), {"b": 2., "d": .5})


def test_incomplete_events_are_not_simulated():
    r = fixture_record()
    r["pjc"]["status"] = "UNKNOWN"
    with pytest.raises(ValueError, match="No certified"):
        v.event_system(r)


def test_unknown_coordinates_and_noninteger_jumps_rejected():
    for delta in ({"Q": 1}, {"X": .5}):
        r = fixture_record()
        r["pjc"]["operator"]["terms"][0]["delta"] = delta
        with pytest.raises(ValueError):
            v.event_system(r)


def test_polynomial_text_is_never_executed(tmp_path):
    r = fixture_record()
    r["pjc"]["operator"]["terms"][0]["propensity_polynomial"]["text"] = "raise Exception('not executed')"
    assert v.event_system(r).activities[0] == 1


def test_parameter_validation():
    for parameters in ({"b": 1.}, {"b": 1., "d": -1.}, {"b": 1., "d": np.inf}):
        with pytest.raises(ValueError):
            affine_generator(fixture_record(), parameters)


def test_ssa_zero_rate_absorbing_state():
    times = np.linspace(0, 2, 11)
    states, trace = v.simulate_events(fixture_record(), {"X": 0}, {"b": 0., "d": 1.}, times, np.random.default_rng(4))
    assert not states.any()
    assert trace.empty
    assert trace.attrs["event_count"] == 0


def test_ssa_total_event_count_is_not_truncated_with_the_trace():
    states, trace = v.simulate_events(
        fixture_record(), {"X": 0}, {"b": 100., "d": 0.},
        [0., 1.], np.random.default_rng(4),
    )
    # Pure birth: each firing adds exactly one X, even after the trace is full.
    assert len(trace) == 25
    assert trace.attrs["event_count"] == states[-1, 0] > len(trace)


def test_ssa_repeatable_seed_and_nonnegative_counts():
    args = (fixture_record(), {"X": 3}, {"b": 2., "d": .8}, np.linspace(0, 2, 21))
    a, _ = v.simulate_events(*args, np.random.default_rng(4))
    b, _ = v.simulate_events(*args, np.random.default_rng(4))
    assert np.array_equal(a, b)
    assert (a >= 0).all()


def test_ssa_budget_never_returns_truncated_path():
    with pytest.raises(RuntimeError, match="budget"):
        v.simulate_events(fixture_record(), {"X": 0}, {"b": 100., "d": 0.}, [0, 10], np.random.default_rng(1), max_events=1)


RULE = '''rule [ ruleID "noedge" context [
node [ id 0 label "C" ] node [ id 1 label "O" ]
] right [ edge [ source 0 target 1 label "-" ] ]
constrainNoEdge [ source 0 target 1 ]
constrainAdj [ id 0 op "=" count 0 nodeLabels [ label "C" ] edgeLabels [ label "-" ] ]
]'''


def test_rule_input_interface_output_and_conditions():
    data = v.load_rule(RULE)
    assert [(edge.leftLabel, edge.rightLabel) for edge in data.edges] == [(None, "-")]
    assert all(vertex.membership.value == "context" for vertex in data.vertices)
    assert [(c.source, c.target) for c in data.matchConstraintsNoEdge] == [(0, 1)]
    assert data.matchConstraintsAdjacency[0].count == 0


def test_renderer_preserves_condition_blocks(monkeypatch):
    outputs = []
    monkeypatch.setattr(v, "display", lambda obj, **kwargs: outputs.append(
        obj['image/svg+xml'] if isinstance(obj, dict) else obj.data))
    v.show_molfoundry_rule(RULE)
    text = ''.join(outputs)
    assert 'NoEdge(0,1)' in text and 'constrainAdj' in text


def test_rule_views_read_the_public_span_and_nested_condition():
    rule = v.load_rule('''rule [
      left [ node [ id 0 label "C" ] ]
      right [ node [ id 0 label "N" ] ]
      condition [ not [ exists [ node [ id 7 label "O" ]
        edge [ source 0 target 7 label "-" ]
        not [ exists [ node [ id 8 label "H" ]
          edge [ source 7 target 8 label "-" ] ] ] ] ] ]
    ]''')
    vertex = rule.getVertex(0)
    assert vertex.leftLabel == 'C' and vertex.rightLabel == 'N'
    assert vertex.membership.value == 'changed'
    text = v.rule_conditions_html(rule)
    assert 'Nested application condition' in text
    assert text.count('not [') == text.count('exists [') == 2


def test_created_edge_admissibility_is_distinguished_from_input_provenance():
    rule = v.load_rule('''rule [ context [
      node [ id 0 label "C" ] node [ id 1 label "O" ] ]
      right [ edge [ source 0 target 1 label "-" ] ] ]''')
    assert rule.matchConstraintsNoEdge[0].implied
    assert 'Created-edge admissibility' in v.rule_conditions_html(rule)
    assert 'may also be written explicitly' in v.rule_conditions_html(rule)


def test_notebook_schema_and_all_python_cells_compile():
    import re
    root = Path(__file__).parent
    names = ['paper_examples', *[p.relative_to(Path(__file__).parent).with_suffix('').as_posix() for p in (Path(__file__).parent/'examples').glob('[0-9][0-9]_*.ipynb')]]
    for name in names:
        notebook = nbformat.read(root / (name + '.ipynb'), as_version=4)
        nbformat.validate(notebook)
        script = (root / (name + '.py')).read_text(encoding='utf-8')
        blocks = re.findall(r'^# %%\n(.*?)(?=^# %%|\Z)', script, flags=re.M | re.S)
        cells = [cell for cell in notebook.cells if cell.cell_type == 'code']
        assert len(blocks) == len(cells)
        for index, (cell, block) in enumerate(zip(cells, blocks)):
            compile(cell.source, f'{name}-cell-{index}', 'exec')
            assert block.strip() == cell.source.strip()
        assert '<nav' not in script and 'style="background:' not in script
        if name != 'paper_examples':
            assert 'paper_examples.ipynb' in notebook.cells[0].source
            assert any(c.get('metadata', {}).get('jupyter', {}).get('source_hidden') for c in notebook.cells)
    index = nbformat.read(root / 'paper_examples.ipynb', as_version=4)
    assert all(cell.cell_type == 'markdown' for cell in index.cells)
    assert '](visual/CORE.md)' in index.cells[0].source
    for chapter in sorted((root/'examples').glob('[0-9][0-9]_*.ipynb')):
        assert f'](examples/{chapter.name})' in index.cells[0].source


def test_negative_rate_rejected_at_runtime():
    r = fixture_record()
    r['pjc']['operator']['terms'][0]['propensity_polynomial'] = poly((-1, {}))
    with pytest.raises(ValueError, match='Invalid rates'):
        v.simulate_events(r, {'X': 1}, {'b': 1., 'd': 1.}, [0, 1], np.random.default_rng(0))


def test_named_rule_channels_map_to_rate_parameters():
    r = fixture_record()
    r['contract'] = {'normalized_rules': [{'id': 'feed', 'rate': 'b'}]}
    r['pjc']['operator']['terms'][0]['rate'] = 'feed'
    assert v.event_system(r).channels[0] == 'b'


def test_word_coordinates_preserve_counts_and_validate_inputs():
    r = fixture_record()
    r['pjc']['observable_coordinates'] = ['M_A', 'N_A']
    r['pjc']['operator']['terms'] = [{
        'event_class': 'feed', 'rate': 'b', 'delta': {'M_A': 1},
        'propensity_polynomial': poly((1, {}))}]
    r['support_discovery'] = {'patterns': {
        'M_A': {'kind': 'pool', 'word': 'A'},
        'N_A': {'kind': 'motif', 'word': 'A'}}}
    assert v.word_coordinates(r, 'AAA', {'A': 2}) == {'M_A': 2, 'N_A': 3}
    for word, pools, gate in [('A', {'A': 1.5}, True),
                             ('AB', {'A': 2}, True),
                             ('A', {'A': 2}, 7),
                             ('A', {'A': True}, True)]:
        with pytest.raises(ValueError):
            v.word_coordinates(r, word, pools, gate)


def test_total_rate_overflow_is_rejected():
    r = fixture_record()
    with np.errstate(over='ignore'):
        with pytest.raises(ValueError, match='Total event rate'):
            v.simulate_events(r, {'X': 1}, {'b': 1e308, 'd': 1e308},
                              [0, 1], np.random.default_rng(2))


def test_native_aldol_drawing_and_conditions(monkeypatch):
    import xml.etree.ElementTree as ET
    from urllib.parse import quote
    gml = (Path(__file__).parents[1] / 'examples/tcs/aldol.gml').read_text(encoding='utf-8')
    outputs = []
    monkeypatch.setattr(v, 'display', lambda obj, **kwargs: outputs.append(
        obj['image/svg+xml'] if isinstance(obj, dict) else obj.data))
    v.show_molfoundry_rule(gml, '<Aldol>')
    assert '&lt;Aldol&gt;' in outputs[0]
    assert ET.fromstring(outputs[1]).tag.endswith('svg')
    assert 'NoEdge(0,5)' in outputs[2] and 'NoEdge(3,4)' in outputs[2]
    assert 'Created-edge admissibility' in outputs[2]
    assert quote(gml, safe='') in outputs[2]
    assert not any('<script' in text or '<iframe' in text for text in outputs)


def test_native_drawing_keeps_nested_condition_visible(monkeypatch):
    from mfobs.models.chlorination.inputs import RULE_GML
    outputs = []
    monkeypatch.setattr(v, 'display', lambda obj, **kwargs: outputs.append(
        obj['image/svg+xml'] if isinstance(obj, dict) else obj.data))
    v.show_molfoundry_rule(RULE_GML)
    assert '<svg' in outputs[1]
    assert 'Nested application condition' in outputs[2]
    assert 'exists' in outputs[2]
    assert 'https://molfoundry.de/rule-editor?gml=' in outputs[2]


def test_single_panel_requires_identity():
    from mfobs.rule_drawing import rule_svg
    identity = v.load_rule('rule [ context [ node [ id 0 label "C" ] ] ]')
    assert '<svg' in rule_svg(identity, identity=True)
    deletion = v.load_rule('rule [ left [ node [ id 0 label "C" ] ] ]')
    with pytest.raises(ValueError, match='identity rule'):
        rule_svg(deletion, identity=True)


def test_unknown_operator_overview_is_not_zero_events():
    row = v.overview({'linear': {'status': 'UNKNOWN'}, 'pjc': {'status': 'UNKNOWN', 'operator': None}}).iloc[0]
    assert row['Event rows'] is None and row['Coordinates'] is None


def test_selected_equations_reject_unknown_coordinates(monkeypatch):
    outputs = []
    monkeypatch.setattr(v, 'display', lambda obj, **kwargs: outputs.append(
        obj['image/svg+xml'] if isinstance(obj, dict) else obj.data))
    v.show_equations(fixture_record(), names=('X',))
    assert len(outputs) == 2
    with pytest.raises(ValueError, match='Unknown requested'):
        v.show_equations(fixture_record(), names=('missing',))


def test_native_drawings_show_named_edge_types():
    import xml.etree.ElementTree as ET
    from mfobs.rule_drawing import rule_svg
    from mfobs.rule_drawing import path_rule_svg

    gml = '''rule [ context [
      node [ id 0 label "Position" ] node [ id 1 label "Next" ]
      node [ id 2 label "Position" ] node [ id 3 label "Symbol:A" ]
      edge [ source 0 target 1 label "src" ]
      edge [ source 1 target 2 label "tgt" ]
      edge [ source 2 target 3 label "kind" ]
    ] ]'''
    from mfobs.visual_links import compact_rule_gml
    compact = compact_rule_gml(gml)
    drawings = (rule_svg(v.load_rule(compact)),
                rule_svg(v.load_rule(compact), identity=True),
                path_rule_svg(gml, identity=True))
    for drawing in drawings:
        labels = {text.text for text in ET.fromstring(drawing).iter()
                  if text.tag.endswith('}text')}
        assert {'src', 'tgt', 'kind'} <= labels


def test_native_chemical_bonds_keep_their_strokes():
    import re
    import molfoundry as mf
    from mfobs.rule_drawing import rule_svg

    rule = v.load_rule('''rule [ context [ node [ id 0 label "C" ]
        node [ id 1 label "O" ] edge [ source 0 target 1 label "=" ] ] ]''')
    printer = mf.GraphPrinter()
    printer.darkMode = mf.DarkMode.Off
    printer.withIndex = True
    printer.collapseHydrogens = False
    printer.drawRuleConditions = True
    before = rule.svg(printer)
    after = rule_svg(rule)
    # MolFoundry assigns a fresh SVG/CSS id on every drawing.
    normalize = lambda svg: re.sub(r'mf-[a-z0-9]+', 'mf-render', svg)
    assert normalize(after) == normalize(before)


def test_graph_drawing_preserves_named_edges_and_source_ids(monkeypatch):
    import xml.etree.ElementTree as ET
    outputs = []
    monkeypatch.setattr(v, 'display', lambda obj, **kwargs: outputs.append(
        obj['image/svg+xml'] if isinstance(obj, dict) else obj.data))
    v.show_graph('''graph [ node [ id 9 label "P" ] node [ id 33 label "T" ]
        edge [ source 9 target 33 label "src" ] ]''', '<Host>')
    assert '&lt;Host&gt;' in outputs[0]
    labels = {node.text for node in ET.fromstring(outputs[1]).iter()
              if node.tag.endswith('}text')}
    assert {'P:9', 'T:33', 'src'} <= labels


def test_notebook_molecule_uses_native_graph_and_source_ids(monkeypatch):
    import xml.etree.ElementTree as ET
    outputs = []
    monkeypatch.setattr(v, 'display', lambda obj, **kwargs: outputs.append(
        obj['image/svg+xml'] if isinstance(obj, dict) else obj.data))
    v.show_molecule('''graph [ node [ id 9 label "C" ] node [ id 33 label "O" ]
        edge [ source 9 target 33 label "=" ] ]''', indices=True)
    labels = {''.join(node.itertext()) for node in ET.fromstring(outputs[0]).iter()
              if node.tag.endswith('}text')}
    assert {'C:9', 'O:33'} <= labels


def test_notebook_svg_survives_stylesheet_removal(monkeypatch):
    import xml.etree.ElementTree as ET
    gml = (Path(__file__).parents[1] / 'examples/tcs/aldol.gml').read_text(encoding='utf-8')
    before = v.rule_svg(v.load_rule(gml))
    outputs = []
    monkeypatch.setattr(v, 'display', lambda obj, **kwargs: outputs.append(
        obj['image/svg+xml'] if isinstance(obj, dict) else obj.data))
    v.show_molfoundry_rule(gml)
    assert outputs[1].isascii() and '&#8743;' in outputs[1] and '&#8708;' in outputs[1]
    root = ET.fromstring(outputs[1])
    assert not any(node.tag.endswith('}style') for node in root.iter())
    backgrounds = [node for node in root.iter() if 'bg' in node.get('class', '').split()]
    assert backgrounds and all(node.get('fill') == '#FFF' for node in backgrounds)
    labels = [node for node in root.iter() if node.tag.endswith('}text')]
    assert any(node.text == 'O' and node.get('fill') == '#000' for node in labels)
    assert all(node.get('font-family') == 'Calibri, Helvetica, sans-serif' for node in labels)
    # Presentation attributes must not change the geometry, labels or conditions.
    original = ET.fromstring(before)
    original_labels = [node.text for node in original.iter() if node.tag.endswith('}text')]
    assert [node.text for node in labels] == original_labels
    assert root.get('viewBox') == original.get('viewBox')
    assert v.notebook_svg(outputs[1]) == outputs[1]


def test_notebook_svg_keeps_nested_graph_palettes_separate():
    import xml.etree.ElementTree as ET
    drawing = '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 40 20">' + ''.join(
        f'<svg id="{name}"><style>#{name} .bg {{ fill: {colour}; }}</style>'
        '<rect class="bg" width="20" height="20"/></svg>'
        for name, colour in [('first', '#FFF'), ('second', '#000')]) + '</svg>'
    root = ET.fromstring(v.notebook_svg(drawing))
    assert [node.get('fill') for node in root.iter() if node.tag.endswith('}rect')] == ['#FFF', '#000']
