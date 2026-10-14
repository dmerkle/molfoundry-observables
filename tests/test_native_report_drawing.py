"""Literal reports retain mapped atoms, labels and conditions in native drawings."""
from mfobs.core.graph import make_graph
import xml.etree.ElementTree as ET


def test_plain_molecules_and_vertex_markers_do_not_require_rdkit(monkeypatch):
    import builtins
    import molfoundry as mf
    from mfobs.chemical_depiction import molecular_svg
    original = builtins.__import__
    def without_rdkit(name, *args, **kwargs):
        assert not name.startswith('rdkit'), 'Ordinary views must use MolFoundry'
        return original(name, *args, **kwargs)
    graph = mf.smiles('CC', add=False)
    carbons = [v.id for v in graph.vertices if v.stringLabel == 'C']
    monkeypatch.setattr(builtins, '__import__', without_rdkit)
    assert '<svg' in molecular_svg(graph)
    drawing = molecular_svg(graph, highlights=carbons)
    markers = [n for n in ET.fromstring(drawing).iter()
               if n.tag.endswith('}circle') and n.get('fill') == '#b8e6e8']
    assert len(markers) == 2


def test_bond_highlight_retains_only_the_selected_pattern_edge():
    from mfobs.chemical_depiction import molecular_svg
    from mfobs.core.molecules import complete_with_h
    graph = complete_with_h(make_graph({0: 'C', 1: 'C', 2: 'C'},
                                      [(0, 1, '-'), (1, 2, '-')]))
    drawing = molecular_svg(graph, highlights=(0, 1, 2), highlight_edges=((0, 1),))
    highlighted = {frozenset(n.get('class').split()[1:]) for n in ET.fromstring(drawing).iter()
                   if n.get('class', '').startswith('bond-')
                   and '#B7E5E8' in n.get('style', '')}
    assert highlighted == {frozenset(('atom-0', 'atom-1'))}
