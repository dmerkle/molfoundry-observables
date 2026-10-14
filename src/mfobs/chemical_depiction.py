"""Molecule views: MolFoundry by default; RDKit for exact bond highlights.

The explicit-atom conversion is also used by independent chemical tests.
No matching, rule application or observable analysis uses RDKit.
"""
import xml.etree.ElementTree as ET
import molfoundry as mf
from mfobs.core.gml import component_graphs
from mfobs.core.graph import vertex_ids

ET.register_namespace('', 'http://www.w3.org/2000/svg')

def to_rdkit(graph):
    """Read MolFoundry's explicit mapped SMILES for highlights and chemical checks."""
    from rdkit import Chem
    if not graph.isMolecule:
        raise ValueError('Chemical atom and bond labels are required')
    parameters = Chem.SmilesParserParams()
    parameters.removeHs = False
    parameters.sanitize = False
    mapped = '.'.join(part.smilesWithExternalIds for part in component_graphs(graph))
    molecule = Chem.MolFromSmiles(mapped, parameters)
    if molecule is None or molecule.GetNumAtoms() != graph.numVertices:
        raise ValueError('Mapped SMILES did not preserve the explicit source atoms')
    for atom in molecule.GetAtoms():
        atom.SetIntProp('source_atom_id', atom.GetAtomMapNum())
        atom.SetAtomMapNum(0)
        atom.SetNoImplicit(True)
        atom.SetNumRadicalElectrons(0)
    Chem.SanitizeMol(molecule, sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_FINDRADICALS)
    return Chem.RenumberAtoms(molecule, sorted(range(molecule.GetNumAtoms()),
                             key=lambda i: molecule.GetAtomWithIdx(i).GetIntProp('source_atom_id')))


def molecular_svg(graph, title='', *, width=520, height=260, highlights=(),
                  highlight_edges=(), labels=None, indices=False,
                  preserve_hydrogens=(), kekulize=False):
    """Draw a complete molecule, retaining exactly the requested mapped bonds."""
    selected = set(highlights)
    labels = labels or {}
    # The native printer currently marks charged aromatic carbon as a radical.
    charged = {v.id for v in graph.vertices if v.stringLabel.endswith(('+', '-'))}
    charged_aromatic = any(e.stringLabel == ':' and
                           (e.source.id in charged or e.target.id in charged)
                           for e in graph.edges)
    if not highlight_edges and not preserve_hydrogens and not charged_aromatic:
        printer = mf.GraphPrinter()
        printer.darkMode = mf.DarkMode.Off
        printer.collapseHydrogens = True
        printer.withExternalId = indices
        source_ids = vertex_ids(graph)
        def overlay(vertex):
            sid = source_ids[vertex.id]
            if sid in selected or sid in labels:
                return mf.VertexOverlay(
                    color='#b8e6e8' if sid in selected else None,
                    label=labels.get(sid, vertex.stringLabel if sid in selected else None),
                    shape=mf.OverlayShape.Circle if sid in selected else mf.OverlayShape.none)
        if selected or labels:
            printer.pushVertexOverlay(overlay)
        drawing = ET.fromstring(graph.svg(printer))
        drawing.set('width', str(width))
        drawing.set('height', str(height - 24 if title else height))
        root = ET.Element('{http://www.w3.org/2000/svg}svg',
                          viewBox=f'0 0 {width} {height}', width=str(width), height=str(height))
        root.append(drawing)
        if title:
            legend = ET.SubElement(root, '{http://www.w3.org/2000/svg}text', {
                'x': str(width/2), 'y': str(height-5), 'text-anchor': 'middle',
                'font-family': 'sans-serif', 'font-size': '14'})
            legend.text = title
    else:
        from rdkit import Chem
        from rdkit.Chem.Draw import rdMolDraw2D
        full = to_rdkit(graph)
        for atom in full.GetAtoms():
            if atom.GetIntProp('source_atom_id') in preserve_hydrogens:
                atom.SetAtomMapNum(atom.GetIdx()+1)
        parameters = Chem.RemoveHsParameters()
        parameters.removeMapped = False
        molecule = Chem.RemoveHs(full, parameters)
        displayed = {a.GetIntProp('source_atom_id'): a.GetIdx() for a in molecule.GetAtoms()}
        for atom in molecule.GetAtoms():
            atom.SetAtomMapNum(0)
            sid = atom.GetIntProp('source_atom_id')
            if sid in labels or indices:
                atom.SetProp('atomNote', str(labels.get(sid, sid)))
        bonds = []
        for u, v, *_ in highlight_edges:
            if u in displayed and v in displayed:
                bond = molecule.GetBondBetweenAtoms(displayed[u], displayed[v])
                if bond is None:
                    raise ValueError('Highlighted pattern bond is absent from the host')
                bonds.append(bond.GetIdx())
        molecule = rdMolDraw2D.PrepareMolForDrawing(molecule, kekulize=kekulize)
        drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
        options = drawer.drawOptions()
        options.prepareMolsBeforeDrawing = False
        options.padding = .12
        options.bondLineWidth = 1.6
        options.legendFontSize = 15
        options.minFontSize = 13
        atoms = [displayed[u] for u in selected if u in displayed]
        # Always pass exact bonds: RDKit otherwise colours edges between selected atoms.
        drawer.DrawMolecule(molecule, legend=title, highlightAtoms=atoms, highlightBonds=bonds,
                            highlightAtomColors={i: (.72, .90, .91) for i in atoms},
                            highlightBondColors={i: (.72, .90, .91) for i in bonds})
        drawer.FinishDrawing()
        root = ET.fromstring(drawer.GetDrawingText())
    heading = ET.Element('{http://www.w3.org/2000/svg}title')
    heading.text = title
    root.insert(0, heading)
    return ET.tostring(root, encoding='unicode')
