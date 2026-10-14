"""Explicit atom semantics of the shared RDKit depiction adapter."""
from mfobs.chemical_depiction import to_rdkit
import pytest

Chem = pytest.importorskip('rdkit.Chem')

from mfobs.core.graph import make_graph
from mfobs.models.chlorination import inputs as chlorination


def test_triple_bond_and_source_ids_survive_hydrogen_folding():
    graph = make_graph({3: 'C', 7: 'C', 11: 'H', 19: 'H'},
                       [(3, 7, '#'), (3, 11, '-'), (7, 19, '-')])
    molecule = to_rdkit(graph)
    assert [atom.GetIntProp('source_atom_id') for atom in molecule.GetAtoms()] == [3, 7, 11, 19]
    assert Chem.MolToSmiles(Chem.RemoveHs(molecule)) == 'C#C'
    folded = Chem.RemoveHs(molecule)
    assert [atom.GetIntProp('source_atom_id') for atom in folded.GetAtoms()] == [3, 7]
    assert all(atom.GetNumImplicitHs() == 0 and atom.GetNumRadicalElectrons() == 0
               for atom in molecule.GetAtoms())


@pytest.mark.parametrize('label,charge', [('N+', 1), ('N2+', 2), ('Mg2+', 2), ('C-', -1)])
def test_formal_charge_does_not_infer_missing_hydrogens_or_radicals(label, charge):
    molecule = to_rdkit(make_graph({12: label}))
    atom = molecule.GetAtomWithIdx(0)
    assert atom.GetFormalCharge() == charge
    assert atom.GetIntProp('source_atom_id') == 12
    assert atom.GetNumImplicitHs() == atom.GetNumExplicitHs() == 0
    assert atom.GetNumRadicalElectrons() == 0


def test_chlorination_uses_shared_conversion_and_retains_aromatic_bonds():
    import molfoundry as mf
    from mfobs.core.gml import graph_gml

    graph = chlorination.molecule(2, 5)
    native = mf.Graph.fromGMLString(graph_gml(graph), add=False)
    for source in (graph, native):
        molecule = to_rdkit(source)
        by_id = {atom.GetIntProp('source_atom_id'): atom for atom in molecule.GetAtoms()}
        assert by_id[0].GetFormalCharge() == -1
        assert all(by_id[u].GetIsAromatic() for u in range(5))
        assert sum(bond.GetBondType() == Chem.BondType.AROMATIC
                   for bond in molecule.GetBonds()) == 5
        assert all(atom.GetNumImplicitHs() == 0 and atom.GetNumRadicalElectrons() == 0
                   for atom in molecule.GetAtoms())


def test_nonchemical_labels_fail_instead_of_changing_the_graph():
    with pytest.raises(ValueError, match='Chemical atom and bond labels'):
        to_rdkit(make_graph({0: 'Position'}))
    with pytest.raises(ValueError, match='Chemical atom and bond labels'):
        to_rdkit(make_graph({0: 'C', 1: 'C'}, [(0, 1, 'src')]))


def test_native_graph_preserves_nonconsecutive_gml_ids_and_bond_labels():
    import molfoundry as mf

    graph = mf.Graph.fromGMLString(
        'graph [ node [ id 9 label "C" ] node [ id 33 label "O" ] '
        'edge [ source 9 target 33 label "=" ] ]', add=False)
    molecule = to_rdkit(graph)
    assert [atom.GetIntProp('source_atom_id') for atom in molecule.GetAtoms()] == [9, 33]
    assert molecule.GetBondWithIdx(0).GetBondType() == Chem.BondType.DOUBLE
    assert all(atom.GetNumImplicitHs() == 0 and atom.GetNumRadicalElectrons() == 0
               for atom in molecule.GetAtoms())


def test_native_graph_without_gml_ids_uses_distinct_native_ids():
    import molfoundry as mf

    graph = mf.Graph.fromSMILES('C')
    molecule = to_rdkit(graph)
    assert [atom.GetIntProp('source_atom_id') for atom in molecule.GetAtoms()] == [
        vertex.id for vertex in graph.vertices]


def test_native_smiles_retains_mixture_multiplicity_and_empty_state():
    import molfoundry as mf
    from mfobs.core.gml import disjoint_union
    from mfobs.core.molecules import smiles
    graph = mf.smiles('C=C', add=False)
    assert smiles(disjoint_union([graph]*3)) == 'C=C.C=C.C=C'
    assert smiles(make_graph({})) == ''
