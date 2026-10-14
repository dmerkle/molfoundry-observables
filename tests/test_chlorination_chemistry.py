"""Independent full-molecule checks of the synthetic charged Cp variant."""
from mfobs.chemical_depiction import to_rdkit
from mfobs.core.graph import edge_rows, node_labels
import pytest
from rdkit import Chem

from mfobs.core.conditions import matches
from mfobs.models.chlorination import inputs as chemistry


def molecular_key(graph):
    return Chem.MolToSmiles(Chem.RemoveHs(to_rdkit(graph)))


def test_all_masks_lengths_have_explicit_h_and_valid_charged_aromaticity():
    for n in (0, 1, 2, 3):
        for mask in range(16):
            graph = chemistry.molecule(n, mask)
            occupied = mask.bit_count()
            expected = {'C': 6+n, 'H': 7+2*n-occupied}
            if occupied:
                expected['Cl'] = occupied
            assert chemistry.formula(graph) == expected
            assert chemistry.formal_charge(graph) == -1
            assert graph.numVertices == 13+3*n
            assert node_labels(graph)[0] == 'C-'
            assert {label for _, _, label in edge_rows(graph)} == {'-', ':'}
            assert sum(label == ':' for _, _, label in edge_rows(graph)) == 5
            assert not any(node_labels(graph)[v if u == 0 else u] == 'H'
                           for u, v, _ in edge_rows(graph) if 0 in (u, v))
            mol = to_rdkit(graph)
            assert mol.GetNumAtoms() == graph.numVertices
            assert sum(atom.GetFormalCharge() for atom in mol.GetAtoms()) == -1
            assert all(atom.GetNumRadicalElectrons() == 0 and atom.GetNumImplicitHs() == 0
                       for atom in mol.GetAtoms())
            for atom in mol.GetAtoms():
                expected_valence = 3 if atom.GetFormalCharge() == -1 else {
                    'C': 4, 'H': 1, 'Cl': 1}[atom.GetSymbol()]
                assert atom.GetTotalValence() == expected_valence
            # Aromatic graph weights are bookkeeping, not the C- chemical valence.
            assert chemistry.bond_order_sums(graph)[0] == 4
            assert len(matches(chemistry.seed_graph(), graph)) == occupied
    assert molecular_key(chemistry.molecule()) == 'C[c-]1cccc1'
    assert chemistry.formula(chemistry.methane()) == {'C': 1, 'H': 4}
    assert chemistry.formal_charge(chemistry.hcl()) == 0


def test_reflection_yields_ten_distinct_head_molecules():
    for mask in range(16):
        reverse = int(f'{mask:04b}'[::-1], 2)
        assert molecular_key(chemistry.molecule(2, mask)) == molecular_key(chemistry.molecule(2, reverse))
    assert len({molecular_key(chemistry.molecule(0, mask)) for mask in range(16)}) == 10


@pytest.mark.parametrize('n,mask', [(-1, 0), (0, 16), (0, -1), (True, 0), (0, False)])
def test_invalid_constructor_inputs(n, mask):
    with pytest.raises(ValueError):
        chemistry.molecule(n, mask)
