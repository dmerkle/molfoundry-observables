"""Shared, application-independent input contracts.

The chemical specification is a permissive local signature/upper-valence
theory.  A domain may narrow it and add structural axioms; neither layer is
advertised as a complete theory of chemistry.  Concrete initial conditions are
lists of graph components with multiplicities, separate from symbolic closure.
"""
from __future__ import annotations

from fractions import Fraction

from mfobs.inputs import InputError


CHEMISTRY_FIELDS = {
    'schema', 'name', 'scope', 'atomTypes', 'bondTypes', 'allowedBondPairs',
}


def _pair(value):
    if not isinstance(value, list) or len(value) != 2 or not all(
            isinstance(item, str) for item in value):
        raise InputError('CHEMISTRY_CONTRACT: bond pairs must be two labels')
    return tuple(sorted(value))


def validate_chemical_specification(data):
    if not isinstance(data, dict) or set(data) != CHEMISTRY_FIELDS:
        raise InputError('CHEMISTRY_CONTRACT: unexpected fields')
    if data['schema'] != 'local-covalent-theory-v1':
        raise InputError('CHEMISTRY_CONTRACT: unsupported schema')
    atoms = data['atomTypes']; bonds = data['bondTypes']
    if not isinstance(atoms, dict) or not atoms:
        raise InputError('CHEMISTRY_CONTRACT: nonempty atomTypes required')
    for label, atom in atoms.items():
        if set(atom) != {'element', 'formalCharge', 'maxValence'}:
            raise InputError(f'CHEMISTRY_CONTRACT: bad atom type {label}')
        if not isinstance(atom['element'], str) or type(atom['formalCharge']) is not int:
            raise InputError(f'CHEMISTRY_CONTRACT: bad atom metadata {label}')
        if Fraction(str(atom['maxValence'])) <= 0:
            raise InputError(f'CHEMISTRY_CONTRACT: bad maxValence {label}')
    if not isinstance(bonds, dict) or not bonds:
        raise InputError('CHEMISTRY_CONTRACT: nonempty bondTypes required')
    for label, bond in bonds.items():
        if set(bond) != {'order'} or Fraction(str(bond['order'])) <= 0:
            raise InputError(f'CHEMISTRY_CONTRACT: bad bond type {label}')
    if set(data['allowedBondPairs']) != set(bonds):
        raise InputError('CHEMISTRY_CONTRACT: allowedBondPairs must cover bondTypes')
    for bond, values in data['allowedBondPairs'].items():
        for value in values:
            pair = _pair(value)
            if not set(pair) <= set(atoms):
                raise InputError(f'CHEMISTRY_CONTRACT: unknown atom in {bond} pair')
    return data


def refinement_certificate(chemistry, forbidden):
    """Check that an engine pruning theory only narrows shared chemistry."""
    validate_chemical_specification(chemistry)
    compiled = forbidden.engine_payload()
    labels = set(compiled['atomLabels'])
    if not labels <= set(chemistry['atomTypes']):
        raise InputError('CHEMISTRY_REFINEMENT: compiled theory has unknown atoms')
    ratios = set()
    for bond, order in compiled['bondOrders'].items():
        if bond not in chemistry['bondTypes']:
            raise InputError(f'CHEMISTRY_REFINEMENT: unknown bond {bond}')
        ratios.add(Fraction(str(order)) /
                   Fraction(str(chemistry['bondTypes'][bond]['order'])))
    if len(ratios) != 1 or next(iter(ratios)) <= 0:
        raise InputError('CHEMISTRY_REFINEMENT: inconsistent bond-order scale')
    scale = next(iter(ratios))
    for label, maximum in compiled['maxValence'].items():
        allowed = Fraction(str(chemistry['atomTypes'][label]['maxValence'])) * scale
        if Fraction(str(maximum)) > allowed:
            raise InputError(f'CHEMISTRY_REFINEMENT: valence of {label} is too permissive')
    for bond, pairs in compiled['allowedBondPairs'].items():
        permitted = {_pair(value) for value in chemistry['allowedBondPairs'][bond]}
        actual = {_pair(list(value)) for value in pairs}
        if not actual <= permitted:
            raise InputError(f'CHEMISTRY_REFINEMENT: {bond} admits a non-base pair')
    return {
        'status': 'VERIFIED_REFINEMENT',
        'shared_chemical_theory': chemistry['name'],
        'compiled_domain_theory': compiled['name'],
        'bond_order_integer_scale': str(scale),
        'scope': ('Signature, bond-pair and upper-valence compatibility only; '
                  'domain structure and rule preservation are separate obligations.'),
    }


def validate_initial_mixture(inputs, data):
    required = {'schema', 'description', 'components'}
    if not isinstance(data, dict) or set(data) != required:
        raise InputError('INITIAL_MIXTURE_CONTRACT: unexpected fields')
    if data['schema'] != 'explicit-component-mixture-v1':
        raise InputError('INITIAL_MIXTURE_CONTRACT: unsupported schema')
    if not isinstance(data['components'], list) or not data['components']:
        raise InputError('INITIAL_MIXTURE_CONTRACT: nonempty components required')
    seen = set()
    for component in data['components']:
        if set(component) != {'id', 'file', 'multiplicity'}:
            raise InputError('INITIAL_MIXTURE_CONTRACT: bad component fields')
        if component['id'] in seen or not isinstance(component['id'], str):
            raise InputError('INITIAL_MIXTURE_CONTRACT: component ids must be unique strings')
        seen.add(component['id'])
        if type(component['multiplicity']) is not int or component['multiplicity'] <= 0:
            raise InputError('INITIAL_MIXTURE_CONTRACT: multiplicities must be positive integers')
        path = inputs.file(component['file'])
        from mfobs.core.engine import mf
        try:
            mf.Graph.fromGMLString(path.read_text(encoding='utf-8'), add=False)
        except Exception as exc:
            raise InputError(f'INITIAL_MIXTURE_CONTRACT: invalid graph {path}') from exc
    return data

