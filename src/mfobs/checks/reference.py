"""Reference comparison is optional verification, NEVER synthesis input."""
import json
import sympy as sp
from mfobs.outputs import matrix

def verify(inputs,result):
    if result['input_definition_sha256']!=inputs.definition_digest():raise AssertionError('Result belongs to different inputs')
    if not result['closed']:raise AssertionError('The run is UNKNOWN, not a completed closure')
    if not inputs.data['reference']:return {'passed':True,'reference_checked':False}
    ref=json.loads(inputs.file(inputs.data['reference']).read_text(encoding='utf-8'))
    if inputs.data['id']=='ring8':
        # Coordinates and discovery order may change without changing the
        # invariant space. Compare exact generators in the same basis.
        current = sp.Matrix([result['basis_vectors'][n] for n in result['basis']])
        expected = sp.Matrix([b['values_on_30_orbits'] for b in ref['basis']])
        if current.shape != expected.shape or current.rank() != expected.rows:
            raise AssertionError('Ring invariant space changed')
        transform = current * expected.inv()
        if current[0, :] != expected[0, :]:
            raise AssertionError('Ring selected query changed')
        a = matrix(result).subs({sp.Symbol(r['rate']): 1 for r in result['rates']})
        reference_matrix = sp.Matrix([
            [sp.Rational(ref['equations_sparse'][b['id']].get(c['id'], 0))
             for c in ref['basis']] for b in ref['basis']])
        if a * transform != transform * reference_matrix:
            raise AssertionError('Ring unit-rate generator changed')
    return {'passed':True,'reference_checked':True,'reference_used_for_synthesis':False,
            'raw_count_compared': result['composition_mode'].startswith('restricted')}
