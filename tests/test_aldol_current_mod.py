from copy import deepcopy

import pytest

from examples.tcs.compare_conditions import compare_conditions
from examples.tcs.current_mod import compare_recorded_mod


def test_current_mod_plain_and_translated_noedge_match_fresh_composition(tmp_path):
    fresh = compare_conditions(tmp_path)
    result = compare_recorded_mod(fresh)
    plain = result['runs']['plain']
    assert plain['totals'] == {'after': 13, 'before': 18}
    assert len(plain['rows']) == 22
    assert next(r for r in plain['rows'] if r['id'] == 978)['plain_class'] == 'new-f7c8786fa6b1'
    conditioned = result['runs']['reaction_and_query_noedge']
    assert conditioned['totals'] == {'after': 10, 'before': 18}
    assert len(conditioned['rows']) == 23
    assert sum(r['after'] != r['before'] for r in conditioned['rows']) == 18

    # Equal totals alone must not pass: move one contribution to another class.
    wrong = deepcopy(fresh)
    rows = wrong['variants']['reaction_and_query_noedge']['rows']
    donor = next(r for r in rows if r['after'])
    receiver = next(r for r in rows if r['key_sha256'] != donor['key_sha256'])
    donor['after'] -= 1
    receiver['after'] += 1
    with pytest.raises(ValueError, match='multiplicities differ'):
        compare_recorded_mod(wrong)
