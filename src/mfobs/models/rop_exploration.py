"""Bounded literal graph expansion of a discovered terminal-path model.

Direct MolFoundry matching and DPO application check each displayed transition
against the exported event model. This expansion is not used for synthesis.
"""
from mfobs.core.graph import edge_rows, node_labels
from collections import Counter, deque

from mfobs.core.molfoundry_execution import count_matches, matches
from mfobs.core.engine import mf
from mfobs.core.gml import graph_gml
from mfobs.core.graph import apply_match
from mfobs.core.recursive_words import TerminalGrowthGrammar, WordState
from mfobs.models.aba.projection import TerminalWordGraphModel
from mfobs.models.aba.analysis import word_coordinates
from mfobs.numerics.events import event_system


def _operand(gml):
    return mf.Rule.fromGMLString(gml, add=False)


def _state(graph, alphabet, gated):
    """Read a literal successor's path and free tokens, without a rewrite formula."""
    labels = node_labels(graph)
    neighbours = {v: [] for v in labels}
    for u, v, bond in edge_rows(graph):
        neighbours[u].append((v, bond)); neighbours[v].append((u, bond))
    root = next(v for v in labels if labels[v] == 'Root')
    current = next(v for v, bond in neighbours[root] if bond == 'kind')
    letters, attached = [], set()
    while True:
        attached.add(current)
        link = next(v for v, bond in neighbours[current] if bond == 'src')
        target = next(v for v, bond in neighbours[link] if bond == 'tgt')
        if labels[link] == 'Terminal':
            assert labels[target] == 'End'
            break
        assert target not in attached
        letters.append(next(labels[v][7:] for v, bond in neighbours[target] if bond == 'kind'))
        current = target
    pools = Counter(next(labels[v][7:] for v, bond in neighbours[p] if bond == 'kind')
                    for p in labels if labels[p] == 'Position' and p not in attached)
    opened = not any(bond == 'bound' for _, _, bond in edge_rows(graph)) if gated else True
    return WordState(tuple(letters), tuple(pools[a] for a in alphabet), opened)


def expand(record, *, word='ABAB', pools=None, opened=True, depth=2):
    if not isinstance(depth, int) or not 0 <= depth <= 3:
        raise ValueError('Expansion depth must be between 0 and 3')
    pools = {'A': 2, 'B': 1} if pools is None else pools
    word_coordinates(record, word, pools, opened)  # check the supplied state
    support = record['support_discovery']
    alphabet = tuple(pools)
    gated = any(p['kind'] == 'gate' for p in support['patterns'].values())
    # No motif or suffix classes enter host construction.
    host_model = TerminalWordGraphModel(TerminalGrowthGrammar('host', alphabet, (), gated, ''))
    initial = WordState(tuple(word), tuple(pools.values()), opened)
    hosts = [host_model.state_graph(initial)]
    states = [initial]; indices = {initial: 0}; depths = [0]; queue = deque([0])
    rules = {name: _operand(r['gml']) for name, r in support['source_rules'].items()}
    patterns = {name: _operand(p['gml']) for name, p in support['patterns'].items()}
    events = event_system(record)
    terms = record['pjc']['operator']['terms']
    transitions, counts = [], []
    def coordinates(host, state):
        values = {n: count_matches(p, host) for n, p in patterns.items()}
        assert values == word_coordinates(record, ''.join(state.word), dict(zip(alphabet, state.pools)), state.gate_open)
        return values
    counts.append(coordinates(hosts[0], initial))
    while queue:
        i = queue.popleft()
        if depths[i] == depth:
            continue
        for name, rule in rules.items():
            observed, targets, maps = Counter(), Counter(), {}
            for mapping in matches(rule, hosts[i]):
                target = apply_match(rule, hosts[i], mapping)
                if target is None:
                    continue
                state = _state(target, alphabet, gated)
                if state not in indices:
                    indices[state] = len(states); states.append(state); hosts.append(target)
                    depths.append(depths[i]+1); queue.append(indices[state])
                    counts.append(coordinates(target, state))
                j = indices[state]; targets[j] += 1; maps[j] = mapping
                delta = tuple(counts[j][n]-counts[i][n] for n in events.names)
                if any(delta): observed[delta] += 1
            expected = Counter()
            for k, term in enumerate(terms):
                if term['rule'] == name:
                    weight = events.activities[k].subs(dict(zip(events.symbols, (counts[i][n] for n in events.names))))
                    if weight: expected[tuple(map(int, events.deltas[k]))] += weight
            assert observed == expected, (i, name, observed, expected)
            transitions.extend({'source': i, 'target': j, 'rule': name, 'multiplicity': count,
                'example_match': maps[j], 'rate': support['source_rules'][name]['rate']}
                for j, count in targets.items())
    return {'depth': depth, 'initial': 0, 'status': 'PASS',
        'method': 'literal conditional graph matching and DPO application; event distribution checked afterwards',
        'states': [{'id': i, 'depth': depths[i], 'word': ''.join(s.word),
                    'pools': dict(zip(alphabet, s.pools)), 'gate_open': s.gate_open if gated else None,
                    'coordinates': counts[i], 'gml': graph_gml(hosts[i])} for i, s in enumerate(states)],
        'transitions': transitions}
