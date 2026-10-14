"""Check the all-size local-theory contract of a recursive state language.

This is NOT a sampled-chain test, and it does NOT prove rule preservation.
The negative theory consists only of labelled edges and one-centre stars.
For this grammar every complete vertex neighbourhood has one of finitely
many tile contexts. Checking those contexts decides language membership.
"""
from __future__ import annotations

from collections import Counter
import hashlib
import json

from mfobs.core.graph import make_graph, node_labels, node_rows, edge_rows
from mfobs.invariants.chain import ChainLanguage
from mfobs.inputs import InputError


def _graph(data, name):
    """Reject malformed descriptions rather than silently deduplicating nodes."""
    if not isinstance(data, dict) or not {'nodes', 'edges'} <= set(data):
        raise InputError(f'LANGUAGE_CONTRACT: {name} needs nodes and edges')
    if not isinstance(data['nodes'], list) or not isinstance(data['edges'], list):
        raise InputError(f'LANGUAGE_CONTRACT: {name} nodes and edges must be lists')
    seen = set()
    for item in data['nodes']:
        if not isinstance(item, (tuple, list)) or len(item) != 2:
            raise InputError(f'LANGUAGE_CONTRACT: invalid vertex in {name}')
        u, label = item
        if type(u) is not int or u in seen or not isinstance(label, str) or not label:
            raise InputError(f'LANGUAGE_CONTRACT: duplicate/invalid vertex in {name}')
        seen.add(u)
    edges = set()
    for item in data['edges']:
        if not isinstance(item, (tuple, list)) or len(item) != 3:
            raise InputError(f'LANGUAGE_CONTRACT: invalid edge in {name}')
        u, v, label = item
        if (type(u) is not int or type(v) is not int or u not in seen or v not in seen
                or u == v or not isinstance(label, str) or not label):
            raise InputError(f'LANGUAGE_CONTRACT: invalid edge endpoints/label in {name}')
        pair = tuple(sorted((u, v)))
        if pair in edges:
            raise InputError(f'LANGUAGE_CONTRACT: duplicate or parallel edge in {name}')
        edges.add(pair)
    return make_graph(data['nodes'], data['edges'])


def _profiles(graph, context, additions=()):
    """Full stars of target-tile vertices, including distinct external ports.

    Additions are (target vertex, edge label, external-neighbour label). Even
    when both ports use the SAME target vertex, predecessor and successor
    vertices are distinct and contribute separately.
    """
    stars = {u: Counter() for u in node_labels(graph)}
    for u, v, edge in edge_rows(graph):
        stars[u][edge, node_labels(graph)[v]] += 1
        stars[v][edge, node_labels(graph)[u]] += 1
    for u, edge, label in additions:
        stars[u][edge, label] += 1
    for u, label in node_rows(graph):
        yield {'context': context, 'vertex': u, 'label': label,
               'ports': [{'edge': e, 'neighbour': l, 'count': count}
                         for (e, l), count in sorted(stars[u].items())]}


def _check_profile(profile, forbidden):
    label = profile['label']
    at = f"{profile['context']} vertex {profile['vertex']} ({label})"
    if forbidden.strict_signature and label not in forbidden.atom_labels:
        raise InputError(f'LANGUAGE_FORBIDDEN_MISMATCH: {at}: undeclared atom label')
    valence = 0
    for port in profile['ports']:
        edge, neighbour, count = port['edge'], port['neighbour'], port['count']
        pair = tuple(sorted((label, neighbour)))
        if edge not in forbidden.bond_orders:
            raise InputError(f'LANGUAGE_FORBIDDEN_MISMATCH: {at}: undeclared edge {edge}')
        if pair not in forbidden.allowed_bond_pairs[edge]:
            raise InputError(f'LANGUAGE_FORBIDDEN_MISMATCH: {at}: forbidden bond {edge}{pair}')
        valence += count * forbidden.bond_orders[edge]
    maximum = forbidden.max_valence.get(label)
    if maximum is not None and valence > maximum:
        raise InputError(f'LANGUAGE_FORBIDDEN_MISMATCH: {at}: valence {valence} > {maximum}')
    for bound in forbidden.incident_bounds:
        if bound.center != label:
            continue
        count = sum(p['count'] for p in profile['ports']
                    if p['edge'] == bound.edge and p['neighbour'] == bound.neighbour)
        if count > bound.maximum:
            raise InputError(f'LANGUAGE_FORBIDDEN_MISMATCH: {at}: incident bound {bound.name}: {count} > {bound.maximum}')
    return {**profile, 'valence': valence}


def validate_language(language, forbidden=None):
    """Return an all-size membership certificate, or reject before synthesis.

    The language has one C_n, one of two fixed singleton graphs, and M disjoint
    copies of a fixed repeated graph. C_0=zero; C_n=L B^(n-1) R for n>=1.
    The finite context list is induced by this grammar, not a chosen cutoff.
    """
    if not isinstance(language, dict) or not {'chain', 'open', 'closed', 'monomer'} <= set(language):
        raise InputError('LANGUAGE_CONTRACT: chain, open, closed, and monomer required')
    spec = language['chain']
    if not isinstance(spec, dict) or not {'zero', 'left', 'repeat', 'right', 'join_bond'} <= set(spec):
        raise InputError('LANGUAGE_CONTRACT: malformed recursive chain grammar')
    graphs = {name: _graph(spec[name], f'chain.{name}') for name in ['zero', 'left', 'repeat', 'right']}
    finite = {name: _graph(language[name], name) for name in ['open', 'closed', 'monomer']}
    for name in ['left', 'repeat', 'right']:
        tile = spec[name]
        for role in ['incoming', 'outgoing']:
            if role not in tile:
                raise InputError(f'LANGUAGE_CONTRACT: chain.{name} missing {role} port')
            port = tile[role]
            boundary = (name, role) in [('left', 'incoming'), ('right', 'outgoing')]
            if boundary:
                if port is not None:
                    raise InputError(f'LANGUAGE_CONTRACT: chain.{name}.{role} must be null')
            elif type(port) is not int or port not in node_labels(graphs[name]):
                raise InputError(f'LANGUAGE_CONTRACT: chain.{name}.{role} must identify a tile vertex')
    try:
        grammar = ChainLanguage(spec)
    except (KeyError, TypeError, ValueError) as exc:
        raise InputError(f'LANGUAGE_CONTRACT: {exc}') from exc

    profiles = list(_profiles(graphs['zero'], 'chain.zero'))
    for name, graph in finite.items():
        profiles.extend(_profiles(graph, name))
    contexts = []
    # These are precisely the contexts reachable in L B* R. For a repeat
    # both neighbours must be supplied simultaneously (including B-B-B).
    for tile, predecessors, successors in [
        ('left', [None], ['right', 'repeat']),
        ('repeat', ['left', 'repeat'], ['right', 'repeat']),
        ('right', ['left', 'repeat'], [None]),
    ]:
        for previous in predecessors:
            for following in successors:
                name = f'chain.{previous or "boundary"}|{tile}|{following or "boundary"}'
                additions = []
                if previous is not None:
                    additions.append((spec[tile]['incoming'], grammar.bond,
                                      node_labels(graphs[previous])[spec[previous]['outgoing']]))
                if following is not None:
                    additions.append((spec[tile]['outgoing'], grammar.bond,
                                      node_labels(graphs[following])[spec[following]['incoming']]))
                profiles.extend(_profiles(graphs[tile], name, additions))
                contexts.append({'previous': previous, 'tile': tile, 'next': following})
    if forbidden is not None:
        try:
            forbidden.validate()
        except (TypeError, ValueError) as exc:
            raise InputError(f'FORBIDDEN_CONTRACT: {exc}') from exc
        profiles = [_check_profile(profile, forbidden) for profile in profiles]
    digest = lambda value: hashlib.sha256(json.dumps(value, sort_keys=True).encode()).hexdigest()
    return {
        'status': 'VERIFIED_ALL_SIZES' if forbidden is not None else 'NO_PRUNING_THEORY',
        'method': 'exhaustive-complete-local-neighbourhood-contexts',
        'grammar_sha256': grammar.digest,
        'language_sha256': digest(language),
        'forbidden_sha256': digest(forbidden.engine_payload()) if forbidden is not None else None,
        'quantification': 'all n>=0, all M>=0, both finite switch states',
        'chain_contexts': contexts,
        'checked_vertex_profiles': len(profiles),
        'profiles': profiles,
        'proof_document': 'docs/contracts.md',
        'does_not_prove': 'rule preservation of the state language; this remains a separate model obligation',
    }
