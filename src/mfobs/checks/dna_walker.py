"""Independent literal multigraph rewrites; no composition or reduction oracle.

These finite checks test the derived identities and incidence encoding. They
are regression evidence, not the proof of closure on every scaffold length.
"""
from collections import Counter
from hashlib import sha256
import json
from fractions import Fraction
from itertools import product
from mfobs.core.graph import matches as incidence_matches
from mfobs.core.rules import rule_graph
from mfobs.models.dna_walker.inputs import identity
from mfobs.models.dna_walker.counting import input_graph
from mfobs.core.engine import mf


def matches(pattern, host):
    """Inject both vertices AND edges, retaining parallel-edge multiplicities."""
    nodes = pattern['nodes']
    candidates = [[n['id'] for n in host['nodes'] if n['label'] == p['label']] for p in nodes]
    for images in product(*candidates):
        if len(set(images)) != len(images):
            continue
        node_map = dict(zip((n['id'] for n in nodes),images))
        edge_candidates = [[e['id'] for e in host['edges']
                            if e['source'] == node_map[p['source']] and e['target'] == node_map[p['target']]]
                           for p in pattern['edges']]
        for edge_images in product(*edge_candidates):
            if len(set(edge_images)) == len(edge_images):
                yield node_map,dict(zip((e['id'] for e in pattern['edges']),edge_images))


def apply_rule(rule, host, match):
    """Literal directed multigraph update using the supplied preserved interface."""
    node_map, edge_map = match
    if len(rule['preserved_nodes']) != len(rule['left']['nodes']) or len(rule['preserved_nodes']) != len(rule['right']['nodes']):
        raise ValueError('Validator only supports edge updates')
    right_nodes = {r:node_map[l] for l,r in rule['preserved_nodes']}
    kept = dict(rule['preserved_edges'])
    removed = {edge_map[e['id']] for e in rule['left']['edges'] if e['id'] not in kept}
    edges = [dict(e) for e in host['edges'] if e['id'] not in removed]
    next_id = max(e['id'] for e in host['edges'])+1
    for edge in rule['right']['edges']:
        if edge['id'] not in kept.values():
            edges.append({'id':next_id,'source':right_nodes[edge['source']], 'target':right_nodes[edge['target']]})
            next_id += 1
    return {'nodes':host['nodes'], 'edges':edges}


def hosts(lengths):
    for n in lengths:
        for site in range(n):
            for extended in (False,True):
                yield {'nodes':[{'id':'b','label':'bimotor'}]+[{'id':i,'label':'chain'} for i in range(n)],
                       'edges':[{'id':i,'source':i,'target':(i+1)%n} for i in range(n)] +
                               [{'id':n,'source':'b','target':site},
                                {'id':n+1,'source':'b','target':(site+int(extended))%n}]}


def validate(result, rules, queries, lengths=range(5,13)):
    lengths = list(lengths)
    if any(n < 5 for n in lengths):
        raise ValueError('Host language requires cycle length at least five')
    # Decode computed patterns once; use plain incidence matching for their counts.
    patterns = {key:input_graph(mf.Rule.fromGMLString(gml,add=False)) for key,gml in result['patterns'].items()}
    reductions = {r['input_gml']:r['reduced_key'] for step in result['steps'] for r in step['reductions']}
    inputs = {gml:input_graph(mf.Rule.fromGMLString(gml,add=False)) for gml in reductions}
    counts = {'hosts':0, 'literal_multigraph_rule_applications':0, 'generator_images_checked':0,
              'original_vs_incidence_observable_counts':0, 'individual_projected_pattern_checks':0}
    for host in hosts(lengths):
        encoded = rule_graph(identity(host,'host'))
        basis = {key:len(incidence_matches(p,encoded)) for key,p in patterns.items()}
        native = {name:sum(1 for _ in matches(query,host)) for name,query in queries.items()}
        for name,query in queries.items():
            assert native[name] == len(incidence_matches(rule_graph(identity(query,name)),encoded)), 'Incidence count mismatch'
            counts['original_vs_incidence_observable_counts'] += 1
        for gml,pattern in inputs.items():
            expected = basis[reductions[gml]] if reductions[gml] is not None else 0
            assert len(incidence_matches(pattern,encoded)) == expected, 'Invalid structural reduction'
            counts['individual_projected_pattern_checks'] += 1
        for rule in rules:
            increments = {name:0 for name in queries}
            for match in matches(rule['left'],host):
                after = apply_rule(rule,host,match)
                for name,query in queries.items():
                    increments[name] += sum(1 for _ in matches(query,after))-native[name]
                counts['literal_multigraph_rule_applications'] += 1
            for name,actual in increments.items():
                image = result['raw_images'][rule['name']+'/'+name]
                predicted = sum(Fraction(c)*native[n] for n,c in image.items())
                assert actual == predicted, (rule['name'],name,actual,predicted)
                counts['generator_images_checked'] += 1
        counts['hosts'] += 1
    return {'status':'PASS','cycle_lengths':lengths,**counts}


def verified_exports(directory, record, provenance):
    expected={**record['fragger_exports_sha256'],**provenance['executed_harness_sha256']}
    for name,digest in expected.items():
        if sha256((directory/name).read_bytes()).hexdigest()!=digest:
            raise ValueError('Recovered evidence hash mismatch: '+name)
    return json.loads((directory/'walker-export.json').read_text(encoding='utf-8')),json.loads((directory/'walker-incidence.json').read_text(encoding='utf-8'))


def compare_incidence(record,incidence):
    """Compare eight retained FRAGGER contributions with the supplied Python result."""
    rows=[]
    for channel in incidence['individual_channels']:
        name=channel['rule']+'/'+channel['coordinate'];counts=Counter()
        for equation in channel['output']['equations']:
            if equation['coordinate']!=channel['coordinate']:raise ValueError('Unexpected left-hand side')
            for term in equation['terms']:
                if term['rates']!=[channel['rule']]:raise ValueError('Unexpected rate monomial')
                counts[term['coordinate']]+=Fraction(str(term['coefficient']))
        computed={k:v for k,v in counts.items() if v}
        expected={k:Fraction(v) for k,v in record['raw_images'][name].items() if v}
        rows.append({'Rule / query':name,'Recovered FRAGGER image':str(computed),'Computed MolFoundry image':str(expected),'Equal':computed==expected})
    if len(rows)!=8 or not all(row['Equal'] for row in rows):raise AssertionError('Recovered contribution mismatch')
    return rows
