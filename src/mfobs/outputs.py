"""Computed analysis records and ODEs from rate-separated coefficients."""
from __future__ import annotations
from pathlib import Path
from fractions import Fraction
import csv,json,re,shutil
from xml.etree import ElementTree as ET
import sympy as sp
from mfobs.core.exact_coefficients import coefficient


def support_record(support, sources=(), *, events=None):
    """Export one support result, with native input rules and optional events."""
    result = {"schema_version": 1, **support.record()}
    if sources:
        result["source_rules"] = {r.identifier: {"rate": r.rate, "gml": r.rule.getGMLString()}
                                  for r in sources}
    if events is not None:
        result["events"] = events
    return result


def matrix(result,normalized=False):
    names=result['basis'];index={n:i for i,n in enumerate(names)};mat=sp.zeros(len(names))
    if 'equations_by_channel' in result:
        equations = result['equations_by_channel']
        rates = {name: sp.Symbol(name) for rows in equations.values() for name in rows}
    else:
        equations = result['equations_by_reaction']
        rates = {r['id']: sp.Symbol(r['rate']) for r in result['rates']}
    for source, rows in equations.items():
        for channel, values in rows.items():
            for target, c in values.items():
                mat[index[source], index[target]] += rates[channel] * coefficient(c, result.get('coefficient_parameters', ()))
    if normalized:
        norm=result.get('normalisation',{})
        scales=[sp.Rational(norm.get(n,{}).get('scale','1')) for n in names]
        mat=sp.diag(*(1/s for s in scales))*mat*sp.diag(*scales)
    return mat


def ode_text(result,normalized=False):
    if not result['closed']:
        return 'UNKNOWN: no closed ODE is certified by this run. Inspect result.json and trace.json.\n'
    names=result['basis'];a=matrix(result,normalized);norm=result.get('normalisation',{})
    symbols=sp.Matrix([sp.Symbol(('z_' if normalized else 'mu_')+n) for n in names])
    combined=result.get('coordinate_kind')=='linear-combination-of-source-observables'
    meaning = ('the exported exact linear observable expression' if result.get('observable_expressions')
               else 'the indicated rational linear combination of source observables' if combined
               else 'the RAW injective pattern count O')
    lines=[f"Model: {result['model']}",f"Status: {result['status']}",
           f"Basis order: {', '.join(names)}",'',f'mu_O(t) = expectation of {meaning}.']
    if result.get('observable_expressions'):
        lines += ['Coordinates are exact linear combinations; see observable_expressions in result.json.']
    if combined:
        lines += [f"{name} = "+' + '.join(f'({c})*{source}' for source,c in combo.items())
                  for name,combo in result['basis_combinations'].items()]
    if normalized:
        lines+=['z_O = mu_O / scale_O; this is only a change of coordinates.']
        lines += [f"z_{n} = mu_{n}/{norm.get(n,{}).get('scale','1')} = E[{norm.get(n,{}).get('meaning',n)}]" for n in names]
    lines+=['']+[f'd{symbols[i]}/dt = {sp.expand((a*symbols)[i])}' for i in range(len(names))]
    lines+=['',result.get('proof_status','')]
    return '\n'.join(lines)+'\n'


def rounds(result):
    rows=[];basis=result['basis'];depths=result['depths'];total=sum(d==0 for d in depths.values())
    rows.append({'round':0,'expanded':[],'before':0,'new':[k for k in basis if depths[k]==0],'after':total})
    completed={s['expanded'] for s in result['steps']}
    for depth in range(max(depths.values(),default=0)+1):
        expanded=[k for k in basis if depths[k]==depth and k in completed]
        if not expanded:break
        new=[k for k in basis if depths[k]==depth+1]
        rows.append({'round':depth+1,'expanded':expanded,'before':total,'new':new,'after':total+len(new)})
        total+=len(new)
    return rows


def ode_json(result):
    """Portable machine-readable form of the computed linear ODE."""
    if not result['closed']:
        return {'schema': 'mfobs-linear-ode-v1', 'status': 'UNKNOWN'}
    names=result['basis'];a=matrix(result)
    y=sp.Matrix([sp.Symbol('mu_'+name) for name in names]);rhs=a*y
    return {
        'schema': 'mfobs-linear-ode-v1', 'model': result['model'],
        'independent_variable': 't',
        'states': [{'id':'mu_'+name,'observable':name,
                    'initial_value':result.get('initial_values',{}).get(name)}
                   for name in names],
        'parameters': sorted({r['rate'] for r in result['rates']} |
                             set(result.get('coefficient_parameters',()))),
        'rhs': {name:str(sp.expand(rhs[index])) for index,name in enumerate(names)},
    }


def _sbml_id(value):
    value=re.sub(r'[^A-Za-z0-9_]', '_', value)
    return value if value and (value[0].isalpha() or value[0]=='_') else '_'+value


def _mathml(expression):
    """Small exact Content-MathML printer for the exported linear systems."""
    ns='http://www.w3.org/1998/Math/MathML'
    tag=lambda name:f'{{{ns}}}{name}'
    def emit(value):
        value=sp.sympify(value)
        if isinstance(value,sp.Symbol):
            node=ET.Element(tag('ci'));node.text=str(value);return node
        if value.is_Integer:
            node=ET.Element(tag('cn'),type='integer');node.text=str(value);return node
        if value.is_Rational:
            node=ET.Element(tag('apply'));ET.SubElement(node,tag('divide'))
            top=ET.SubElement(node,tag('cn'),type='integer');top.text=str(value.p)
            bottom=ET.SubElement(node,tag('cn'),type='integer');bottom.text=str(value.q)
            return node
        operations=((sp.Add,'plus'),(sp.Mul,'times'),(sp.Pow,'power'))
        for kind,name in operations:
            if isinstance(value,kind):
                node=ET.Element(tag('apply'));ET.SubElement(node,tag(name))
                for argument in value.args:node.append(emit(argument))
                return node
        raise ValueError(f'Unsupported SBML MathML expression: {value}')
    root=ET.Element(tag('math'));root.append(emit(expression));return root


def sbml(result):
    """SBML Level 3 Core export using parameters and exact RateRules.

    Observable expectations are mathematical state variables rather than
    chemical species, so SBML parameters with ``constant=false`` are the
    faithful Core representation.
    """
    if not result['closed']:
        return None
    ns='http://www.sbml.org/sbml/level3/version2/core'
    root=ET.Element('sbml',xmlns=ns,level='3',version='2')
    model=ET.SubElement(root,'model',id=_sbml_id(result['model']))
    params=ET.SubElement(model,'listOfParameters')
    names=result['basis'];state_ids={name:_sbml_id('mu_'+name) for name in names}
    initial=result.get('initial_values',{})
    for name in names:
        attrs={'id':state_ids[name],'name':f'E[{name}]','constant':'false'}
        if name in initial:
            value=sp.sympify(initial[name])
            if not value.free_symbols:attrs['value']=str(float(value))
        ET.SubElement(params,'parameter',attrs)
    parameter_names=sorted({r['rate'] for r in result['rates']} |
                           set(result.get('coefficient_parameters',())))
    initial_state=result.get('input_initial_state',{})
    for name in parameter_names:
        attrs={'id':_sbml_id(name),'name':name,'constant':'true'}
        if name in initial_state and isinstance(initial_state[name],(int,float)):
            attrs['value']=str(initial_state[name])
        ET.SubElement(params,'parameter',attrs)
    rules=ET.SubElement(model,'listOfRules')
    a=matrix(result);symbols=sp.Matrix([sp.Symbol(state_ids[n]) for n in names])
    substitutions={sp.Symbol(n):sp.Symbol(_sbml_id(n)) for n in parameter_names}
    for index,name in enumerate(names):
        rate=ET.SubElement(rules,'rateRule',variable=state_ids[name])
        expression=sp.expand((a*symbols)[index]).xreplace(substitutions)
        rate.append(_mathml(expression))
    ET.indent(root,space='  ')
    return '<?xml version="1.0" encoding="UTF-8"?>\n'+ET.tostring(root,encoding='unicode')+'\n'


def export(result,directory):
    linear_system = 'closed' in result and 'basis' in result
    serialized = None if linear_system else json.dumps(result, indent=2) + '\n'
    out=Path(directory);out.mkdir(parents=True,exist_ok=True)
    # Never leave a stale closed ODE after overwriting it with UNKNOWN.
    for name in ['ode.txt','ode-normalized.txt','ode.tex','ode.json','model.sbml.xml','matrix.csv','coefficients-by-rate.csv','verification.json','trajectory.csv',
                 'trace.json','iterations.csv','observable-index.json',
                 'basis-combinations.json','transformation.csv','seed-reconstruction.json']:
        (out/name).unlink(missing_ok=True)
    patterns=out/'identity-patterns'
    if patterns.is_symlink():raise ValueError('Refuse to overwrite a symlinked pattern output')
    if patterns.exists():shutil.rmtree(patterns)
    observables=out/'observables'
    if observables.is_symlink():raise ValueError('Refuse to overwrite a symlinked observable output')
    if observables.exists():shutil.rmtree(observables)
    # Support/event analyses and comparisons retain their own result schema.
    if not linear_system:
        (out/'result.json').write_text(serialized, encoding='utf-8')
        return out
    result['rounds']=rounds(result)
    (out/'result.json').write_text(json.dumps(result,indent=2)+'\n', encoding='utf-8')
    (out/'ode.txt').write_text(ode_text(result), encoding='utf-8')
    (out/'ode-normalized.txt').write_text(ode_text(result,True), encoding='utf-8')
    trace={k:result[k] for k in ['model','closed','basis','depths','steps','rounds','raw_compositions','ordered_composition_calls','expansions','unknown']}
    (out/'trace.json').write_text(json.dumps(trace,indent=2)+'\n', encoding='utf-8')
    with (out/'iterations.csv').open('w',newline='', encoding='utf-8') as f:
        writer=csv.writer(f);writer.writerow(['round','expanded','before','new','after'])
        for r in result['rounds']:writer.writerow([r['round'],','.join(r['expanded']),r['before'],','.join(r['new']),r['after']])
    if result['closed']:
        (out/'ode.json').write_text(json.dumps(ode_json(result),indent=2)+'\n', encoding='utf-8')
        (out/'model.sbml.xml').write_text(sbml(result), encoding='utf-8')
        a=matrix(result);names=result['basis']
        with (out/'matrix.csv').open('w',newline='', encoding='utf-8') as f:
            writer=csv.writer(f);writer.writerow(['row/column',*names])
            for i,n in enumerate(names):writer.writerow([n,*[str(v) for v in a.row(i)]])
        with (out/'coefficients-by-rate.csv').open('w',newline='', encoding='utf-8') as f:
            by_channel = 'equations_by_channel' in result
            writer=csv.writer(f);writer.writerow(['source','rate_channel' if by_channel else 'reaction','target','coefficient_after_per_match_factor'])
            equations = result['equations_by_channel'] if by_channel else result['equations_by_reaction']
            for source,rows in equations.items():
                for channel,vals in rows.items():
                    for target,c in vals.items():writer.writerow([source,channel,target,c])
        ys=sp.Matrix([sp.Symbol('mu_'+n) for n in names]);rhs=a*ys
        tex=['% Include with amsmath. Generated from the computed coefficients.','\\begin{align*}']
        tex += [f'\\frac{{d}}{{dt}} {sp.latex(ys[i])} &= {sp.latex(sp.expand(rhs[i]))} '+r'\\' for i in range(len(names))]
        tex+=['\\end{align*}'];(out/'ode.tex').write_text('\n'.join(tex)+'\n', encoding='utf-8')
    if result.get('coordinate_kind')=='linear-combination-of-source-observables':
        (out/'basis-combinations.json').write_text(json.dumps(result['basis_combinations'],indent=2)+'\n', encoding='utf-8')
        (out/'seed-reconstruction.json').write_text(json.dumps(result['seed_reconstruction'],indent=2)+'\n', encoding='utf-8')
        with (out/'transformation.csv').open('w',newline='', encoding='utf-8') as f:
            writer=csv.writer(f);writer.writerow(['compressed/source',*result['source_basis']])
            for name,row in zip(result['basis'],result['transformation_matrix']):writer.writerow([name,*row])
        return out
    if 'observable_expressions' in result:
        # No fictitious single GML file for a linear combination.
        observables.mkdir(parents=True, exist_ok=True)
        for name, terms in result['observable_expressions'].items():
            index = []
            for i, term in enumerate(terms):
                filename = f'{name}-{i}.gml'
                (observables / filename).write_text(term['gml'], encoding='utf-8')
                index.append({'coefficient': term['coefficient'], 'file': filename})
            (observables / (name + '.json')).write_text(json.dumps(index, indent=2) + '\n', encoding='utf-8')
        return out
    observables.mkdir(parents=True,exist_ok=True)
    observable_index=[]
    for name in result['basis']:
        filename=name+'.gml'
        (observables/filename).write_text(result['identity_patterns'][name]['gml'], encoding='utf-8')
        observable_index.append({'id':name,'file':'observables/'+filename,
                                 'depth':result['depths'][name],
                                 'origin':result.get('basis_origins',{}).get(name),
                                 'meaning':result['identity_patterns'][name]['meaning'],
                                 'semantics':'identity rule id(O): O <- O -> O; its left object O is the counted observable'})
    (out/'observable-index.json').write_text(json.dumps(observable_index,indent=2)+'\n', encoding='utf-8')
    for round_ in result['rounds']:
        d=out/'identity-patterns'/f"round-{round_['round']}";d.mkdir(parents=True,exist_ok=True)
        for name in result['basis']:
            if result['depths'][name]<=round_['round']:
                (d/(name+'.gml')).write_text(result['identity_patterns'][name]['gml'], encoding='utf-8')
    return out


def export_summary(result, directory, *, name=None, details=False):
    """Save the readable calculation; retain full records only when requested."""
    from mfobs.formatting import observable_rows, render

    out = Path(directory)
    out.mkdir(parents=True, exist_ok=True)
    if details:
        export(result, out)
    (out / 'summary.txt').write_text(render(result, name=name) + '\n', encoding='utf-8')
    for group in ('variants', 'queries'):
        if group in result:
            for index, (label, child) in enumerate(result[group].items()):
                folder = label if re.fullmatch(r'[A-Za-z0-9_-]{1,80}', label) else f'calculation-{index}'
                export_summary(child, out / folder, name=label)
            return out
    rows = observable_rows(result)
    if rows:
        lines = ['Observable definitions', 'count(file) = raw injective count of the conditional pattern in file.', '']
        patterns = out / 'observables'
        for index, row in enumerate(rows):
            label = row['name']
            stem = label if re.fullmatch(r'[A-Za-z0-9_-]{1,80}', label) else f'O{index}'
            summands = []
            for i, term in enumerate(row['terms']):
                patterns.mkdir(exist_ok=True)
                filename = f'{stem}-{i}.gml'
                (patterns / filename).write_text(term['gml'], encoding='utf-8')
                summands.append(f"({term['coefficient']}) count(observables/{filename})")
            expression = ' + '.join(summands) or row['text']
            if expression == label:
                lines.append('Expression: ' + expression)
            elif expression:
                lines.append(f'{label} = {expression}')
            else:
                lines.append(f'{label}: no literal definition retained in this calculation')
            if row['meaning']:
                lines.append('  ' + row['meaning'])
        (out / 'observables.txt').write_text('\n'.join(lines) + '\n', encoding='utf-8')
    linear = result.get('linear') or result
    if linear.get('closed') or linear.get('status') == 'CLOSED':
        text = (ode_text(linear) if 'closed' in linear else render(result, name=name, section='ode') + '\n')
        (out / 'ode.txt').write_text(text, encoding='utf-8')
    else:
        (out / 'ode.txt').unlink(missing_ok=True)
        (out / 'ode-normalized.txt').unlink(missing_ok=True)
    if 'measurements' in result and result['measurements']:
        with (out / 'measurements.csv').open('w', newline='', encoding='utf-8') as stream:
            writer = csv.DictWriter(stream, fieldnames=list(result['measurements'][0]))
            writer.writeheader()
            writer.writerows(result['measurements'])
    if 'compositions' in result:
        patterns = out / 'compositions'
        patterns.mkdir(exist_ok=True)
        lines = []
        for index, (label, terms) in enumerate(result['compositions'].items()):
            stem = label if re.fullmatch(r'[A-Za-z0-9_-]{1,80}', label) else f'C{index}'
            summands = []
            for i, term in enumerate(terms):
                filename = f'{stem}-{i}.gml'
                (patterns / filename).write_text(term['gml'], encoding='utf-8')
                summands.append(f"({term['coefficient']}) compositions/{filename}")
            lines.append(f"{label} = " + (' + '.join(summands) or '0'))
        (out / 'compositions.txt').write_text('\n'.join(lines) + '\n', encoding='utf-8')
    return out


def dna_walker_readout(result, rules):
    """Name discovered counts AFTER closure; normalize by graph automorphisms.

    This presentation layer is not called by the image oracle or worklist.
    The half-step displacement is the change of mean leg position along each
    rule's preserved directed scaffold edge. Rates are per physical move.
    """
    from mfobs.core.engine import mf
    from mfobs.models.dna_walker.inputs import convert_rule
    from mfobs.models.dna_walker.counting import input_graph
    from mfobs.core.graph import matches, node_rows
    from mfobs.core.rules import rule_graph
    meanings, symmetries, scales = {}, {}, {}
    for name, expression in result['basis'].items():
        if len(expression) != 1:
            raise ValueError('Configuration readout needs single-pattern basis elements')
        key, scale = next(iter(expression.items()))
        graph = input_graph(mf.Rule.fromGMLString(result['patterns'][key], add=False))
        kind = 'Craw' if sum(l == 'chain' for _,l in node_rows(graph)) == 1 else 'E'
        meanings[name], scales[name] = kind, Fraction(scale)
        symmetries[kind] = len(matches(graph,graph))
    if set(meanings.values()) != {'E','Craw'}:
        raise ValueError('Both configurations must have been discovered')
    raw = {s['channel']+'/'+meanings[s['source']]:
           {meanings[n]: str(Fraction(c)*scales[n]/scales[s['source']])
            for n,c in s['coordinates'].items()} for s in result['steps']}
    movement, multipliers = {}, {}
    for definition in rules:
        rule = convert_rule(definition)
        multipliers[rule.name] = Fraction(1,len(matches(rule_graph(rule),rule_graph(rule))))
        labels = {n['id']:n['label'] for n in definition['left']['nodes']}
        scaffold, = [e for e in definition['left']['edges'] if labels[e['source']] == 'chain']
        positions = {scaffold['source']:0, scaffold['target']:1}
        right_to_left = {r:l for l,r in definition['preserved_nodes']}
        def mean_leg_position(side):
            graph = definition[side]
            node_labels = {n['id']:n['label'] for n in graph['nodes']}
            targets = [e['target'] for e in graph['edges'] if node_labels[e['source']] == 'bimotor']
            if side == 'right':
                targets = [right_to_left[t] for t in targets]
            return Fraction(sum(positions[t] for t in targets),len(targets)),len(set(targets))
        before, distinct = mean_leg_position('left')
        after, _ = mean_leg_position('right')
        movement[rule.name] = {'displacement':str(after-before),
                               'normalized_activity':'C' if distinct == 1 else 'E'}
    physical, equations = {}, {'E':sp.Integer(0), 'C':sp.Integer(0)}
    rename = lambda name:'C' if name == 'Craw' else name
    for key, vector in raw.items():
        channel, source = key.split('/')
        coefficients = {rename(n):sp.Rational(c)*multipliers[channel]*symmetries[n]/symmetries[source]
                        for n,c in vector.items()}
        physical[channel+'/'+rename(source)] = {n:str(c) for n,c in coefficients.items()}
        equations[rename(source)] += sp.Symbol(channel)*sum(c*sp.Symbol(n) for n,c in coefficients.items())
    E,C = sp.symbols('E C')
    stationary, = sp.solve([equations['E'], equations['C'], E+C-1], (E,C), dict=True)
    velocity = sum(sp.Rational(r['displacement'])*sp.Symbol(k)*sp.Symbol(r['normalized_activity'])
                   for k,r in movement.items())
    return {'raw_images':raw, 'display_names':meanings,
            'comparison':{'physical_images':physical,
                'physical_equations':{n:str(sp.expand(e)) for n,e in equations.items()},
                'stationary_probabilities':{str(n):str(sp.factor(v)) for n,v in stationary.items()},
                'stationary_velocity':str(sp.factor(velocity.subs(stationary))),
                'movement_readout_from_rules':movement,
                'per_match_multipliers':{n:str(v) for n,v in multipliers.items()}}}
