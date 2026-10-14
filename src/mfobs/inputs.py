"""Declarative model input, explicit file hashes, strict validation.

Importing an adapter is executing trusted Python code, not sandboxed input.
Numerical rates and initial states are separate from the symbolic synthesis.
"""
from __future__ import annotations
from dataclasses import dataclass
from copy import deepcopy
from fractions import Fraction
from pathlib import Path
import hashlib, importlib, json, re

class InputError(ValueError): pass

@dataclass
class ModelManifest:
    path: Path
    data: dict
    @property
    def directory(self): return self.path.parent
    def file(self, value):
        p=(self.directory/value).resolve()
        if not p.is_file(): raise InputError(f'Input file does not exist: {p}')
        return p
    def read(self, key): return json.loads(self.file(self.data[key]).read_text(encoding='utf-8'))
    def definition_digest(self):
        # Numerical rates and initial state change initial trajectories, not the
        # generator's symbolic closure. Reference output NEVER feeds synthesis.
        fields={k:self.data[k] for k in ('adapter','rules','observables','initial_observables','state_language','reduction','forbidden','normalisation')}
        fields['overlap_policy']={k:self.data['composition'][k] for k in ['connected','include_empty']}
        files=[v['file'] for v in self.data['rules']+self.data['observables']]
        files += [self.data['state_language'],self.data['reduction']]
        if self.data['forbidden']:files.append(self.data['forbidden'])
        if self.data.get('specifications'):
            fields['specifications']=self.data['specifications']
            files += list(self.data['specifications'].values())
        hashes={v:hashlib.sha256(self.file(v).read_bytes()).hexdigest() for v in files}
        return hashlib.sha256(json.dumps({'fields':fields,'files':hashes},sort_keys=True).encode()).hexdigest()
    def backend_digest(self):
        """Hash rules and language proof inputs, independently of the query.

        A language adapter may be reviewed for a fixed grammar while accepting
        different identity observables.  Such an adapter binds to this digest;
        complete analysis records still use ``definition_digest()``, which also
        binds the selected observable files and normalization.
        """
        fields={k:self.data[k] for k in (
            'adapter','rules','state_language','reduction','forbidden'
        )}
        fields['overlap_policy']={
            k:self.data['composition'][k] for k in ['connected','include_empty']
        }
        files=[v['file'] for v in self.data['rules']]
        files += [self.data['state_language'],self.data['reduction']]
        if self.data['forbidden']:files.append(self.data['forbidden'])
        if self.data.get('specifications'):
            fields['specifications']=self.data['specifications']
            files += list(self.data['specifications'].values())
        hashes={v:hashlib.sha256(self.file(v).read_bytes()).hexdigest() for v in files}
        return hashlib.sha256(json.dumps(
            {'fields':fields,'files':hashes},sort_keys=True
        ).encode()).hexdigest()
    def specification_certificate(self):
        refs=self.data.get('specifications')
        if not refs:return None
        from mfobs.specifications import validate_chemical_specification,refinement_certificate
        chemistry=json.loads(self.file(refs['chemistry']).read_text(encoding='utf-8'))
        validate_chemical_specification(chemistry)
        # The domain is deliberately adapter-defined, but it must be an actual
        # JSON object and must be the same file used as the state language.
        domain=json.loads(self.file(refs['domain']).read_text(encoding='utf-8'))
        if not isinstance(domain,dict):raise InputError('DOMAIN_CONTRACT: JSON object required')
        if self.file(refs['domain'])!=self.file(self.data['state_language']):
            raise InputError('DOMAIN_CONTRACT: specifications.domain must equal state_language')
        if not self.data['forbidden']:
            raise InputError('CHEMISTRY_CONTRACT: shared chemistry requires compiled forbidden theory')
        return refinement_certificate(chemistry,self.forbidden_spec())
    def initial_mixture(self):
        reference=self.data.get('initial_mixture')
        if not reference:return None
        from mfobs.specifications import validate_initial_mixture
        return validate_initial_mixture(self,json.loads(self.file(reference).read_text(encoding='utf-8')))
    def adapter(self):
        self.validate_overlap_policy()
        name,sep,cls=self.data['adapter'].partition(':')
        if not sep:raise InputError('adapter must be importable.module:Class')
        obj=getattr(importlib.import_module(name),cls)(self)
        return obj
    def validate_overlap_policy(self):
        if self.data['composition']['connected'] and self.data['adapter']!='mfobs.models.ring:Adapter':
            raise InputError('UNSOUND_OVERLAP_POLICY: connected=True omits admissible overlaps; exact synthesis requires connected=False. Only the digest-bound Ring8 adapter may use connected overlaps because every derived generator row is independently checked on its entire finite state space.')
    def specification(self):
        """Expose the manifest's Ω, rated native rules and selected raw counts."""
        from mfobs.api import ModelInputs, Observable, RatedRule
        from mfobs.state_space import from_manifest
        rules, queries = self.engine_rules(), self.engine_observables()
        return ModelInputs(
            omega=from_manifest(self),
            rules={row['id']: RatedRule(rules[row['id']], row['rate'], Fraction(row['per_match']))
                   for row in self.data['rules']},
            observables={name: Observable(name, queries[name])
                         for name in self.data['initial_observables']},
        )
    def engine_rules(self):
        from mfobs.core.engine import mf
        return {r['id']:mf.Rule.fromGMLString(self.file(r['file']).read_text(encoding='utf-8'),add=False) for r in self.data['rules']}
    def engine_observables(self):
        from mfobs.core.engine import mf
        return {r['id']:mf.Rule.fromGMLString(self.file(r['file']).read_text(encoding='utf-8'),add=False) for r in self.data['observables']}
    def forbidden_spec(self):
        if not self.data['forbidden']:return None
        from mfobs.core.forbidden import ForbiddenPatternSpecification,IncidentBound
        d=self.read('forbidden')
        required={'name','atomLabels','bondOrders','allowedBondPairs','maxValence','incidentBounds','strictSignature'}
        if not isinstance(d,dict) or set(d)!=required:
            raise InputError('FORBIDDEN_CONTRACT: unexpected fields in forbidden specification')
        try:
            return ForbiddenPatternSpecification(d['name'],d['atomLabels'],d['bondOrders'],
                d['allowedBondPairs'],d['maxValence'],
                tuple(IncidentBound(**v) for v in d['incidentBounds']),d['strictSignature'])
        except (TypeError,ValueError,KeyError) as exc:
            raise InputError(f'FORBIDDEN_CONTRACT: {exc}') from exc


def load(path):
    path=Path(path).resolve()
    try: data=json.loads(path.read_text(encoding='utf-8'))
    except (OSError,json.JSONDecodeError) as e:raise InputError(str(e)) from e
    required={'schema_version','id','title','adapter','rules','observables','initial_observables','composition','forbidden','state_language','reduction','normalisation','initial_state','parameters','reference'}
    optional={'specifications','initial_mixture'}
    if not required<=set(data) or set(data)-required-optional:raise InputError(f'Manifest keys: missing={sorted(required-set(data))}, unknown={sorted(set(data)-required-optional)}')
    if data['schema_version']!=1:raise InputError('Unsupported schema_version')
    identifier=lambda s:isinstance(s,str) and re.fullmatch(r'[A-Za-z_][A-Za-z_0-9]*',s)
    if not identifier(data['id']):raise InputError('Invalid model id')
    if not isinstance(data['rules'],list) or not data['rules']:raise InputError('Nonempty rules list required')
    for key,fields in [('rules',{'id','file','rate','per_match','kind'}),('observables',{'id','file'})]:
        if not isinstance(data[key],list) or not data[key]:raise InputError(f'Nonempty {key} list required')
        ids=[]
        for r in data[key]:
            if set(r)!=fields:raise InputError(f'Unexpected fields in {key}: {r}')
            if not identifier(r['id']):raise InputError(f'Invalid {key} id')
            ids.append(r['id'])
            if key=='rules':
                if not identifier(r['rate']):raise InputError('Rate names must be identifiers')
                if not isinstance(r['per_match'],str) or not re.fullmatch(r'[0-9]+(?:/[1-9][0-9]*)?',r['per_match']):raise InputError('per_match must be a quoted positive integer or fraction, e.g. 1 or 1/2')
                try:factor=Fraction(r['per_match'])
                except (ValueError,TypeError,ZeroDivisionError) as e:raise InputError('Invalid rational per_match multiplier') from e
                if factor<=0:raise InputError('per_match multiplier must be positive')
                if r['kind'] not in ['chemical','boundary','abstract']:raise InputError('Unknown rule kind')
        if len(ids)!=len(set(ids)):raise InputError(f'Duplicate id in {key}')
    if not data['initial_observables'] or len(set(data['initial_observables']))!=len(data['initial_observables']):raise InputError('Initial observables must be nonempty and unique')
    if set(data['initial_observables'])-set(r['id'] for r in data['observables']):raise InputError('Initial observable has no GML definition')
    comp=data['composition']
    if set(comp)!={'mode','connected','include_empty'}:raise InputError('Unexpected composition fields')
    if comp['mode'] not in ['public','restricted']:raise InputError('Unsupported composition mode')
    if any(type(comp[k]) is not bool for k in ['connected','include_empty']):raise InputError('Composition switches must be booleans')
    if comp['mode'].startswith('restricted') and not data['forbidden']:raise InputError('Restricted composition requires forbidden.json')
    obj=ModelManifest(path,data)
    for r in data['rules']+data['observables']:obj.file(r['file'])
    for key in ['state_language','reduction']:
        if data[key] is not None:obj.read(key)
    # Numerical values are not parsed during symbolic synthesis.
    if data['parameters'] is not None:obj.file(data['parameters'])
    if data['forbidden']:obj.forbidden_spec()
    if data.get('specifications'):
        if not isinstance(data['specifications'],dict) or set(data['specifications'])!={'chemistry','domain'}:
            raise InputError('specifications must contain exactly chemistry and domain')
        obj.specification_certificate()
    if data.get('initial_mixture'):obj.initial_mixture()
    if set(data['normalisation'])-set(r['id'] for r in data['observables']):raise InputError('Normalisation of undefined observable')
    for v in data['normalisation'].values():
        if set(v)!={'scale','meaning'} or Fraction(v['scale'])<=0:raise InputError('Invalid normalisation')
    obj.validate_overlap_policy()
    return obj


def with_observable_queries(inputs: ModelManifest, selection) -> ModelManifest:
    """Return an in-memory model whose selected observables come from files.

    This never edits the checked-in manifest.  The classic manifest worklist
    accepts one identity rule per coordinate; expression-valued queries are
    handled by the rule-expression worklist used by models such as Rémy.
    """

    if any(len(query.terms) != 1 or query.terms[0].coefficient != 1
           for query in selection.observables):
        raise InputError(
            "this model backend accepts one unit-coefficient GML identity per "
            "observable; repeated-name expressions require an expression backend"
        )
    data = deepcopy(inputs.data)
    data["observables"] = [
        {"id": query.name, "file": str(query.terms[0].file)}
        for query in selection.observables
    ]
    data["initial_observables"] = list(selection.names)
    data["normalisation"] = {
        query.name: data.get("normalisation", {}).get(query.name, {
            "scale": "1", "meaning": "raw injective identity-pattern count",
        })
        for query in selection.observables
    }
    overridden = ModelManifest(inputs.path, data)
    # Parse eagerly and require literal identity rules.  Backends subsequently
    # decide whether the graph shape belongs to their declared state language.
    from mfobs.core.engine import mf
    from mfobs.core.rules import is_identity
    for query in selection.observables:
        rule = mf.Rule.fromGMLString(query.terms[0].gml, add=False)
        if not is_identity(rule):
            raise InputError(f"{query.name}: selected observable is not an identity rule")
    return overridden
