"""Paper-example presets and their native input declarations."""
from __future__ import annotations

import mfobs.models.dna_walker.analysis as walker_analysis
import mfobs.models.dna_walker.inputs as walker_inputs
import mfobs.outputs as outputs
import mfobs.models.molecular_rop.inputs as molecular_inputs
import mfobs.models.remy.analysis as remy_analysis
import mfobs.models.aba.inputs as path_inputs

from dataclasses import dataclass
from pathlib import Path
from fractions import Fraction

import molfoundry as mf

from mfobs.api import ModelInputs, Observable, RatedRule
from mfobs.models.chlorination.analysis import Chlorination
from mfobs.state_space import from_manifest
from mfobs.inputs import ModelManifest


def _specification(omega, rule_gmls, query_gmls, rates, factors=None, query_factors=None):
    factors = factors or {}
    query_factors = query_factors or {}
    return ModelInputs(
        omega=omega,
        rules={name: RatedRule(mf.Rule.fromGMLString(gml, add=False), rates[name],
                              Fraction(factors.get(name, 1)))
               for name, gml in rule_gmls.items()},
        observables={name: Observable(name, terms=((Fraction(query_factors.get(name, 1)),
                                                  mf.Rule.fromGMLString(gml, add=False)),))
                     for name, gml in query_gmls.items()},
    )


class _ManifestInspection:
    """Rule and query GML of an application backed by a model manifest (``inputs()``)."""

    def specification(self) -> ModelInputs:
        return self.inputs().specification()

    def rule_gmls(self) -> dict[str, str]:
        model = self.inputs()
        return {row["id"]: model.file(row["file"]).read_text(encoding="utf-8") for row in model.data["rules"]}

    def query_gmls(self) -> dict[str, str]:
        model = self.inputs()
        return {row["id"]: model.file(row["file"]).read_text(encoding="utf-8") for row in model.data["observables"]}


@dataclass(frozen=True)
class RemyTrees:
    """Rémy's rooted planar full-binary-tree model."""

    order: int = 3

    def specification(self) -> ModelInputs:
        from mfobs.models.remy.inputs import state_family
        name = f"P{self.order}"
        rules = {key: RatedRule(mf.Rule.fromGMLString(gml, add=False), "rho")
                 for key, gml in self.rule_gmls().items()}
        terms = tuple((Fraction(1), mf.Rule.fromGMLString(gml, add=False))
                      for gml in self.query_gmls().values())
        return ModelInputs(state_family(), rules,
                           {name: Observable(name, terms=terms)})

    def run(
        self,
        inputs: ModelInputs,
        *,
        linear_limit: int = 16,
        record_patterns: bool = False,
    ) -> dict:
        inputs.require(self.specification())
        state_space = inputs.omega
        order = self.order
        result = remy_analysis.common_analysis(
            linear_limit=linear_limit,
            order=order,
            record_patterns=record_patterns,
            state_space=state_space,
        )
        result["observable_designation"] = {
            "selected": f"P{order}",
            "role": "post-input structural interpretation",
            "generic_core_identity": (
                "the worklist compares canonical sparse rule expressions; the P_k "
                "label is not used for equality or span membership"
            ),
        }
        return result

    def rule_gmls(self) -> dict[str, str]:
        from mfobs.models.remy import inputs as remy_trees

        return {
            f"G{side}_{stem}": remy_trees.insertion_rule_gml(side, stem)
            for side in ("L", "R")
            for stem in remy_trees.EDGE_TYPES
        }

    def query_gmls(self) -> dict[str, str]:
        from mfobs.models.remy import inputs as remy_trees

        return {
            stem: remy_trees.observable_gml(f"P{self.order}", stem)
            for stem in remy_trees.EDGE_TYPES
        }


@dataclass(frozen=True)
class ABA:
    """Growing directed paths with free tokens and selected word counts."""

    motifs: tuple[str, ...] = ("ABA",)
    alphabet: tuple[str, ...] = ("A", "B")
    gated: bool = False

    @classmethod
    def motif(cls, word, *, alphabet=("A", "B"), gated=False):
        return cls((word,), alphabet, gated)

    def specification(self) -> ModelInputs:
        omega = path_inputs.state_family(self.alphabet, self.gated)
        sources = path_inputs.source_rules(self.alphabet, self.gated, state_space=omega)
        patterns = [path_inputs.query_pattern(word) for word in self.motifs]
        return ModelInputs(omega, {r.identifier: r for r in sources},
                           {p.name: Observable(p.name, p.rule) for p in patterns})

    def rule_gmls(self):
        return self.specification().rule_gmls()

    def query_gmls(self):
        return {word: path_inputs.query_pattern(word).rule.getGMLString() for word in self.motifs}

    def run(self, inputs, **options):
        from mfobs.models.aba.analysis import analyze
        inputs.require(self.specification())
        return analyze(inputs, **options)

    def linear_closure(self, inputs, *, maximum_dimension=24):
        from mfobs.core.counting_linear_closure import discover_linear_observables
        from mfobs.models.aba.counting import PathCounter
        inputs.require(self.specification())
        return discover_linear_observables(tuple(inputs.rules.values()), inputs.observables,
            counter=PathCounter(state_space=inputs.omega), maximum_dimension=maximum_dimension)

    def jump_closure(self, inputs, **options):
        return self.run(inputs, run_linear=False, **options)


@dataclass(frozen=True)
class MolecularROP:
    """Atom-resolved lactone polymerisation on its covalent constructor family."""

    selection: object = "baseline"

    @classmethod
    def motif(cls, word, *, alphabet=("A", "B"), gated=False):
        return cls(molecular_inputs.configured_variant(motifs=(word,), alphabet=alphabet, gated=gated))

    def configuration(self):
        return molecular_inputs.configuration(self.selection)

    def specification(self) -> ModelInputs:
        selected = self.configuration()
        omega = molecular_inputs.state_family(selected.alphabet, selected.gated)
        sources = molecular_inputs.source_rules(selected, state_space=omega)
        words = ["".join(word) for word in selected.motifs]
        if selected in molecular_inputs.VARIANTS.values():
            words = [*selected.alphabet, *words]
        patterns = [molecular_inputs.query_pattern(word) for word in words]
        return ModelInputs(omega, {r.identifier: r for r in sources},
                           {p.name: Observable(p.name, p.rule) for p in patterns})

    def rule_gmls(self):
        return self.specification().rule_gmls()

    def query_gmls(self):
        selected = self.configuration()
        words = ["".join(word) for word in selected.motifs] or list(selected.alphabet)
        return {word: molecular_inputs.query_pattern(word).rule.getGMLString() for word in words}

    def run(self, inputs, **options):
        from mfobs.models.molecular_rop.analysis import analyze
        inputs.require(self.specification())
        return analyze(self.selection, state_space=inputs.omega, **options)

    def linear_closure(self, inputs, *, maximum_dimension=24, molecular_overlap_limit=2000000):
        from mfobs.core.counting_linear_closure import discover_linear_observables
        from mfobs.models.molecular_rop.counting import MolecularReducer
        inputs.require(self.specification())
        counter = MolecularReducer(self.configuration(), state_space=inputs.omega,
                                   maximum_candidates=molecular_overlap_limit)
        return discover_linear_observables(tuple(inputs.rules.values()), inputs.observables,
            counter=counter, maximum_dimension=maximum_dimension)

    def jump_closure(self, inputs, **options):
        return self.run(inputs, run_linear=False, **options)


@dataclass(frozen=True)
class EightSlots(_ManifestInspection):
    """Symbolic eight-site synthesis shared across the supported topologies."""

    model_path: Path | None = None

    def inputs(self):
        from mfobs.inputs import load

        root = Path(__file__).resolve().parents[2]
        model_path = (
            Path(self.model_path)
            if self.model_path is not None
            else root / "examples/eight-slots/model.json"
        )
        return load(model_path)

    def run(
        self,
        inputs: ModelInputs,
        *,
        max_expansions: int = 100,
        output: Path | None = None,
        checkpoint: Path | None = None,
        progress: bool = False,
    ) -> dict:
        inputs.require(self.specification())
        state_space = inputs.omega
        from mfobs.core.worklist import run
        from mfobs.outputs import export

        output = Path(output) if output is not None else None
        if output is not None:
            output.mkdir(parents=True, exist_ok=True)
        if checkpoint is None and output is not None:
            checkpoint = output / "checkpoint.json"

        from mfobs.provenance import source_hashes, runtime_fingerprint
        sources = source_hashes()
        runtime = runtime_fingerprint()

        model = self.inputs()
        state_space.require(from_manifest(model))
        result = run(
            model,
            progress=progress,
            max_expansions=max_expansions,
            checkpoint=checkpoint,
        )
        result["source_sha256"] = sources
        result["runtime"] = runtime
        if output is not None:
            export(result, output)
        return result


@dataclass(frozen=True)
class Ring8(_ManifestInspection):
    """Complete eight-site ring quotient with composition-derived observables."""

    manifest: ModelManifest | None = None

    def inputs(self):
        if self.manifest is not None:
            return self.manifest
        from mfobs.inputs import load
        from mfobs.models import ring8_common

        return load(ring8_common.MODEL_PATH)

    def run(
        self, inputs: ModelInputs, *, worklist_limit: int = 1000,
        include_evidence: bool = False,
    ) -> dict:
        inputs.require(self.specification())
        state_space = inputs.omega
        from mfobs.core.worklist import run as discover_linear
        from mfobs.models import ring8_common

        model = self.inputs()
        state_space.require(from_manifest(model))
        linear_evidence = discover_linear(model, max_expansions=worklist_limit)
        result = ring8_common.analyze(
            composition_result=linear_evidence,
            inputs=model,
        )
        if include_evidence:
            result["linear_evidence"] = linear_evidence
        return result


@dataclass(frozen=True)
class BrusselatorROP:
    """Brusselator molecular rules coupled to the lactone-ROP coordinate."""

    selected_queries: tuple[str, ...] | None = None

    def specification(self) -> ModelInputs:
        from mfobs.models import brusselator_growth_chemistry as chemistry
        contract = self.contract()
        from mfobs.core.molecules import automorphism_count
        from mfobs.models.molecular_rop import inputs as copolymer_rop
        query_factors = {
            "X": Fraction(1, automorphism_count(chemistry.species("X"))),
            "Y": Fraction(1, automorphism_count(chemistry.species("Y"))),
            "M": Fraction(1, automorphism_count(copolymer_rop.monomer("B"))),
            "N": Fraction(1),
        }
        return _specification(
            chemistry.state_family(), self.rule_gmls(), self.query_gmls(),
            {rule.identifier: rule.rate for rule in contract.rules},
            chemistry.rate_factors(), query_factors,
        )

    def contract(self):
        from mfobs.models import brusselator_pjc

        return brusselator_pjc.analysis_contract()

    def rule_gmls(self) -> dict[str, str]:
        contract = self.contract()
        return {rule.identifier: rule.gml for rule in contract.rules}

    def coordinate_gmls(self) -> dict[str, str]:
        """All declared population/support coordinate identities."""
        from mfobs.models import brusselator_pjc

        names = (*brusselator_pjc.SELECTED_QUERIES, *brusselator_pjc.SUPPORT_COORDINATES)
        return {
            name: brusselator_pjc.coordinate_identity_gml(name)
            for name in names
        }

    def query_gmls(self) -> dict[str, str]:
        selected = (
            self.selected_queries
            if self.selected_queries is not None
            else tuple(self.coordinate_gmls())[:2]
        )
        all_coordinates = self.coordinate_gmls()
        return {name: all_coordinates[name] for name in selected}

    def run(self, inputs: ModelInputs, *, linear_limit: int = 24) -> dict:
        inputs.require(self.specification())
        state_space = inputs.omega
        from mfobs.models import brusselator_pjc

        selected = tuple(inputs.observables)
        return brusselator_pjc.analyze(
            linear_limit=linear_limit,
            selected_queries=selected,
            state_space=state_space,
        )

    def linear_closure(self, inputs: ModelInputs, *, maximum_dimension=24) -> dict:
        from mfobs.models.brusselator_pjc import linear_closure
        inputs.require(self.specification())
        return linear_closure(tuple(inputs.observables), linear_limit=maximum_dimension)

    def jump_closure(self, inputs: ModelInputs) -> dict:
        from mfobs.models.brusselator_pjc import analyze
        inputs.require(self.specification())
        return analyze(selected_queries=tuple(inputs.observables), state_space=inputs.omega,
                       run_linear=False)


@dataclass(frozen=True)
class GrowingRing(_ManifestInspection):
    """Unbounded structural-NAC ring with an honest bounded exact search."""

    manifest: ModelManifest | None = None

    def inputs(self):
        if self.manifest is not None:
            return self.manifest
        from mfobs.inputs import load
        from mfobs.models import structural_nac_common

        return load(structural_nac_common.EXAMPLE_ROOT / "model.json")

    def run(
        self, inputs: ModelInputs, *, worklist_limit: int = 80,
        include_evidence: bool = False,
    ) -> dict:
        inputs.require(self.specification())
        state_space = inputs.omega
        from mfobs.core.worklist import run as discover_linear
        from mfobs.models import structural_nac_common

        model = self.inputs()
        state_space.require(from_manifest(model))
        linear_evidence = discover_linear(
            model,
            max_expansions=worklist_limit,
        )
        profile = {
            "budget": worklist_limit,
            "status": linear_evidence["status"],
            "basis_count": len(linear_evidence["basis"]),
            "pending_count": len(linear_evidence["pending"]),
            "maximum_word_length": linear_evidence[
                "maximum_generated_word_length"
            ],
            "method": (
                "fresh exact local-word FIFO worklist from "
                + ",".join(model.data["initial_observables"])
            ),
            "interpretation": (
                "bounded outcome with retained pending frontier; "
                "not a nonclosure proof"
            ),
        }
        result = structural_nac_common.analyze(
            linear_result=linear_evidence,
            profile=profile,
            inputs=model,
        )
        if include_evidence:
            result["linear_evidence"] = linear_evidence
        return result


@dataclass(frozen=True)
class DNAWalker:
    """Two-legged walker on an arbitrary directed cycle of length at least five."""

    def inputs(self):
        import json
        source = Path(__file__).resolve().parents[2] / "examples/dna-walker/input.json"
        return json.loads(source.read_text(encoding="utf-8"))

    def specification(self) -> ModelInputs:
        data = self.inputs()
        return ModelInputs(
            walker_inputs.state_family(),
            {r["name"]: RatedRule(walker_inputs.convert_rule(r), r["name"])
             for r in data["rules"]},
            {"E": Observable("E", walker_inputs.identity(data["query"], "E"))},
        )

    def run(self, inputs: ModelInputs, *, maximum_dimension=24):
        inputs.require(self.specification())
        state_space = inputs.omega
        data = self.inputs()
        result = walker_analysis.analyze(data["rules"], data["query"], maximum_dimension=maximum_dimension, state_space=state_space)
        result.update(outputs.dna_walker_readout(result, data["rules"]))
        return result


__all__ = [
    "BrusselatorROP",
    "Chlorination",
    "GrowingRing",
    "EightSlots",
    "DNAWalker",
    "ABA",
    "MolecularROP",
    "RemyTrees",
    "Ring8",
]
