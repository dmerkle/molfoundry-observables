"""Architecture guard for maintained public entry points.

The goal is not to ban model modules from notebooks: they are still useful for
input drawings, concrete-state constructors and independent checks.  The guard
only prevents maintained analysis entry points from bypassing the public API
and wiring the internal worklists/compilers directly again.
"""
from pathlib import Path
import ast

ROOT = Path(__file__).resolve().parents[1]

ACTIVE_ANALYSIS_SOURCES = (
    ROOT / "notebooks/examples/02_aba.py",
    ROOT / "notebooks/examples/03_chlorination.py",
    ROOT / "notebooks/examples/04_remy.py",
    ROOT / "notebooks/examples/05_eight_slots.py",
    ROOT / "notebooks/examples/07_molecular_rop.py",
    ROOT / "notebooks/examples/08_query_variation.py",
    ROOT / "notebooks/examples/10_brusselator.py",
    ROOT / "notebooks/examples/11_growing_ring.py",
)

SUPPORTING_ENTRYPOINTS = (
    ROOT / "scripts/run_eight_slots.py",
    ROOT / "src/mfobs/cli.py",
)

EXAMPLE_RUNNERS = (
    ROOT / "examples/remy/run.py",
    ROOT / "examples/aba/run.py",
    ROOT / "examples/ring8/run.py",
    ROOT / "examples/chlorination/run.py",
    ROOT / "examples/brusselator/run.py",
    ROOT / "examples/growing_ring/run.py",
)

FORBIDDEN_ANALYSIS_BYPASSES = (
    "from mfobs.core.worklist import run",
    "common_analysis(",
    "copolymer_rop.analyze(",
    "rop.analyze(",
    "brusselator_pjc.analyze(",
    "structural_nac_common.analyze(",
    "from examples.ring8.run import",
    "from examples.chlorination.run import",
    "from examples.remy.run import",
    "from examples.brusselator.run import",
    "from examples.growing_ring.run import",
    '"run.py", "rop-analysis"',
)


def test_maintained_notebook_analysis_enters_through_public_api():
    for path in ACTIVE_ANALYSIS_SOURCES:
        source = path.read_text(encoding="utf-8")
        assert "from mfobs import" in source, path
        assert "analyze" in source or (".linear_closure(" in source and ".jump_closure(" in source), path
        for forbidden in FORBIDDEN_ANALYSIS_BYPASSES:
            assert forbidden not in source, (path, forbidden)


def test_artifact_runners_enter_through_public_api():
    for path in EXAMPLE_RUNNERS:
        source = path.read_text(encoding="utf-8")
        assert "from mfobs import" in source, path
        assert "analyze" in source, path
        assert "from mfobs.core.worklist import run" not in source, path
        assert "common_analysis(" not in source, path


def test_tcs_is_intentionally_not_forced_through_closure_api():
    source = (ROOT / "examples/tcs/replay.py").read_text(encoding="utf-8")
    assert "rcCommon" in source
    assert "from mfobs import" not in source


def test_supporting_entrypoints_do_not_reintroduce_example_specific_analysis_bypasses():
    eight_slots = (ROOT / "scripts/run_eight_slots.py").read_text(encoding="utf-8")
    assert "from mfobs import EightSlots, analyze" in eight_slots
    assert "from mfobs.core.worklist import run" not in eight_slots

    cli = (ROOT / "src/mfobs/cli.py").read_text(encoding="utf-8")
    rop_section = cli.split("def rop_analysis_command", 1)[1].split(
        "def main", 1
    )[0]
    assert "from mfobs.api import analyze as analyze_model" in rop_section
    assert "MolecularROP" in rop_section
    assert "copolymer_rop.analyze(" not in rop_section


def test_active_notebooks_do_not_reintroduce_model_specific_input_plumbing():
    forbidden_by_file = {
        "02_aba.py": ("terminal_paths",),
        "03_chlorination.py": ("from mfobs.inputs import load",),
        "05_eight_slots.py": ("from mfobs.inputs import load",),
        "07_molecular_rop.py": ("molecular_rop_discovery",),
        "08_query_variation.py": ("terminal_paths",),
        "10_brusselator.py": ("brusselator_pjc",),
        "11_growing_ring.py": ("structural_nac_common", "from mfobs.inputs import load"),
    }
    root = ROOT / "notebooks/examples"
    for name, forbidden in forbidden_by_file.items():
        # Inspect executable Python, not Markdown links to the implementation.
        source = ast.unparse(ast.parse((root / name).read_text(encoding="utf-8")))
        for token in forbidden:
            assert token not in source, (name, token)
