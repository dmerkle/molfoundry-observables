"""Export the atom-resolved ROP analysis for the ABA motif."""
from mfobs import MolecularROP, analyze


def run(*, queries=None):
    if queries is not None:
        raise ValueError("Choose the molecular motif in tutorials/rop_analysis.py")
    backend = MolecularROP.motif("ABA")
    return analyze(backend.specification(), backend=backend, linear_limit=24).record()
