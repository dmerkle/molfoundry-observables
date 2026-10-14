"""Independent population reference for the composition-derived Brusselator model."""
from functools import lru_cache
from mfobs.core.polynomial_jump_closure import (
    JumpEvent, JumpModelProvenance, PolynomialJumpModel, PolynomialRing,
)
from mfobs.models import brusselator_growth_chemistry as chemistry
from mfobs.models.brusselator_pjc import COORDINATES, SELECTED_QUERIES, analysis_contract


@lru_cache(None)
def semantic_oracle_model(
    selected_queries: tuple[str, ...] = SELECTED_QUERIES,
) -> PolynomialJumpModel:
    """Independent population-language oracle; never a production input."""

    contract = analysis_contract(selected_queries)
    ring = PolynomialRing(COORDINATES)
    x, y, monomer = (ring.coordinate(name) for name in ("X", "Y", "M"))
    falling_x2_y = ring.multiply(
        x, ring.subtract(x, ring.constant(1)), y
    )
    events = (
        JumpEvent("source", chemistry.RATE_CHANNELS[chemistry.SOURCE],
                  (1, 0, 0, 0), ring.constant(1),
                  chemistry.SOURCE, "one chemostatted source event"),
        JumpEvent("autocatalysis", chemistry.RATE_CHANNELS[chemistry.AUTOCATALYSIS],
                  (1, -1, 0, 0), falling_x2_y,
                  chemistry.AUTOCATALYSIS, "two distinct X and one Y"),
        JumpEvent("conversion", chemistry.RATE_CHANNELS[chemistry.CONVERSION],
                  (-1, 1, 0, 0), x,
                  chemistry.CONVERSION, "one X; B chemostatted"),
        JumpEvent("sink", chemistry.RATE_CHANNELS[chemistry.SINK],
                  (-1, 0, 0, 0), x,
                  chemistry.SINK, "one X"),
        JumpEvent("monomer-inlet", chemistry.RATE_CHANNELS[chemistry.MONOMER_INLET],
                  (0, 0, 1, 0), ring.constant(1),
                  chemistry.MONOMER_INLET, "one persistent boundary anchor"),
        JumpEvent("growth", chemistry.RATE_CHANNELS[chemistry.GROWTH],
                  (0, 0, -1, 1), ring.multiply(x, monomer),
                  chemistry.GROWTH, "one X, one free lactone and the unique chain end"),
    )
    certificate = {
        "status": "VERIFIED_ALL_STATES",
        "method": "literal stoichiometry plus exact molecular-match normalization",
        "argument": (
            "On every state in the declared language each rule has the recorded "
            "constant population jump. Complete-component matching gives X, Y and M; "
            "injectivity gives (X)_2; the unique chain end gives the growth factor X*M."
        ),
        "enumeration_used": False,
    }
    return PolynomialJumpModel(
        contract.model_id,
        ring,
        events,
        certificate,
        contract.digest,
        JumpModelProvenance(
            backend="independent population-language validation oracle",
            event_partition_origin="literal reaction channels",
            jump_origin="explicit oracle stoichiometry",
            activity_origin="exact molecular component and falling-factorial count",
            jump_projection="semantic population transition",
            composition_scope="none; independent validation oracle",
            all_jumps_composition_derived=False,
            supplied_delta_table=True,
        ),
    )

