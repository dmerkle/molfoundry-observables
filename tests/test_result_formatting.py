"""The human-readable view must preserve the mathematical status of a record."""
from __future__ import annotations

from mfobs.formatting import render


def test_nac_chlorination_prints_its_generated_events(chlorination_records):
    record = chlorination_records["guarded"]
    view = render(record, name="chlorination-guarded")
    assert "Selected query: Z" in view and "Condition: generic-not-exists" in view
    assert "Coordinates: Z, P1, P2, P3" in view
    assert "Event classes: 4" in view and "event-1" in view
    assert "Jump Δ" in view and "Rate × activity" in view
    assert "No finite event table" not in view
    assert "not a closed linear expectation ODE" in render(record, section="ode")


def test_aba_does_not_mislabel_first_moment_identities_as_closed_ode(aba_records):
    record = aba_records["aba"]
    view = render(record, name="aba")
    assert "Selected queries: C_ABA" in view
    assert "No finite closed linear ODE established" in view
    assert "undetermined; dimension budget reached" in view
    linear, events = view.split("Finite events and EMGF (Algorithm 2)")
    assert "d⟨" not in linear
    assert "Exact first-moment identities from these events" in events
    assert "d⟨C_ABA⟩/dt = k_A ⟨M_A S_AB⟩" in events
    assert "Event classes: 8" in view
    assert "Maximum differential order: 2" in view
    assert "beta_A × 1" in view
    assert "beta_A 1" not in view


def test_worker_wrapper_and_section_selection(aba_records):
    analysis = aba_records["aba"]
    record = {"analysis": analysis}
    ode = render(record, section="ode")
    jumps = render(record, section="jumps")
    assert "No finite closed linear ODE established" in ode
    assert "Event classes" not in ode
    assert "Exact identities computed by Algorithm 1" in ode
    assert "d⟨M_A*S_AB⟩/dt" in ode
    assert ode == render({"linear": analysis["linear"], "contract": analysis["contract"]}, section="ode")
    assert "Event classes: 8" in jumps
    assert "d⟨C_ABA⟩" not in jumps
    operator = render(record, section="operator")
    assert "Example: Growing paths" in operator
    assert "Explicit terms of D" in operator
    assert "feed_A-support-0: beta_A [exp(λ_M_A) - 1] (1)" in operator
    assert "∂_{λ_M_A} ∂_{λ_S_AB}" in operator


def test_local_tcs_claims_no_dynamics():
    view = render({"status": "PASS", "host_witness": {}}, name="tcs")
    assert "no ODE or jump closure is claimed" in view


def test_growing_ring_bounded_search_is_not_reported_as_a_result():
    record = {
        "linear": {"status": "UNKNOWN_EXPANSION_LIMIT", "pending": ["next"]},
        "pjc": {"status": "UNKNOWN_UNBOUNDED_COORDINATE_LIFT", "operator": None},
    }
    view = render(record, name="growing-ring")
    assert "undetermined; expansion budget reached" in view
    assert "not certified for the selected coordinates" in view
    assert "No finite event table or EMGF operator" in view


def test_missing_input_does_not_claim_a_search_ran():
    record = {
        "linear": {"status": "UNKNOWN_PREREQUISITE"},
        "pjc": {"status": "UNKNOWN_PREREQUISITE"},
    }
    view = render(record)
    assert view.count("not run; required input unavailable") == 2
    assert "Exact finite ODE" not in view


def test_remy_pjc_only_record_is_inspectable(remy_record):
    pjc = remy_record["pjc"]
    view = render(pjc, name="remy", section="jumps")
    assert "Event classes: 4" in view
    assert "J0" in view and "J3" in view
    assert "d⟨P3⟩" not in view


def test_ring8_branches_produce_closed_ode_and_jumps(ring_composition):
    from mfobs.models import ring8_common
    linear_result = ring_composition
    pjc = ring8_common.analyze(composition_result=linear_result)["pjc"]
    record = {"linear": {"status": linear_result["algorithm1_status"],
                         "basis": linear_result["basis"],
                         "steps": linear_result["algorithm1_steps"]}, "pjc": pjc}
    view = render(record, name="ring8")
    assert "Exact finite ODE" in view
    assert "Event classes:" in view
    assert "d⟨P0⟩/dt" in view
    assert "use --section ode for all" in view
    assert "use --section jumps for all" in view
    assert "basal_bind-jump-1" not in view
    assert "d⟨P29⟩/dt" in render(record, section="ode")
    assert render(record, section="ode") == render({"linear": record["linear"]}, section="ode")


def test_incomplete_rate_channels_are_not_printed_as_a_complete_derivative():
    record = {"linear": {
        "status": "UNKNOWN_DIMENSION_LIMIT", "basis": ["X"],
        "steps": [{"source": "X", "channel": "birth", "coordinates": {"X": "1"}}],
        "obligations": [{"source": "X", "channel": "death"}],
    }}
    assert "d⟨X⟩/dt" not in render(record, section="ode")


def test_flat_exact_linear_record_uses_the_recorded_equations(ring_composition):
    record = ring_composition
    view = render(record, name='Eight sites')
    assert 'Selected queries: P0' in view
    assert 'Exact finite ODE; 30 retained observables' in view
    assert 'dmu_P0/dt =' in view
    assert 'Statistics: compositions:' in view
    assert len(view.splitlines()) < 30
    full = render(record, section='ode')
    assert 'dmu_P29/dt =' in full


def test_terminal_definitions_show_support_meanings_without_gml(aba_records):
    record = aba_records["aba"]
    view = render(record, name='ABA')
    assert 'C_ABA: occurrences in the oriented chain backbone' in view
    assert 'M_A: number of free monomer copies' in view
    assert 'rule [' not in view
