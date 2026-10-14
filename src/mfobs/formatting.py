"""Readable view of a computed paper record; never reruns or changes an analysis.

From the implementation root: PYTHONPATH=src python -m mfobs.formatting result.json
The JSON remains the complete derivation record; this view prints its mathematical
claims without exposing backend diagnostics as if they were equations.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from fractions import Fraction
import json
from pathlib import Path
from tabulate import tabulate


PREVIEW_EQUATIONS = 6
PREVIEW_EVENTS = 8


def _analysis(record: dict) -> dict:
    if isinstance(record.get("analysis"), dict):
        record = record["analysis"]
    if "operator" in record and "observable_coordinates" in record:
        return {"pjc": record}
    if "basis" in record and "status" in record and "closed" not in record and "linear" not in record:
        record = {**record, "linear": record}
    linear = record.get("linear") or {}
    basis = linear.get("basis", [])
    if isinstance(basis, list) and basis and isinstance(basis[0], dict):
        record = {**record, "linear": {**linear, "basis": [row["name"] for row in basis],
                                       "origins": {row["name"]: row for row in basis}}}
    return record


def observable_rows(record: dict) -> list[dict]:
    """Read definitions already present in a result, without reconstructing counts."""
    result = _analysis(record)
    linear = result.get("linear") or result
    pjc = result.get("pjc") or {}
    support = result.get("support_discovery") or result.get("synthesis", {}).get("support", {})
    patterns = {**result.get("identity_patterns", {}), **support.get("patterns", {})}
    expressions = {**linear.get("observable_expressions", {}),
                   **result.get("observable_expressions", {})}
    for name, expression in linear.get("raw_rule_expressions", {}).items():
        if expression.get("gml_terms"):
            expressions.setdefault(name, expression["gml_terms"])
    if isinstance(linear.get("basis"), dict) and linear.get("patterns"):
        expressions.update({name: [{"coefficient": str(value), "gml": linear["patterns"][key]}
                                  for key, value in expression.items()]
                            for name, expression in linear["basis"].items()})
    meanings = {row["name"]: row.get("meaning", "")
                for row in result.get("contract", {}).get("observables", [])}
    names = dict.fromkeys([*linear.get("basis", []), *patterns, *expressions,
                           *pjc.get("observable_coordinates", []),
                           *pjc.get("coordinates", [])])
    rows = []
    for name in names:
        pattern = patterns.get(name, {})
        meaning = (meanings.get(name) or pattern.get("meaning")
                   or linear.get("origins", {}).get(name, {}).get("meaning"))
        if meaning and (str(meaning).startswith("{") or "Fraction(" in str(meaning)):
            meaning = None
        if not meaning:
            meaning = " ".join(str(pattern.get(key, "")) for key in ("kind", "word")).strip()
        terms = expressions.get(name)
        if terms is None and pattern.get("gml"):
            terms = [{"coefficient": "1", "gml": pattern["gml"]}]
        text = linear.get("basis_expressions", {}).get(name, {}).get("text")
        if text == name and name.isidentifier():
            text = None
        rows.append({"name": name, "meaning": meaning or "", "terms": terms or [], "text": text})
    return rows


def _preview(lines: list[str], record: dict) -> None:
    rows = observable_rows(record)
    if rows:
        lines.append(f"Discovered observables: {len(rows)}")
        descriptive = [row for row in rows if row["meaning"]]
        for row in descriptive[:PREVIEW_EQUATIONS]:
            description = row["meaning"] or row["text"]
            if description:
                lines.append(f"  {row['name']}: {description[:180]}")
    stats = []
    record = {**(record.get("linear") or {}), **record}
    if "raw_compositions" not in record and "raw_composites" in record:
        record["raw_compositions"] = record["raw_composites"]
    for key, label in (("raw_compositions", "compositions"),
                       ("ordered_composition_calls", "ordered calls")):
        if key in record:
            stats.append(f"{label}: {record[key]}")
    if isinstance(record.get("seconds"), (int, float)):
        stats.append(f"analysis: {record['seconds']:.3f} s")
    if stats:
        lines.append("Statistics: " + "; ".join(stats))


def _flat_linear(record: dict, name: str | None, section: str) -> str:
    from mfobs.outputs import ode_text

    if section in {"jumps", "operator"}:
        return f"Example: {name or record.get('title') or record['model']}\nNo finite event table or EMGF operator is established in this record."
    lines = [f"Example: {name or record.get('title') or record['model']}",
             "Selected queries: " + ", ".join(record.get("initial_observables", [])),
             "\nLinear expectations (Algorithm 1)", "Status: " + _status(record["status"])]
    if record["closed"]:
        lines.append(f"Exact finite ODE; {len(record['basis'])} retained observables.")
        equations = [line for line in ode_text(record).splitlines() if line.startswith("dmu_")]
        shown = equations if section == "ode" else equations[:3]
        lines.extend("  " + (line if section == "ode" else line[:240] + (" …" if len(line) > 240 else ""))
                     for line in shown)
        if len(shown) < len(equations):
            lines.append(f"  … {len(equations) - len(shown)} more equations.")
    else:
        lines.append("No finite closed linear ODE established by this run.")
    if section == "all":
        _preview(lines, record)
    return "\n".join(lines)


def _monomial(factors: dict) -> str:
    if not factors:
        return "1"
    return " ".join(
        name if power == 1 else f"{name}^{power}"
        for name, power in sorted(factors.items())
    )


def _polynomial(poly: dict, variable) -> str:
    """Render exact recorded coefficients with a chosen monomial notation."""
    pieces = []
    for term in poly.get("terms", []):
        coefficient = Fraction(term["coefficient"])
        if not coefficient:
            continue
        monomial = variable(term["monomial"])
        magnitude = abs(coefficient)
        body = monomial if magnitude == 1 else f"{magnitude} {monomial}"
        pieces.append(("-" if coefficient < 0 else "+", body))
    if not pieces:
        return "0"
    first_sign, first_body = pieces[0]
    return ("-" if first_sign == "-" else "") + first_body + "".join(
        f" {sign} {body}" for sign, body in pieces[1:]
    )


def _expectation(poly: dict) -> str:
    """Expectation of a recorded generator polynomial, not a fitted ODE."""
    return _polynomial(poly, lambda factors: f"⟨{_monomial(factors)}⟩" if factors else "1")


def _differential(poly: dict) -> str:
    def derivative(factors: dict) -> str:
        return " ".join(
            f"∂_{{λ_{name}}}" + (f"^{power}" if power != 1 else "")
            for name, power in sorted(factors.items())
        ) or "1"

    return _polynomial(poly, derivative)


def _linear_images(linear: dict) -> dict:
    images = defaultdict(dict)
    for step in linear.get("steps", []):
        if "generator_image" in step:
            images[step["source"]][step["channel"]] = step["generator_image"]
        elif "coordinates" in step:
            images[step["source"]][step["channel"]] = {"terms": [
                {"coefficient": str(value), "monomial": {name: 1}}
                for name, value in step["coordinates"].items() if Fraction(value)]}
    channels = {step["channel"] for step in linear.get("steps", [])}
    channels.update(row["channel"] for row in linear.get("obligations", []))
    # A stopped worklist may have computed only some channels for its last source.
    return {source: rows for source, rows in images.items() if channels <= rows.keys()}


def _rhs(channels: dict, rate_symbols: dict[str, str]) -> str:
    summands = []
    for channel, polynomial in channels.items():
        terms = polynomial.get("terms", [])
        if not terms:
            continue
        rate = rate_symbols.get(channel, channel)
        if len(terms) > 1:
            summands.append(f"{rate} ({_expectation(polynomial)})")
        else:
            term = terms[0]
            coefficient = Fraction(term["coefficient"])
            if not coefficient:
                continue
            factors = term["monomial"]
            parts = ([] if abs(coefficient) == 1 else [str(abs(coefficient))]) + [rate]
            if factors:
                parts.append(f"⟨{_monomial(factors)}⟩")
            summands.append(("-" if coefficient < 0 else "") + " ".join(parts))
    return " + ".join(summands).replace(" + -", " - ") or "0"


def _jump(delta: dict) -> str:
    if not delta:
        return "0"
    return ", ".join(f"{name}:{amount:+d}" for name, amount in sorted(delta.items()))


def _event_table(rows) -> str:
    return tabulate(rows, headers=("Event class", "Jump Δ", "Rate × activity"),
                    tablefmt="simple", disable_numparse=True,
                    colalign=("left", "left", "left"), maxcolwidths=[None, 36, 48])


def _status(status: str) -> str:
    return {
        "CLOSED": "closed",
        "UNKNOWN_DIMENSION_LIMIT": "undetermined; dimension budget reached",
        "UNKNOWN_EXPANSION_LIMIT": "undetermined; expansion budget reached",
        "UNKNOWN_UNBOUNDED_COORDINATE_LIFT": "not certified for the selected coordinates",
        "UNKNOWN_PREREQUISITE": "not run; required input unavailable",
    }.get(status, status.replace("_", " ").lower())


def _names(values: list[str], *, full: bool) -> str:
    shown = values if full or len(values) <= 8 else values[:PREVIEW_EQUATIONS]
    suffix = f", … ({len(values)} total)" if len(shown) < len(values) else ""
    return ", ".join(shown) + suffix


def render(record: dict, *, name: str | None = None, section: str = "all") -> str:
    """Format only claims actually present in the supplied result schema."""
    if section not in {"all", "ode", "jumps", "operator"}:
        raise ValueError(f"unknown section {section!r}")
    result = _analysis(record)
    for group in ("variants", "queries"):
        if group in result:
            if section == "all":
                lines = [f"Example: {name or group}", "Calculation | linear closure | pJC | counts | events"]
                for label, child in result[group].items():
                    child = _analysis(child)
                    linear, pjc = child.get("linear") or {}, child.get("pjc") or {}
                    support = child.get("support_discovery") or child.get("synthesis", {}).get("support", {})
                    coordinates = list(support.get("patterns", {})) or pjc.get("observable_coordinates") or pjc.get("coordinates", [])
                    count = pjc.get("event_class_count", pjc.get("event_count"))
                    if count is None:
                        terms = (pjc.get("operator") or {}).get("terms")
                        count = len(terms) if terms is not None else len(child["events"]) if "events" in child else "—"
                    lines.append(f"  {label} | {_status(linear.get('status', 'not run'))} | "
                                 f"{_status(pjc.get('status', 'not run'))} | {_names(coordinates, full=False)} | {count}")
                return "\n".join(lines)
            return "\n\n".join(render(child, name=child_name, section=section)
                                 for child_name, child in result[group].items())
    if "closed" in result and "basis" in result:
        return _flat_linear(result, name, section)
    title = (name or result.get("contract", {}).get("title")
             or result.get("configuration", {}).get("model", {}).get("title")
             or result.get("variant", "paper analysis"))
    lines = [f"Example: {title}"]
    if "measurements" in result:
        lines.append("Word | MolFoundry seconds | FRAGGER seconds")
        lines.extend(f"  {row['word']} | {row['molfoundry_median_seconds']} | {row['fragger_median_seconds']}"
                     for row in result["measurements"])
        return "\n".join(lines)
    if "compositions" in result:
        lines.append("Conditional compositions; no ODE or jump closure is claimed.")
        lines.extend(f"  {name}: {len(terms)} signed terms"
                     for name, terms in result["compositions"].items())
        return "\n".join(lines)
    if "events" in result and "synthesis" in result:
        lines += ["Selected query: " + result["synthesis"]["selected_query"],
                  "Condition: " + result["rule"]["condition_kind"],
                  "Polynomial jump closure: " + _status(result["pjc"]["status"]),
                  "Coordinates: " + ", ".join(result["pjc"]["coordinates"])]
        if section == "ode":
            lines.append("This run computes an event model, not a closed linear expectation ODE.")
        else:
            lines.append(f"Event classes: {len(result['events'])}")
            if result["events"]:
                lines.append(_event_table([
                    [event["id"], _jump(event["delta"]), f"k × ({event['activity']['text']})"]
                    for event in result["events"]]))
            if section == "operator":
                lines.append("D = k Σ_e (exp(λ·Δ_e) − 1) a_e(∂λ)")
        if section == "all":
            _preview(lines, result)
        return "\n".join(lines)
    selected = [observable["name"] for observable in
                result.get("contract", {}).get("observables", [])
                if observable.get("role") == "selected"]
    if not selected:
        selected = result.get("configuration", {}).get("selected_queries", [])
    if not selected:
        selected = result.get("support_discovery", {}).get("selected_queries", [])
    if not selected:
        selected = (result.get("linear") or {}).get("seeds", [])
    if not selected:
        selected = [name for name, origin in (result.get("linear") or {}).get("origins", {}).items()
                    if origin.get("origin") == "input-seed"]
    if selected:
        lines.append("Selected queries: " + ", ".join(selected))
    if "linear" not in result and "pjc" not in result:
        if result.get("status") == "PASS" and "host_witness" in result:
            lines.append("Local conditional composition: PASS; no ODE or jump closure is claimed.")
        else:
            lines.append("This record contains no linear or jump-closure result.")
        return "\n".join(lines)

    linear = result.get("linear") or {}
    pjc = result.get("pjc") or {}
    rate_symbols = {
        rule["id"]: rule["rate"]
        for rule in result.get("contract", {}).get("normalized_rules", [])
    }
    if section in {"all", "ode"}:
        lines.append("\nLinear expectations (Algorithm 1)")
        status = linear.get("status", "not run")
        lines.append(f"Status: {_status(status)}")
        if status == "CLOSED":
            if linear.get("basis"):
                basis = list(linear["basis"])
                lines.append("Retained basis: " + _names(basis, full=section == "ode"))
            lines.append("Exact finite ODE for expectations of the retained basis:")
        else:
            lines.append("No finite closed linear ODE established by this run.")
            pending = linear.get("pending") or []
            if pending:
                lines.append(f"Pending work: {len(pending)} direction(s).")
        images = _linear_images(linear)
        if images and (status == "CLOSED" or section == "ode"):
            if status != "CLOSED":
                lines.append("Exact identities computed by Algorithm 1 (not a closed ODE):")
            coordinates = [coordinate for coordinate in (linear.get("basis") or images)
                           if coordinate in images or (coordinate == "1" and status == "CLOSED")]
            limit = 3 if len(coordinates) > 20 else PREVIEW_EQUATIONS
            shown = coordinates if section == "ode" else coordinates[:limit]
            for coordinate in shown:
                if coordinate in images:
                    lines.append(f"  d⟨{coordinate}⟩/dt = {_rhs(images[coordinate], rate_symbols)}")
                elif coordinate == "1" and status == "CLOSED":
                    lines.append("  d⟨1⟩/dt = 0")
            if len(shown) < len(coordinates):
                lines.append(f"  … {len(coordinates) - len(shown)} more equations; use --section ode for all.")

    if section in {"all", "jumps", "operator"}:
        lines.append("\nFinite events and EMGF (Algorithm 2)")
        lines.append(f"Status: {_status(pjc.get('status', 'not run'))}")
        operator = pjc.get("operator") or {}
        terms = operator.get("terms") or []
        if terms:
            lines.append("Coordinates: " + _names(pjc.get("observable_coordinates") or [], full=section != "all"))
            lines.append(f"Event classes: {len(terms)}")
            if section in {"all", "jumps"}:
                shown = terms if section == "jumps" else (terms[:PREVIEW_EVENTS] if len(terms) <= 20 else [])
                rows = []
                for term in shown:
                    activity = term["propensity_polynomial"]
                    activity_text = activity["text"]
                    if len(activity.get("terms", [])) > 1:
                        activity_text = f"({activity_text})"
                    rows.append([term["event_class"], _jump(term["delta"]),
                                 f"{rate_symbols.get(term['rate'], term['rate'])} × {activity_text}"])
                if rows:
                    lines.append(_event_table(rows))
                if len(shown) < len(terms):
                    lines.append(f"  … {len(terms) - len(shown)} event rows omitted; use --section jumps for all.")
            moment_images = pjc.get("first_moment_generator_images") or {}
            if section == "all" and moment_images and linear.get("status") != "CLOSED":
                lines.append("Exact first-moment identities from these events:")
                for coordinate, channels in list(moment_images.items())[:PREVIEW_EQUATIONS]:
                    lines.append(f"  d⟨{coordinate}⟩/dt = {_rhs(channels, rate_symbols)}")
            if section in {"all", "operator"}:
                lines.append("D = Σ_e rate_e (exp(λ·Δ_e) − 1) a_e(∂λ)")
                lines.append("Each event row supplies Δ_e and a_e; replace x by ∂λ in its activity.")
                lines.append(f"Maximum differential order: {operator.get('maximum_differential_order', 'not recorded')}")
                if section == "operator":
                    lines.append("Explicit terms of D (one per derived event class):")
                    for term in terms:
                        rate = rate_symbols.get(term["rate"], term["rate"])
                        factor = term["factor"].replace("lambda_", "λ_")
                        lines.append(
                            f"  {term['event_class']}: {rate} [{factor}] "
                            f"({_differential(term['propensity_polynomial'])})"
                        )
        elif pjc.get("status") == "CLOSED" and operator.get("terms") == []:
            lines.extend(["Event classes: 0", "D = 0; all selected coordinates are invariant."])
        else:
            lines.append("No finite event table or EMGF operator is established in this record.")
    if section == "all":
        _preview(lines, result)
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("record", type=Path, help="JSON from mfobs.example_cli or run.py rop-analysis")
    parser.add_argument("--section", choices=("all", "ode", "jumps", "operator"), default="all")
    args = parser.parse_args()
    print(render(json.loads(args.record.read_text(encoding="utf-8")), section=args.section))


if __name__ == "__main__":
    main()
