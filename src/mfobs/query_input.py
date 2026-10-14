"""Common, file-based observable-query input for every analysis backend.

One observable may be a finite rational linear combination of identity rules.
Repeated ``NAME=FILE`` arguments therefore denote terms of the same observable,
not several accidentally colliding queries.  Model backends may support a
strict subset of these expressions, but they all receive the same object and
must reject unsupported shapes explicitly.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import re
from typing import Iterable


_NAME = re.compile(r"[A-Za-z_][A-Za-z_0-9]*\Z")


@dataclass(frozen=True)
class QueryTerm:
    file: Path
    coefficient: Fraction = Fraction(1)

    @property
    def gml(self) -> str:
        return self.file.read_text(encoding="utf-8")

    def record(self) -> dict:
        text = self.gml
        return {
            "file": str(self.file),
            "coefficient": str(self.coefficient),
            "sha256": sha256(text.encode()).hexdigest(),
        }


@dataclass(frozen=True)
class ObservableQuery:
    name: str
    terms: tuple[QueryTerm, ...]

    def __post_init__(self) -> None:
        if not _NAME.fullmatch(self.name):
            raise ValueError(f"invalid observable name {self.name!r}")
        if not self.terms or any(not term.coefficient for term in self.terms):
            raise ValueError(f"observable {self.name!r} needs nonzero terms")

    def record(self) -> dict:
        return {"name": self.name, "terms": [term.record() for term in self.terms]}


@dataclass(frozen=True)
class QuerySelection:
    observables: tuple[ObservableQuery, ...]

    def __post_init__(self) -> None:
        names = [query.name for query in self.observables]
        if not names or len(names) != len(set(names)):
            raise ValueError("query observable names must be nonempty and unique")
        for query in self.observables:
            for term in query.terms:
                if not term.file.is_file():
                    raise ValueError(f"observable file does not exist: {term.file}")

    @property
    def names(self) -> tuple[str, ...]:
        return tuple(query.name for query in self.observables)

    def record(self) -> dict:
        payload = {"schema_version": 1, "observables": [q.record() for q in self.observables]}
        digest_payload = [{"name": q.name, "terms": [
            {"coefficient": str(t.coefficient), "sha256": sha256(t.gml.encode()).hexdigest()}
            for t in q.terms
        ]} for q in self.observables]
        payload["definition_sha256"] = sha256(json.dumps(
            digest_payload, sort_keys=True, separators=(",", ":")
        ).encode()).hexdigest()
        return payload


def _term(value: str, *, base: Path) -> tuple[str, QueryTerm]:
    try:
        name, raw_file = value.split("=", 1)
    except ValueError as exc:
        raise ValueError("--observable requires NAME=FILE") from exc
    coefficient = Fraction(1)
    if "@" in raw_file:
        raw_file, raw_coefficient = raw_file.rsplit("@", 1)
        coefficient = Fraction(raw_coefficient)
    path = Path(raw_file).expanduser()
    if not path.is_absolute():
        path = base / path
    return name, QueryTerm(path.resolve(), coefficient)


def from_arguments(values: Iterable[str], *, base: Path | None = None) -> QuerySelection:
    """Parse repeatable ``NAME=FILE[@COEFFICIENT]`` arguments."""

    base = Path.cwd() if base is None else Path(base)
    grouped: dict[str, list[QueryTerm]] = {}
    for value in values:
        name, term = _term(value, base=base)
        grouped.setdefault(name, []).append(term)
    return QuerySelection(tuple(
        ObservableQuery(name, tuple(terms)) for name, terms in grouped.items()
    ))


def from_file(path: str | Path) -> QuerySelection:
    """Load the same representation from a small explicit JSON file."""

    source = Path(path).expanduser().resolve()
    data = json.loads(source.read_text(encoding="utf-8"))
    if set(data) != {"schema_version", "observables"} or data["schema_version"] != 1:
        raise ValueError("query file must contain schema_version=1 and observables")
    values = []
    for observable in data["observables"]:
        if set(observable) != {"name", "terms"}:
            raise ValueError("each query observable needs exactly name and terms")
        for term in observable["terms"]:
            if set(term) - {"file", "coefficient"} or "file" not in term:
                raise ValueError("each query term needs file and optional coefficient")
            value = f"{observable['name']}={term['file']}"
            if "coefficient" in term:
                value += f"@{term['coefficient']}"
            values.append(value)
    return from_arguments(values, base=source.parent)
