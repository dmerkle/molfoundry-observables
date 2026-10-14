# Manifest and recursive-language contracts

This page describes manifest validation and the recursive-chain language
certificate. The general backend interface is documented in the
[Python API](public-api.md) and [algorithm map](../notebooks/visual/CORE.md).

`CLOSED` requires exact composition and reduction on the declared graph family.
The generic and invariant adapters reject connected-only overlap enumeration
and a negative theory incompatible with that family.

## Full overlap enumeration

Exact conditional-DPO commutators require every admissible overlap, including
disconnected overlaps. The generic and recursive-language adapters reject
`connected=True` before composition. This is an invalid-input error, not a
proof of nonclosure. A future optimization may use a restricted overlap set
only with a separate completeness argument.

The **digest-bound Ring8 adapter** also supports connected-only composition. Its final
result compares every rate-separated derived row with an independent generator
on every dihedral orbit of its entire finite state space. The digest prevents
silently applying that contract to changed rules, language, or overlap policy;
failed comparison raises an error and cannot emit `CLOSED`. This exception is
not exported as a generic user-provided `trust_connected` switch.

## Numerical and signature validation

Forbidden theories declare a finite atom/edge signature, nonnegative integer
bond orders, nonnegative integer valence maxima, and nonnegative integer
incidence maxima. Booleans, negative numbers, floats, strings, and integers
outside JavaScript's exact integer range are rejected. Zero-order relations are
allowed; they contribute zero covalent valence and can have incidence bounds.
All allowed pairs and incidence bounds must refer to declared labels; edge
orders and pair tables must have identical keys. Atom labels are unique and
pair endpoints are canonicalized as unordered pairs. `|` is excluded from
labels because it is the backend's pair-key delimiter.

The manifest loader checks this schema. The dataclasses also validate direct
Python construction, and payload/certificate creation revalidates mutable
mapping values. Runtime adapter construction rechecks the overlap setting, so
editing an already-loaded input in memory does not bypass the guard.

## All-size language inclusion in the negative theory

Let the declared states be disjoint unions

    C_n + S_x + M R,        n,M >= 0, x in {0,1},

where `S_0,S_1,R` are fixed finite graphs and

    C_0 = Z,               C_n = L B^(n-1) T for n>=1.

The port grammar joins only consecutive tiles, using one outgoing and one
incoming vertex; it requires uniform incoming and outgoing port labels. A
port is a tile vertex, not a synthetic atom label. Both ports of a repeat tile
may be the **same** vertex. Predecessor and successor vertices remain distinct.

The supported forbidden family consists of labelled vertices, labelled edges,
overvalent stars and incidence-bound stars. Whether such a pattern occurs is
determined by a vertex's label and its **complete incident neighbourhood**.
Consequently its membership on the whole infinite language is decidable by
this finite context exhaustion:

1. Check every vertex star in `Z`, both fixed singleton states, and `R`.
2. For a vertex in `L`, add its outgoing port neighbour from `T` or `B`.
3. For a vertex in `T`, add its incoming port neighbour from `L` or `B`.
4. For a vertex in `B`, add both neighbours for all four contexts
   `L-B-T`, `L-B-B`, `B-B-T`, and `B-B-B`.

These eight tile contexts exhaust the adjacent tile types permitted by the
grammar. Internal edges are retained, and external port edges are added
**simultaneously**, including the two contributions when the ports coincide.
Every star of every `C_n` is one of these contexts; conversely every listed
context is realizable in a finite grammar word. Their validity therefore
implies, and is implied by, validity for all lengths. Disjoint union introduces
no new incident edges. Every forbidden pattern is connected, so arbitrary
repeated-component multiplicities introduce no new violation either.

The implementation constructs the local stars directly. It does not enumerate
chain lengths or assume that a chosen finite cutoff is representative. The
result includes `language_forbidden_contract`, with input digests, all contexts,
their checked stars, and `VERIFIED_ALL_SIZES`. A test disables chain instantiation
entirely while producing this certificate. A separate long-chain test is merely
an independent regression check, not part of the proof.

## What this does not prove

Language inclusion in the local negative theory is **not** rule preservation
of the language. Preservation remains a separate supplied model argument:
every enabled rule application from a declared state must again yield a
declared state, with the documented component-coordinate update. No finite
host test alone proves this unbounded claim. The result metadata explicitly
retains this obligation. Arbitrary nested conditions, arbitrary graph grammars,
termination/nonclosure, chemical realism and experimental rate calibration
are also not certified by this checker.
