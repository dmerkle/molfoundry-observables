"""Independent finite eight-slot semantics; NEVER imported by synthesis."""
from itertools import product
from mfobs.checks import local_slots_reference as reference


TOPOLOGIES = {'8': (8,), '4+4': (4, 4), '3+5': (3, 5)}


def topology_edges(lengths):
    edges = []; offset = 0
    for n in lengths:
        edges.extend((a+offset, b+offset) for a, b in reference.cycle(n))
        offset += n
    return tuple(edges)


def assignments(lengths):
    return product((0, 1), repeat=sum(lengths))
