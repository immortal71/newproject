"""Dinucleotide-preserving sequence shuffling -- the matched null model.

A mononucleotide (per-base) shuffle only matches composition; it does not
preserve dinucleotide frequencies, so it can systematically over- or
under-estimate how foldable a random sequence of the same base composition
"should" be (stacking energies are a dinucleotide-level property). The
standard fix, used throughout this kind of analysis (and named explicitly in
the project brief), is an Altschul-Erikson-style dinucleotide shuffle: model
the sequence as a walk over a multigraph whose nodes are the four bases and
whose edges are the sequence's own dinucleotides, then generate a *different*
Eulerian path through that same multigraph (same start, same end, every edge
used exactly once). Any such path is, by construction, a sequence with
identical length, identical base composition, and identical dinucleotide
counts to the original -- it only reshuffles which occurrences of each
dinucleotide follow which.
"""

from __future__ import annotations

import random


def dinucleotide_shuffle(seq: str, rng: random.Random) -> str:
    """Return a random sequence with the same length, base composition, and
    dinucleotide composition as `seq`."""
    seq = seq.upper()
    n = len(seq)
    if n < 3:
        return seq

    first_char = seq[0]
    last_char = seq[-1]

    edges: dict[str, list[str]] = {}
    for i in range(n - 1):
        a, b = seq[i], seq[i + 1]
        edges.setdefault(a, []).append(b)

    # Shuffle each node's out-edge order. A stack-based Eulerian-path walk
    # (Hierholzer's algorithm) is correct regardless of the order edges are
    # tried at each node, *provided* an Eulerian path from first_char to
    # last_char exists in the graph -- and it must, since the original
    # sequence is itself a witness that one does.
    remaining = {node: out.copy() for node, out in edges.items()}
    for out in remaining.values():
        rng.shuffle(out)

    stack = [first_char]
    path: list[str] = []
    while stack:
        node = stack[-1]
        out = remaining.get(node)
        if out:
            stack.append(out.pop())
        else:
            path.append(stack.pop())
    path.reverse()

    shuffled = "".join(path)
    assert shuffled[0] == first_char and shuffled[-1] == last_char
    assert len(shuffled) == n
    return shuffled
