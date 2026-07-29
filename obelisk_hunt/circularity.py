"""Sequence-level circular-contig detection.

Linear assemblers (rnaSPAdes, MEGAHIT, ...) report a circular molecule as a
linear contig broken at an arbitrary point. The signature this leaves behind
is a terminal redundancy: the end of the contig re-reads the beginning,
because the assembler's path walked all the way around the cycle and back
onto its own start k-mer. This module looks for that redundancy directly in
the contig sequence, independent of any assembly graph.

Reference for the heuristic: the same end-vs-start overlap trim used by
Recycler/Unicycler for circular plasmid contigs, and by circlator-style
finishing tools for circularizing genome assemblies.

Detection is two-pass:
  1. Scan candidate overlap lengths k (from longest to shortest) and score
     seq[:k] against seq[-k:] by Hamming distance, vectorized with numpy.
     This is the "does a border of length k exist" question -- equivalent to
     what the KMP failure function gives exactly for k=0 mismatches, extended
     here to tolerate a few substitutions. We scan from the longest candidate
     down because a k a little larger than the true overlap collapses to
     near-random identity within a few extra bases (4-letter alphabet), so
     the first k that clears the identity threshold is effectively the true
     overlap length, not a coincidental shorter sub-border.
  2. Refine that single best candidate with an edlib edit-distance alignment,
     which (unlike Hamming) tolerates indels in the join -- assembly errors
     at a circular join are not always clean substitutions.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import edlib
import numpy as np


@dataclass
class SelfOverlap:
    overlap_len: int
    edit_distance: int
    identity: float
    trimmed_length: int


def detect_self_overlap(
    seq: str,
    min_overlap: int = 20,
    max_overlap_frac: float = 0.5,
    max_overlap_abs: int = 300,
    max_edit_frac: float = 0.10,
) -> Optional[SelfOverlap]:
    """Test whether `seq`'s end redundantly re-reads its own start.

    Returns None if the sequence is too short to test, or no qualifying
    overlap is found.
    """
    seq = seq.upper()
    length = len(seq)

    max_overlap = int(min(max_overlap_abs, max_overlap_frac * length, length // 2))
    if max_overlap < min_overlap:
        return None

    arr = np.frombuffer(seq.encode("ascii"), dtype=np.uint8)

    candidate_k = None
    for k in range(max_overlap, min_overlap - 1, -1):
        mismatches = int(np.count_nonzero(arr[:k] != arr[-k:]))
        identity = 1.0 - mismatches / k
        if identity >= 1.0 - max_edit_frac:
            candidate_k = k
            break

    if candidate_k is None:
        return None

    # Refine with an indel-tolerant alignment on just the candidate window.
    head = seq[:candidate_k]
    tail = seq[-candidate_k:]
    aligned = edlib.align(head, tail, mode="NW", task="distance")
    edit_distance = aligned["editDistance"]
    identity = 1.0 - edit_distance / candidate_k
    if identity < 1.0 - max_edit_frac:
        return None

    return SelfOverlap(
        overlap_len=candidate_k,
        edit_distance=edit_distance,
        identity=identity,
        trimmed_length=length - candidate_k,
    )
