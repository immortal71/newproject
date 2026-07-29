"""The core test: is this sequence's secondary structure more stable than its
own dinucleotide-shuffled controls?

A random sequence folds into something mediocre. A sequence under selection
to hold a shape folds far better than its own shuffles, which share its
length, base composition, and dinucleotide composition but not its base-pair
register. This module folds the candidate once with ViennaRNA, folds many
dinucleotide-shuffled versions of it, and reports where the real sequence
sits in that null distribution.

Sign convention: MFE is kcal/mol, more negative = more stable. We define
z_score so that a *positive* z means "the real sequence is more stable than
its shuffled controls" (the signature we're looking for):

    z = (mean(shuffled_mfe) - real_mfe) / std(shuffled_mfe)

empirical_p is the one-sided permutation p-value for "the real sequence is at
least this stable by chance": the fraction of shuffles whose MFE is at least
as negative as the real sequence's, with the standard +1/+1 correction.
"""

from __future__ import annotations

import random
import statistics
from dataclasses import dataclass

import RNA

from .shuffle import dinucleotide_shuffle


@dataclass
class StabilityResult:
    length: int
    gc_content: float
    real_mfe: float
    shuffled_mean_mfe: float
    shuffled_std_mfe: float
    z_score: float
    empirical_p: float
    n_shuffles: int


def fold_mfe(seq: str) -> float:
    rna_seq = seq.upper().replace("T", "U")
    _structure, mfe = RNA.fold(rna_seq)
    return mfe


def gc_content(seq: str) -> float:
    seq = seq.upper()
    gc = seq.count("G") + seq.count("C")
    return gc / len(seq) if seq else 0.0


def score_stability(seq: str, n_shuffles: int = 100, seed: int | None = None) -> StabilityResult:
    rng = random.Random(seed)
    real_mfe = fold_mfe(seq)

    shuffled_mfes = [fold_mfe(dinucleotide_shuffle(seq, rng)) for _ in range(n_shuffles)]

    mean_shuffled = statistics.mean(shuffled_mfes)
    std_shuffled = statistics.pstdev(shuffled_mfes)

    z_score = (mean_shuffled - real_mfe) / std_shuffled if std_shuffled > 0 else float("inf")
    n_at_least_as_stable = sum(1 for m in shuffled_mfes if m <= real_mfe)
    empirical_p = (n_at_least_as_stable + 1) / (n_shuffles + 1)

    return StabilityResult(
        length=len(seq),
        gc_content=gc_content(seq),
        real_mfe=real_mfe,
        shuffled_mean_mfe=mean_shuffled,
        shuffled_std_mfe=std_shuffled,
        z_score=z_score,
        empirical_p=empirical_p,
        n_shuffles=n_shuffles,
    )
