"""The third clause of the signature: zero database homology.

A candidate that's circular and structurally stable but turns out to be a
known rRNA fragment, a known plasmid, or a known viroid isn't dark matter --
it's a rediscovery. This module checks each candidate against a reference
FASTA/database with minimap2 and reports whether it has a significant hit.

minimap2 (not BLASTN/mmseqs2) is the backend here because it's a single
static binary installable straight from its GitHub releases with no
dependency chain -- the same install path already used in this project for
rnaSPAdes/MEGAHIT. For a real deployment against nt/RefSeq at scale, a
proper BLASTN/mmseqs2 search is more standard and more sensitive to short,
divergent hits; swapping the backend later means changing only
`run_minimap2`, not anything that calls `check_homology`.
"""

from __future__ import annotations

import shutil
import subprocess
from dataclasses import dataclass


class Minimap2NotFound(RuntimeError):
    pass


@dataclass
class Hit:
    query: str
    query_len: int
    query_start: int
    query_end: int
    target: str
    target_len: int
    matches: int
    aln_len: int

    @property
    def identity(self) -> float:
        return self.matches / self.aln_len if self.aln_len else 0.0

    @property
    def query_coverage(self) -> float:
        return (self.query_end - self.query_start) / self.query_len if self.query_len else 0.0


@dataclass
class HomologyResult:
    has_hit: bool
    best_identity: float
    best_coverage: float
    best_target: str | None


def find_minimap2(minimap2_path: str | None = None) -> str:
    path = minimap2_path or shutil.which("minimap2")
    if not path:
        raise Minimap2NotFound(
            "minimap2 not found. Install it (bioconda `minimap2`, or the static "
            "binary from https://github.com/lh3/minimap2/releases) and either "
            "put it on PATH or pass --minimap2-path."
        )
    return path


def run_minimap2(
    query_fasta: str,
    reference_fasta: str,
    preset: str | None = None,
    minimap2_path: str | None = None,
) -> list[Hit]:
    """Align every sequence in `query_fasta` against `reference_fasta` and
    return one Hit per PAF alignment record (a query can have zero, one, or
    several)."""
    minimap2 = find_minimap2(minimap2_path)

    cmd = [minimap2, "-c"]
    if preset:
        cmd += ["-x", preset]
    cmd += [reference_fasta, query_fasta]

    result = subprocess.run(cmd, capture_output=True, text=True, check=True)

    hits = []
    for line in result.stdout.splitlines():
        fields = line.split("\t")
        hits.append(Hit(
            query=fields[0],
            query_len=int(fields[1]),
            query_start=int(fields[2]),
            query_end=int(fields[3]),
            target=fields[5],
            target_len=int(fields[6]),
            matches=int(fields[9]),
            aln_len=int(fields[10]),
        ))
    return hits


def check_homology(
    candidates_fasta: str,
    reference_fasta: str,
    min_identity: float = 0.80,
    min_coverage: float = 0.50,
    preset: str | None = None,
    minimap2_path: str | None = None,
) -> dict[str, HomologyResult]:
    """Return one HomologyResult per sequence in `candidates_fasta`. A
    sequence absent from the dict's hits or below both thresholds on every
    hit has "zero database homology" against this reference."""
    from Bio import SeqIO

    query_ids = [r.id for r in SeqIO.parse(candidates_fasta, "fasta")]
    hits = run_minimap2(candidates_fasta, reference_fasta, preset=preset, minimap2_path=minimap2_path)

    best: dict[str, Hit] = {}
    for hit in hits:
        current = best.get(hit.query)
        if current is None or hit.identity * hit.query_coverage > current.identity * current.query_coverage:
            best[hit.query] = hit

    results = {}
    for qid in query_ids:
        hit = best.get(qid)
        if hit is None:
            results[qid] = HomologyResult(has_hit=False, best_identity=0.0, best_coverage=0.0, best_target=None)
            continue
        has_hit = hit.identity >= min_identity and hit.query_coverage >= min_coverage
        results[qid] = HomologyResult(
            has_hit=has_hit,
            best_identity=hit.identity,
            best_coverage=hit.query_coverage,
            best_target=hit.target,
        )
    return results
