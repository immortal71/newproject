"""CLI: score each sequence in a FASTA against its own dinucleotide-shuffled
null (the structure-stability half of the project's signature).

    python -m obelisk_hunt score-stability trimmed_circular.fasta --n-shuffles 100 --out stability.tsv

Typical input is the --trimmed-fasta output of `detect-circular`. This step
does not filter by database homology -- that's still a separate stage.
"""

from __future__ import annotations

import argparse
import csv
import sys

from Bio import SeqIO

from .structure import score_stability


def run(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fasta", help="Candidate sequences (e.g. detect-circular's --trimmed-fasta)")
    parser.add_argument("--out", default=None, help="Output TSV path (default: stdout)")
    parser.add_argument("--n-shuffles", type=int, default=100,
                         help="Dinucleotide-shuffled controls per sequence (default: 100)")
    parser.add_argument("--seed", type=int, default=None, help="RNG seed for reproducibility")
    args = parser.parse_args(argv)

    out_fh = open(args.out, "w", newline="") if args.out else sys.stdout
    writer = csv.writer(out_fh, delimiter="\t")
    writer.writerow([
        "contig", "length", "gc_content", "real_mfe",
        "shuffled_mean_mfe", "shuffled_std_mfe", "z_score", "empirical_p", "n_shuffles",
    ])

    n_significant = 0
    n_total = 0
    for record in SeqIO.parse(args.fasta, "fasta"):
        n_total += 1
        result = score_stability(str(record.seq), n_shuffles=args.n_shuffles, seed=args.seed)
        if result.empirical_p <= 0.05:
            n_significant += 1

        writer.writerow([
            record.id,
            result.length,
            f"{result.gc_content:.4f}",
            f"{result.real_mfe:.2f}",
            f"{result.shuffled_mean_mfe:.2f}",
            f"{result.shuffled_std_mfe:.2f}",
            f"{result.z_score:.3f}",
            f"{result.empirical_p:.4f}",
            result.n_shuffles,
        ])

    if args.out:
        out_fh.close()

    print(f"[obelisk_hunt] {n_significant}/{n_total} sequences more stable than shuffled "
          f"controls at p<=0.05", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(run())
