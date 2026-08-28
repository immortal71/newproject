"""CLI: filter candidates by database homology (the "zero database homology"
half of the signature).

    python -m obelisk_hunt filter-homology candidates.fasta --db reference.fasta \
        --no-hit-fasta survivors.fasta --hits-tsv hits.tsv

Sequences with no hit clearing --min-identity and --min-coverage against
--db pass through to --no-hit-fasta. --hits-tsv records every sequence's
best hit (or lack of one) for audit.
"""

from __future__ import annotations

import argparse
import csv
import sys

from Bio import SeqIO

from .homology import Minimap2NotFound, check_homology


def run(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fasta", help="Candidate sequences")
    parser.add_argument("--db", required=True, help="Reference FASTA to search against")
    parser.add_argument("--no-hit-fasta", default=None,
                         help="Write sequences with zero database homology here")
    parser.add_argument("--hits-tsv", default=None, help="Per-sequence hit report (default: stdout)")
    parser.add_argument("--min-identity", type=float, default=0.80,
                         help="Minimum identity to count as a hit (default: 0.80)")
    parser.add_argument("--min-coverage", type=float, default=0.50,
                         help="Minimum query coverage to count as a hit (default: 0.50)")
    parser.add_argument("--preset", default=None, help="minimap2 -x preset (default: minimap2's own default)")
    parser.add_argument("--minimap2-path", default=None, help="Path to minimap2 binary (default: search PATH)")
    args = parser.parse_args(argv)

    try:
        results = check_homology(
            args.fasta, args.db,
            min_identity=args.min_identity,
            min_coverage=args.min_coverage,
            preset=args.preset,
            minimap2_path=args.minimap2_path,
        )
    except Minimap2NotFound as e:
        print(f"[obelisk_hunt] {e}", file=sys.stderr)
        return 1

    hits_fh = open(args.hits_tsv, "w", newline="") if args.hits_tsv else sys.stdout
    writer = csv.writer(hits_fh, delimiter="\t")
    writer.writerow(["contig", "has_hit", "best_identity", "best_coverage", "best_target"])

    no_hit_fh = open(args.no_hit_fasta, "w") if args.no_hit_fasta else None
    n_no_hit = 0

    for record in SeqIO.parse(args.fasta, "fasta"):
        result = results[record.id]
        writer.writerow([
            record.id, result.has_hit,
            f"{result.best_identity:.4f}", f"{result.best_coverage:.4f}",
            result.best_target or "",
        ])
        if not result.has_hit:
            n_no_hit += 1
            if no_hit_fh:
                no_hit_fh.write(f">{record.id}\n{record.seq}\n")

    if args.hits_tsv:
        hits_fh.close()
    if no_hit_fh:
        no_hit_fh.close()

    print(f"[obelisk_hunt] {n_no_hit} sequences with zero database homology "
          f"(identity>={args.min_identity}, coverage>={args.min_coverage})", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(run())
