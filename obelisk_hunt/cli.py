"""CLI: scan an assembly's contigs for circularity evidence.

    python -m obelisk_hunt contigs.fasta [--gfa assembly_graph.gfa] [--out circular_contigs.tsv]

Writes one row per contig with the sequence-overlap evidence, the
assembly-graph self-loop evidence (if a GFA was given), and a combined
verdict. This is Week-1 scope only: circularity detection. It does not filter
by database homology and does not score structure stability -- those are
separate, later stages of the pipeline.
"""

from __future__ import annotations

import argparse
import csv
import sys

from Bio import SeqIO

from .circularity import detect_self_overlap
from .gfa import parse_gfa


def classify(seq_overlap, graph_self_loop: bool) -> str:
    if seq_overlap and graph_self_loop:
        return "circular_graph+sequence"
    if graph_self_loop:
        return "circular_graph_only"
    if seq_overlap:
        return "circular_sequence_only"
    return "linear"


def run(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fasta", help="Assembly contigs FASTA")
    parser.add_argument("--gfa", default=None, help="Assembly graph GFA1 (optional)")
    parser.add_argument("--out", default=None, help="Output TSV path (default: stdout)")
    parser.add_argument("--min-length", type=int, default=100,
                         help="Skip contigs shorter than this (default: 100)")
    parser.add_argument("--min-overlap", type=int, default=20,
                         help="Minimum terminal-overlap length to call circular (default: 20)")
    parser.add_argument("--max-edit-frac", type=float, default=0.10,
                         help="Max fraction of edits allowed in the terminal overlap (default: 0.10)")
    args = parser.parse_args(argv)

    graph = parse_gfa(args.gfa) if args.gfa else None

    out_fh = open(args.out, "w", newline="") if args.out else sys.stdout
    writer = csv.writer(out_fh, delimiter="\t")
    writer.writerow([
        "contig", "length", "verdict", "overlap_len", "identity",
        "trimmed_length", "graph_self_loop",
    ])

    n_circular = 0
    n_total = 0
    for record in SeqIO.parse(args.fasta, "fasta"):
        length = len(record.seq)
        if length < args.min_length:
            continue
        n_total += 1

        overlap = detect_self_overlap(
            str(record.seq),
            min_overlap=args.min_overlap,
            max_edit_frac=args.max_edit_frac,
        )
        graph_loop = bool(graph and graph.has_self_loop(record.id))
        verdict = classify(overlap, graph_loop)
        if verdict != "linear":
            n_circular += 1

        writer.writerow([
            record.id,
            length,
            verdict,
            overlap.overlap_len if overlap else "",
            f"{overlap.identity:.4f}" if overlap else "",
            overlap.trimmed_length if overlap else "",
            graph_loop,
        ])

    if args.out:
        out_fh.close()

    print(f"[obelisk_hunt] {n_circular}/{n_total} contigs >= {args.min_length}bp "
          f"flagged circular", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(run())
