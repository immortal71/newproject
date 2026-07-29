#!/usr/bin/env bash
# Assemble paired-end metatranscriptome reads with rnaSPAdes, producing both
# a contigs FASTA and the assembly graph (needed for the circularity
# detector's graph-based evidence).
#
# Usage:
#   scripts/assemble.sh <reads_1.fastq.gz> <reads_2.fastq.gz> <output_dir> [threads]
set -euo pipefail

if [[ $# -lt 3 ]]; then
  echo "Usage: $0 <reads_1.fastq.gz> <reads_2.fastq.gz> <output_dir> [threads]" >&2
  exit 1
fi

R1="$1"
R2="$2"
OUTDIR="$3"
THREADS="${4:-4}"

if ! command -v rnaspades.py >/dev/null; then
  echo "[assemble] ERROR: rnaspades.py not found. Install via environment.yml" \
       "(conda env create -f environment.yml)" >&2
  exit 1
fi

mkdir -p "$OUTDIR"
rnaspades.py -1 "$R1" -2 "$R2" -o "$OUTDIR" -t "$THREADS"

echo "[assemble] contigs:        ${OUTDIR}/transcripts.fasta"
echo "[assemble] assembly graph: ${OUTDIR}/assembly_graph_with_scaffolds.gfa"
echo "[assemble] next: python -m obelisk_hunt ${OUTDIR}/transcripts.fasta --gfa ${OUTDIR}/assembly_graph_with_scaffolds.gfa --out ${OUTDIR}/circular_contigs.tsv"
