#!/usr/bin/env bash
# Fetch one iHMP (or any SRA/ENA) metatranscriptome run as paired FASTQ.
#
# NOTE: this needs real outbound access to NCBI SRA and/or EBI ENA. It will
# NOT run inside a network-sandboxed environment that blocks those hosts --
# run it on a workstation, HPC login node, or a cloud VM with normal internet
# access.
#
# Usage:
#   scripts/fetch_sample.sh <SRA_or_ERR_accession> <output_dir>
#
# Tries, in order:
#   1. SRA toolkit (prefetch + fasterq-dump) -- works for any SRA accession.
#   2. Direct ENA FTP fastq.gz download -- works when the run has an ENA
#      mirror, and needs no SRA toolkit at all.
set -euo pipefail

if [[ $# -ne 2 ]]; then
  echo "Usage: $0 <accession> <output_dir>" >&2
  exit 1
fi

ACCESSION="$1"
OUTDIR="$2"
mkdir -p "$OUTDIR"

fetch_via_sra_toolkit() {
  command -v prefetch >/dev/null && command -v fasterq-dump >/dev/null
}

fetch_via_ena() {
  command -v curl >/dev/null
}

if fetch_via_sra_toolkit; then
  echo "[fetch_sample] using SRA toolkit for ${ACCESSION}"
  prefetch --output-directory "$OUTDIR" "$ACCESSION"
  fasterq-dump --split-files --outdir "$OUTDIR" "$OUTDIR/${ACCESSION}/${ACCESSION}.sra"
  gzip -f "$OUTDIR/${ACCESSION}"_*.fastq
  echo "[fetch_sample] wrote ${OUTDIR}/${ACCESSION}_1.fastq.gz (+_2 if paired)"
  exit 0
fi

if fetch_via_ena; then
  echo "[fetch_sample] SRA toolkit not found, falling back to ENA FTP for ${ACCESSION}"
  REPORT_URL="https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${ACCESSION}&result=read_run&fields=fastq_ftp&format=tsv"
  FTP_PATHS=$(curl -sS "$REPORT_URL" | tail -n +2 | cut -f1)
  if [[ -z "$FTP_PATHS" ]]; then
    echo "[fetch_sample] ERROR: ENA has no fastq_ftp entry for ${ACCESSION}" >&2
    exit 1
  fi
  IFS=';' read -ra URLS <<< "$FTP_PATHS"
  for url in "${URLS[@]}"; do
    echo "[fetch_sample] downloading ftp://${url}"
    curl -sS -o "${OUTDIR}/$(basename "$url")" "ftp://${url}"
  done
  exit 0
fi

echo "[fetch_sample] ERROR: neither SRA toolkit nor curl is available" >&2
exit 1
