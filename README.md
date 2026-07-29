# RNA dark matter: a structure-first search

Obelisks were found by defining a signature, not by matching sequences. This
project sweeps for a different signature:

> Circular RNA contigs with zero database homology whose predicted secondary
> structure is significantly more stable than dinucleotide-shuffled controls
> matched on length and GC.

Anything that clears that bar is either a new class of replicon, a new
ribozyme family, or a characterizable artifact.

## Plan

- **Weeks 1-3 -- reproduce.** Run obelisk detection on iHMP gut data. Recover
  known Obelisks as a positive control.
- **Weeks 4-6 -- calibrate.** Inject synthetic structured circular elements
  at known abundance into real data to establish detection limits.
- **Weeks 7-12 -- sweep.** Soil, hot springs, insect/plant hosts, anaerobic
  digesters, glacier meltwater -- datasets in SRA that haven't been searched
  for this signature.
- **Weeks 13+ -- characterize.** Fold structures, look for ORFs, check
  cross-sample co-occurrence, build phylogenies of any protein families that
  turn up.

## Status: Week 1, step 1

Per plan: *"pull one iHMP metatranscriptome, assemble it, and write the
circularity detector. Nothing else."*

**The circularity detector is built and tested.** `obelisk_hunt/` implements
it; `tests/test_circularity.py` validates it against synthetic positive and
negative controls (see below).

**The fetch-and-assemble-from-scratch half did not run for real iHMP data in
this session.** This repo was built inside a network-sandboxed environment
whose proxy rejects connections to NCBI (`www.ncbi.nlm.nih.gov`,
`sra-download.ncbi.nlm.nih.gov`, `trace.ncbi.nlm.nih.gov`, `pmc.ncbi.nlm.nih.gov`)
and EBI (`www.ebi.ac.uk`, `ftp.ebi.ac.uk`) -- verified with direct `curl` and
`WebFetch` probes, all return `403` at the proxy, not at the origin. Those
are the two places iHMP metatranscriptome reads actually live, and I'm not
reporting a fabricated "ran it, got N contigs, M circular" result for data I
couldn't touch.

`scripts/fetch_sample.sh` and `scripts/assemble.sh` are real, runnable
pipeline code, not placeholders. Two things worth knowing before you run them
for real: (1) `rnaspades.py` and `megahit` install cleanly from GitHub release
binaries even inside this sandbox -- `objects.githubusercontent.com` and
`raw.githubusercontent.com` are not blocked, only `github.com`'s HTML/API
surface is; (2) `sra-tools` is *not* a single static binary -- its bioconda
package pulls in `ncbi-vdb`, `ossuuid`, `perl-xml-libxml`, etc., which is real
dependency-resolution work, better done with `conda`/`mamba` proper (as
`environment.yml` does) than hand-extracted package-by-package. `conda.anaconda.org`
(bioconda's actual package host) was reachable from here, for what that's
worth, if a future session wants to push on this further.

### What "reproduce known Obelisks" actually requires

To claim the positive control, you need real iHMP reads run through the
real pipeline below, on a real Obelisk-positive sample. A specific SRA
accession isn't hard-pinned here on purpose: pick one from the iHMP/IBDMDB
metatranscriptome cohort (or from the accession list in the Obelisks paper's
supplement) after confirming it in the SRA/ENA browser -- don't take an
accession number from memory, verify it against the source. The Obelisks
paper's own detection tool, [VNom](https://github.com/Zheludev/VNom), is
public and MIT-licensed; its README documents the same core heuristic this
detector uses ("identify contigs with terminal k-mer repeats, consistent
with circularity") and is worth reading before the real run.

## Running it for real (outside this sandbox)

```bash
conda env create -f environment.yml
conda activate obelisk-hunt

# 1. fetch
scripts/fetch_sample.sh <SRA_accession> data/<accession>

# 2. assemble (rnaSPAdes; produces both contigs and the assembly graph)
scripts/assemble.sh data/<accession>/<accession>_1.fastq.gz \
                    data/<accession>/<accession>_2.fastq.gz \
                    data/<accession>_assembly

# 3. detect circular contigs
python -m obelisk_hunt data/<accession>_assembly/transcripts.fasta \
    --gfa data/<accession>_assembly/assembly_graph_with_scaffolds.gfa \
    --out data/<accession>_assembly/circular_contigs.tsv
```

## The circularity detector

`obelisk_hunt/circularity.py` looks for the terminal-redundancy signature a
linear assembler leaves on a circular molecule: the contig's end re-reads its
own start, because the assembly path walked all the way around the loop.
It scans candidate overlap lengths (longest first, vectorized Hamming
distance) to find the border, then refines with an edlib edit-distance
alignment so a few sequencing/assembly errors in the join don't break
detection.

`obelisk_hunt/gfa.py` adds an independent, stronger signal when a GFA
assembly graph is available: a segment with a link back to itself is a
direct de Bruijn-graph self-loop -- the same "contig links to its own start"
check Recycler/Unicycler use to call circular plasmids.

`obelisk_hunt/cli.py` combines both into one verdict per contig:
`circular_graph+sequence` (both signals agree, highest confidence),
`circular_graph_only`, `circular_sequence_only`, or `linear`.

Explicitly out of scope for this step (later stages of the plan):
database-homology filtering, secondary-structure folding, and the
shuffled-control stability test.

### Validation

Real iHMP data wasn't reachable in this sandbox, so the detector is
validated against synthetic contigs built to mimic exactly what an assembler
does to a circular molecule (duplicate a k-1-ish window across the linearization
break point), plus negative controls:

```bash
pip install -r requirements.txt
python -m pytest tests/ -v
```

10 tests, including a 300-random-contig sweep asserting zero false positives
on pure linear noise. This is a stand-in for a real positive/negative
control, not a replacement for one -- the actual Week 1-3 "recover known
Obelisks in iHMP data" milestone still requires running the real pipeline
above on real reads.
