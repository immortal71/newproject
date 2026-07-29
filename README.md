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

## Status

Week 1 step 1 was *"pull one iHMP metatranscriptome, assemble it, and write
the circularity detector. Nothing else."* That circularity detector is done.
Since then this has grown to cover the other half of the core signature too:
the structure-stability test against a dinucleotide-shuffled null. Database-
homology filtering is the one piece of the signature definition not yet
built.

**Circularity detector, structure-stability scorer: built and tested.**
`obelisk_hunt/circularity.py` + `obelisk_hunt/gfa.py` (circularity, see
below) and `obelisk_hunt/structure.py` + `obelisk_hunt/shuffle.py`
(stability-vs-shuffled-null, see below) all have synthetic-data test suites,
plus one real-data validation (see "Validation").

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

# 3. detect circular contigs, and write out their trimmed sequences
python -m obelisk_hunt detect-circular data/<accession>_assembly/transcripts.fasta \
    --gfa data/<accession>_assembly/assembly_graph_with_scaffolds.gfa \
    --out data/<accession>_assembly/circular_contigs.tsv \
    --trimmed-fasta data/<accession>_assembly/circular_trimmed.fasta

# 4. score each circular candidate against its own dinucleotide-shuffled null
python -m obelisk_hunt score-stability data/<accession>_assembly/circular_trimmed.fasta \
    --n-shuffles 100 --seed 0 \
    --out data/<accession>_assembly/stability.tsv
```

(`detect-circular` is also the default when no subcommand is given, for
backward compatibility with the Week-1 invocation:
`python -m obelisk_hunt contigs.fasta ...`.)

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
`circular_graph_only`, `circular_sequence_only`, or `linear`. Pass
`--trimmed-fasta` to also write out the overlap-trimmed (de-duplicated)
sequence of every circular-flagged contig, ready to feed into
`score-stability`.

## The structure-stability scorer

`obelisk_hunt/shuffle.py` implements an Altschul-Erikson-style dinucleotide
shuffle: it models the sequence as a walk over a multigraph (nodes = the four
bases, edges = the sequence's own dinucleotides) and generates a *different*
Eulerian path through that same multigraph via a randomized Hierholzer walk.
Any such path has identical length, base composition, and dinucleotide
composition to the original -- the matched null the project brief calls for,
not just a same-GC shuffle.

`obelisk_hunt/structure.py` folds a candidate with ViennaRNA (`RNA.fold`),
folds `--n-shuffles` dinucleotide-shuffled versions of it, and reports where
the real sequence's MFE sits in that null distribution: a z-score (positive
= more stable than its shuffles) and a one-sided empirical p-value.

`python -m obelisk_hunt score-stability` runs this over a FASTA of
candidates (typically `detect-circular`'s `--trimmed-fasta` output) and
writes one row per sequence.

This is the actual test in the project's definition: *"does this look
designed by evolution to be a structure, while matching nothing"* -- the
"matching nothing" (homology) half is still not built.

## Validation

Real iHMP data wasn't reachable in this sandbox (see "Status" above), so
both modules are validated against synthetic data built to exercise the
exact failure modes that matter, plus one genuine real-data result:

```bash
pip install -r requirements.txt
python -m pytest tests/ -v
```

20 tests:
- **Circularity** (`test_circularity.py`): synthetic circular/linear
  contigs, including a 300-random-contig sweep asserting zero false
  positives on pure linear noise, and a GFA self-loop parser test.
- **Shuffle** (`test_shuffle.py`): asserts the shuffle exactly preserves
  length, base composition, and dinucleotide composition, and that it
  actually produces different orderings when more than one exists.
- **Stability** (`test_structure.py`): a sequence built as an arbitrary arm
  plus its own reverse complement (an obligate hairpin -- the "designed to
  hold a shape" case) scores z > 10, p <= 0.02 against its shuffles; a
  same-length random sequence scores -3 < z < 3.
- **Real data** (`test_real_data.py`): see below.

### The one real-data result so far

`tests/fixtures/SRR11060618_subset.fasta` (provenance in
`tests/fixtures/README.md`) is a real rnaSPAdes assembly of a real SRA run --
not iHMP, and not necessarily an Obelisk, but genuine sequence, pulled from
a source this sandbox could actually reach. Running the full chain on it:

```bash
python -m obelisk_hunt detect-circular tests/fixtures/SRR11060618_subset.fasta \
    --trimmed-fasta /tmp/trimmed.fasta
python -m obelisk_hunt score-stability /tmp/trimmed.fasta --n-shuffles 200 --seed 0
```

recovers 28/38 contigs as circular, with independently-assembled contigs
converging tightly on two length clusters at a clean ~2x ratio (~337-338nt
monomer, ~673-674nt dimer) -- the known genome length of Peach Latent Mosaic
Viroid and its expected rolling-circle-replication concatemer, on a
peach-tissue sample. `test_real_data.py` asserts that clustering. Whether
those same monomer sequences also clear the stability-vs-shuffled-null bar
is exactly the kind of result worth checking next.
