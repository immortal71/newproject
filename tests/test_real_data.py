"""Real-data validation: run the detector on genuine SRA-derived contigs.

See tests/fixtures/README.md for provenance. This is a stand-in for the
Week 1-3 "recover known Obelisks in iHMP data" milestone -- it's a different
real sample (not iHMP, not necessarily an Obelisk), but it is real assembled
sequence, not synthetic, and it turns up a biologically coherent circular
signal rather than noise.

On this fixture the detector flags a cluster of ~18 independently-assembled
contigs whose trimmed length converges on ~337-338 nt -- the known genome
length of Peach Latent Mosaic Viroid, a real circular RNA plant pathogen (the
underlying sample is peach tissue, per VNom's own usage example). It also
separately flags a handful of ~673-674 nt contigs, almost exactly 2x that
unit length: a head-to-tail dimer, the expected rolling-circle-replication
concatemer intermediate for a viroid. Two independent length clusters landing
in a clean ~2x ratio on real data is exactly the kind of internal consistency
check that is hard to get from noise.
"""

from pathlib import Path

from Bio import SeqIO

from obelisk_hunt.circularity import detect_self_overlap

FIXTURE = Path(__file__).parent / "fixtures" / "SRR11060618_subset.fasta"


def test_detector_recovers_a_consistent_circular_unit_length():
    trimmed_lengths = []
    for record in SeqIO.parse(FIXTURE, "fasta"):
        overlap = detect_self_overlap(str(record.seq))
        if overlap is not None:
            trimmed_lengths.append(overlap.trimmed_length)

    monomer_cluster = [n for n in trimmed_lengths if 330 <= n <= 345]
    dimer_cluster = [n for n in trimmed_lengths if 660 <= n <= 690]

    # Independently-assembled contigs converging tightly on one unit length
    # (and a second cluster at ~2x that length) is the signature of a real
    # circular replicon, not a detector artifact.
    assert len(monomer_cluster) >= 10
    assert len(dimer_cluster) >= 2

    monomer_mean = sum(monomer_cluster) / len(monomer_cluster)
    dimer_mean = sum(dimer_cluster) / len(dimer_cluster)
    assert 1.9 <= dimer_mean / monomer_mean <= 2.1
