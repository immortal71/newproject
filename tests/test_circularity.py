"""Synthetic-data validation of the circularity detector.

We can't reach NCBI/EBI from this sandbox (see README), so this suite is the
stand-in positive/negative control until the detector is run on a real
assembly. It builds contigs that mimic exactly what a linear assembler does
to a circular molecule (duplicate a k-1 window across the break point) and
checks the detector recovers the closure -- plus negative controls that must
NOT be flagged.
"""

import random

import pytest

from obelisk_hunt.circularity import detect_self_overlap
from obelisk_hunt.gfa import parse_gfa


def random_seq(n: int, seed: int) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


def mutate(seq: str, n_subs: int, seed: int) -> str:
    rng = random.Random(seed)
    seq = list(seq)
    positions = rng.sample(range(len(seq)), min(n_subs, len(seq)))
    for p in positions:
        seq[p] = rng.choice([b for b in "ACGT" if b != seq[p]])
    return "".join(seq)


def linearize_circular(core: str, overlap_k: int) -> str:
    """Simulate a de-Bruijn-style linearization of a circular molecule:
    the assembler's path re-reads its own first `overlap_k` bases at the end."""
    return core + core[:overlap_k]


class TestSelfOverlap:
    def test_perfect_circular_closure_detected(self):
        core = random_seq(1000, seed=1)
        contig = linearize_circular(core, overlap_k=63)  # e.g. k=64 de Bruijn overlap

        result = detect_self_overlap(contig)

        assert result is not None
        assert result.overlap_len == 63
        assert result.identity == 1.0
        assert result.trimmed_length == 1000

    def test_circular_closure_with_sequencing_errors_still_detected(self):
        core = random_seq(1000, seed=2)
        contig = linearize_circular(core, overlap_k=80)
        # introduce a few mismatches only within the duplicated join region
        contig = contig[:1000] + mutate(contig[1000:], n_subs=3, seed=3)

        result = detect_self_overlap(contig)

        assert result is not None
        assert result.identity >= 0.9
        assert abs(result.trimmed_length - 1000) <= 5

    def test_pure_linear_sequence_not_flagged(self):
        contig = random_seq(1000, seed=4)

        result = detect_self_overlap(contig)

        assert result is None

    def test_internal_repeat_not_mistaken_for_circularity(self):
        # a repeat that sits in the middle of the contig, not at the termini,
        # must not be called circular
        core = random_seq(1000, seed=5)
        repeat = core[100:150]
        contig = core[:400] + repeat + core[400:]

        result = detect_self_overlap(contig)

        assert result is None

    def test_short_contig_below_window_is_skipped(self):
        contig = random_seq(30, seed=6)

        result = detect_self_overlap(contig, min_overlap=20)

        assert result is None

    def test_overlap_shorter_than_min_overlap_not_flagged(self):
        core = random_seq(1000, seed=7)
        contig = linearize_circular(core, overlap_k=5)

        result = detect_self_overlap(contig, min_overlap=20)

        assert result is None

    def test_too_many_mismatches_in_join_rejected(self):
        core = random_seq(1000, seed=8)
        contig = linearize_circular(core, overlap_k=100)
        tail = mutate(contig[1000:], n_subs=40, seed=9)  # ~40% mismatches
        contig = contig[:1000] + tail

        result = detect_self_overlap(contig, max_edit_frac=0.10)

        assert result is None


class TestFalsePositiveRate:
    def test_no_false_positives_across_many_random_linear_contigs(self):
        """A metatranscriptome assembly has thousands of linear contigs; the
        detector's false-positive rate on pure noise needs to be ~0, not just
        0 on one lucky seed."""
        false_positives = 0
        n = 300
        for seed in range(n):
            length = random.Random(seed).randint(150, 3000)
            contig = random_seq(length, seed=100_000 + seed)
            if detect_self_overlap(contig) is not None:
                false_positives += 1

        assert false_positives == 0, f"{false_positives}/{n} random linear contigs flagged circular"


class TestGfaSelfLoop:
    def test_self_loop_detected(self, tmp_path):
        gfa = tmp_path / "graph.gfa"
        gfa.write_text(
            "H\tVN:Z:1.0\n"
            "S\tNODE_1\t*\tLN:i:1000\n"
            "S\tNODE_2\t*\tLN:i:500\n"
            "L\tNODE_1\t+\tNODE_1\t+\t63M\n"
            "L\tNODE_1\t+\tNODE_2\t+\t0M\n"
        )

        graph = parse_gfa(str(gfa))

        assert graph.has_self_loop("NODE_1")
        assert not graph.has_self_loop("NODE_2")

    def test_no_self_loops_in_linear_graph(self, tmp_path):
        gfa = tmp_path / "graph.gfa"
        gfa.write_text(
            "H\tVN:Z:1.0\n"
            "S\tNODE_1\t*\tLN:i:1000\n"
            "S\tNODE_2\t*\tLN:i:500\n"
            "L\tNODE_1\t+\tNODE_2\t+\t0M\n"
        )

        graph = parse_gfa(str(gfa))

        assert not graph.has_self_loop("NODE_1")
        assert not graph.has_self_loop("NODE_2")


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-v"]))
