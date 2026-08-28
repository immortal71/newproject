import random
import shutil

import pytest

from obelisk_hunt.homology import check_homology

pytestmark = pytest.mark.skipif(
    shutil.which("minimap2") is None,
    reason="minimap2 not installed (bioconda `minimap2`, or a static binary from "
           "https://github.com/lh3/minimap2/releases)",
)


def random_seq(n: int, seed: int) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


def write_fasta(path, records):
    with open(path, "w") as fh:
        for name, seq in records:
            fh.write(f">{name}\n{seq}\n")


class TestCheckHomology:
    def test_exact_match_is_flagged_as_homologous(self, tmp_path):
        known = random_seq(500, seed=1)
        ref = tmp_path / "ref.fasta"
        write_fasta(ref, [("known_element", known)])

        candidates = tmp_path / "candidates.fasta"
        write_fasta(candidates, [("hit_candidate", known)])

        results = check_homology(str(candidates), str(ref))

        assert results["hit_candidate"].has_hit
        assert results["hit_candidate"].best_identity > 0.99
        assert results["hit_candidate"].best_target == "known_element"

    def test_unrelated_sequence_has_zero_homology(self, tmp_path):
        ref = tmp_path / "ref.fasta"
        write_fasta(ref, [("known_element", random_seq(500, seed=1))])

        candidates = tmp_path / "candidates.fasta"
        write_fasta(candidates, [("novel_candidate", random_seq(400, seed=2))])

        results = check_homology(str(candidates), str(ref))

        assert not results["novel_candidate"].has_hit

    def test_partial_short_overlap_below_coverage_threshold_not_flagged(self, tmp_path):
        # a candidate that only shares a short coincidental stretch with the
        # reference (well below --min-coverage) should not count as a hit
        ref_seq = random_seq(600, seed=3)
        ref = tmp_path / "ref.fasta"
        write_fasta(ref, [("known_element", ref_seq)])

        # 40 bases borrowed from the reference embedded in an otherwise
        # unrelated 500bp sequence: real but tiny fractional coverage
        candidate_seq = random_seq(230, seed=4) + ref_seq[100:140] + random_seq(230, seed=5)
        candidates = tmp_path / "candidates.fasta"
        write_fasta(candidates, [("mostly_novel", candidate_seq)])

        results = check_homology(str(candidates), str(ref), min_coverage=0.5)

        assert not results["mostly_novel"].has_hit

    def test_mixed_batch_reports_each_independently(self, tmp_path):
        known = random_seq(400, seed=6)
        ref = tmp_path / "ref.fasta"
        write_fasta(ref, [("known_element", known)])

        candidates = tmp_path / "candidates.fasta"
        write_fasta(candidates, [
            ("matches_known", known),
            ("totally_novel", random_seq(400, seed=7)),
        ])

        results = check_homology(str(candidates), str(ref))

        assert results["matches_known"].has_hit
        assert not results["totally_novel"].has_hit
