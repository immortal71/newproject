import random

from obelisk_hunt.structure import fold_mfe, score_stability


def revcomp(seq: str) -> str:
    comp = {"A": "T", "T": "A", "C": "G", "G": "C"}
    return "".join(comp[c] for c in reversed(seq))


class TestFoldMfe:
    def test_known_hairpin_folds_as_expected(self):
        mfe = fold_mfe("GGGGAAACCCC")
        assert mfe == -4.5

    def test_accepts_dna_alphabet(self):
        # contigs come out of assemblers as ACGT; the biological molecule is
        # RNA. Folding must not silently mis-handle a literal 'T'.
        assert fold_mfe("GGGGAAACCCC") == fold_mfe("GGGGAAACCCC".replace("T", "U"))


class TestScoreStability:
    def test_obligate_hairpin_is_far_more_stable_than_its_shuffles(self):
        """A sequence built from an arbitrary arm + its exact reverse
        complement is exactly the "designed to hold a shape" case: any
        dinucleotide shuffle keeps the same composition but (near-certainly)
        destroys the exact positional complementarity, so it should fold far
        worse than the real sequence."""
        rng = random.Random(42)
        arm = "".join(rng.choice("ACGT") for _ in range(40))
        hairpin = arm + "AAAAAA" + revcomp(arm)

        result = score_stability(hairpin, n_shuffles=100, seed=1)

        assert result.z_score > 10
        assert result.empirical_p <= 0.02

    def test_unstructured_random_sequence_looks_like_its_shuffles(self):
        rng = random.Random(43)
        seq = "".join(rng.choice("ACGT") for _ in range(86))

        result = score_stability(seq, n_shuffles=100, seed=1)

        assert -3 < result.z_score < 3
