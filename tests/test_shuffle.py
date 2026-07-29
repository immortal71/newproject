import random
from collections import Counter

from obelisk_hunt.shuffle import dinucleotide_shuffle


def dinuc_counts(seq: str) -> Counter:
    return Counter(seq[i:i + 2] for i in range(len(seq) - 1))


class TestDinucleotideShuffle:
    def test_preserves_length_composition_and_dinucleotide_counts(self):
        rng = random.Random(0)
        seq = "".join(rng.choice("ACGT") for _ in range(2000))

        for trial_seed in range(20):
            shuffled = dinucleotide_shuffle(seq, random.Random(trial_seed))

            assert len(shuffled) == len(seq)
            assert Counter(shuffled) == Counter(seq)
            assert dinuc_counts(shuffled) == dinuc_counts(seq)

    def test_preserves_first_and_last_base(self):
        seq = "ACGACGTACGGGT"
        for trial_seed in range(20):
            shuffled = dinucleotide_shuffle(seq, random.Random(trial_seed))
            assert shuffled[0] == seq[0]
            assert shuffled[-1] == seq[-1]

    def test_actually_reshuffles_when_multiple_paths_exist(self):
        # a sequence with real branching in its dinucleotide graph should
        # produce more than one distinct shuffle across many trials
        seq = "ACGACGTACGACGTACGACGTACGT" * 3
        outputs = {dinucleotide_shuffle(seq, random.Random(s)) for s in range(50)}
        assert len(outputs) > 1

    def test_handles_repetitive_sequence(self):
        seq = "AAAAAAAAAA"
        shuffled = dinucleotide_shuffle(seq, random.Random(1))
        assert shuffled == seq  # only one possible Eulerian path

    def test_handles_short_sequences(self):
        assert dinucleotide_shuffle("A", random.Random(0)) == "A"
        assert dinucleotide_shuffle("AC", random.Random(0)) == "AC"
