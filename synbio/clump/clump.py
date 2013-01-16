"""Raw model of a "clump" of DNA."""

from . import priv

class Annealment:
    def __init__(self, sequences, starts, length):
        self.sequences = tuple(sequences)
        self.starts = tuple(starts)
        self.length = length
        self._validate()

    def _validate(self):
        f0, f1 = tuple(seq.fragment(start, self.length)
                       for seq, start in zip(self.sequences, self.starts))
        if not f0.reverse_complement().has_same_nucleotides(f1):
            raise ValueError("sequences must be reverse complementary over " \
                             "annealment region")

    def ends(self):
        return tuple(s + self.length for s in self.starts)

    def has_overlap(self, other):
        # TODO: there is room for optimization
        for s, o in (0, 0), (0, 1), (1, 0), (1, 1):
            common_sequence = self.sequences[s]
            if common_sequence is not other.sequences[o]:
                continue
            common_starts_ends = ((self.starts[s], self.ends()[s]),
                                  (other.starts[o], other.ends()[o]))
            if priv.sequence_has_overlap(common_sequence, common_starts_ends):
                return True
        return False

class Clump:
    def __init__(self):
        self.sequences = list()
        self.annealments = list()

    def _validate_pre_add_annealment(self, sequences, starts, length):
        if not len(sequences) == len(starts) == 2:
            raise ValueError("must have two sequences and two start points")
        if not length > 0:
            raise ValueError("length must be positive")
        if not all(s >= 0 for s in starts):
            raise ValueError("start points must be non-negative")
        if not all(s in self.sequences for s in sequences):
            raise ValueError("clump must contain the input sequences")

    def add_sequence(self, new):
        if new in self.sequences:
            raise ValueError("clump already contains this sequence")
        self.sequences.append(new)

    def add_annealment(self, sequences, starts, length, *, overwrite=False):
        """
        Declares annealment between nucleotides of two sequences, with
        two given start positions.

        length    -- (int) number of nucleotides annealed
        overwrite -- (boolean) allow changing of existing annealments
        """
        sequences = tuple(sequences)
        starts = tuple(starts)
        self._validate_pre_add_annealment(sequences, starts, length)
        new = Annealment(sequences, starts, length)
        if overwrite:
            raise NotImplementedError
        if any(new.has_overlap(old) for old in self.annealments):
            raise ValueError("refusing to overwrite existing annealment")
        self.annealments.append(new)

    def strip_annealments(self, sequences):
        """
        Strip all annealments referring to any sequence in `sequences'.

        sequences -- (iterable)
        """
        S = tuple(sequences)
        self.annealments = [ann for ann in self.annealments
                            if not any(seq in ann.sequences for seq in S)]

    def remove_sequence(self, old):
        """
        Remove `old' from this Clump, and strip all annealments
        referring to `old'.

        old -- sequence to be removed

        Raises ValueError if `old' is not in this Clump.
        """
        self.strip_annealments((old,))
        self.sequences.remove(old)
        assert old not in self.sequences
