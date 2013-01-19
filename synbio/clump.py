"""Raw model of a "clump" of DNA."""

def _half_annealments(annealments):
    for ann in annealments:
        for i in 0, 1:
            yield ann.sequences[i], ann.starts[i], ann.length

def _has_overlap_halves(seqs, starts, lengths):
    common_seq = seqs[0]
    if seqs[1] is not common_seq:
        return False
    diff = starts[0] - starts[1]
    if not common_seq.is_circular:
        return -lengths[0] < diff < lengths[1]
    L = len(common_seq)
    assert all(0 < x <= L for x in lengths)
    return diff % L < lengths[1] or diff % L > L - lengths[0]

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

    def has_overlap(self, others):
        return any(_has_overlap_halves(*zip(self_half, other_half))
                   for self_half in _half_annealments((self,))
                   for other_half in _half_annealments(others))

class Clump:
    def __init__(self):
        self.sequences = tuple()
        self.annealments = tuple()

    def __replace(self, *, seq=None, ann=None):
        new = Clump()
        new.sequences = self.sequences if seq is None else tuple(seq)
        new.annealments = self.annealments if ann is None else tuple(ann)
        return new

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
        return self.__replace(seq=self.sequences + (new,))

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
        if new.has_overlap(self.annealments):
            raise ValueError("refusing to overwrite existing annealment")
        return self.__replace(ann=self.annealments + (new,))

    def strip_annealments(self, sequences):
        """Strip all annealments referring to any sequence in `sequences'.

        sequences -- (iterable)
        """
        S = tuple(sequences)
        i = (ann for ann in self.annealments
             if not any(seq in ann.sequences for seq in S))
        return self.__replace(ann=i)

    def remove_sequence(self, old):
        """Remove sequence `old' and strip all annealments referring to `old'.

        old -- sequence to be removed

        Raises ValueError if `old' is not in this Clump.
        """
        new = self.strip_annealments((old,))
        return new.__replace(seq=(s for s in new.sequences if s is not old))
