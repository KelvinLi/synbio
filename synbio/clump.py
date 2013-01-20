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

def _query_pair_match(q, x):
    assert len(q) == len(x) == 2
    q = list(q)
    for i in 0, 1:
        if q[i] is None:
            q[i] = x[i]
    return x if tuple(q) == tuple(x) else None

class AnnealmentQuery:
    def __init__(self, sequences, starts, length):
        self.sequences = tuple(sequences)
        self.starts = tuple(starts)
        self.length = length
        if not len(self.sequences) == len(self.starts) == 2:
            raise ValueError

    def _reversed(self):
        return AnnealmentQuery(reversed(self.sequences),
                               reversed(self.starts),
                               self.length)

    def _match_try(self, annealment):
        seqs_match = _query_pair_match(self.sequences, annealment.sequences)
        if seqs_match is None:
            return None
        starts_match = _query_pair_match(self.starts, annealment.starts)
        if starts_match is None:
            return None
        if self.length is not None and self.length != annealment.length:
            return None
        return AnnealmentQuery(seqs_match, starts_match, annealment.length)

    def _match_one(self, annealment):
        m = self._match_try(annealment)
        if m is None:
            return self._reversed()._match_try(annealment)
        return m

    def match_all(self, annealments):
        return filter(None, (self._match_one(ann) for ann in annealments))

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

    def __replace_sequences(self, replacement):
        new = Clump()
        new.sequences = tuple(replacement)
        new.annealments = self.annealments
        return new

    def __replace_annealments(self, replacement):
        new = Clump()
        new.sequences = self.sequences
        new.annealments = tuple(replacement)
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
        return self.__replace_sequences(self.sequences + (new,))

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
        return self.__replace_annealments(self.annealments + (new,))

    def strip_annealments(self, sequences):
        """Strip all annealments referring to any sequence in `sequences'.

        sequences -- (iterable)
        """
        return self.__replace_annealments(
            ann for ann in self.annealments
            if not any(seq in ann.sequences for seq in tuple(sequences)))

    def remove_sequence(self, old):
        """Remove sequence `old' and strip all annealments referring to `old'.

        old -- sequence to be removed

        Raises ValueError if `old' is not in this Clump.
        """
        new = self.strip_annealments((old,))
        return new.__replace_sequences(s for s in new.sequences
                                       if s is not old)

    def query_annealments(self, query):
        return query.match_all(self.annealments)
