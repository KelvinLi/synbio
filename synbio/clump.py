"""Raw model of a "clump" of DNA."""

def _has_overlap(common_seq, starts, lengths):
    diff = starts[0] - starts[1]
    if not common_seq.is_circular:
        return -lengths[0] < diff < lengths[1]
    L = len(common_seq)
    assert all(0 < x <= L for x in lengths)
    return diff % L < lengths[1] or diff % L > L - lengths[0]

class _Annealment:
    def __init__(self, keys, starts, length):
        self.keys = keys
        self.starts = starts
        self.length = length

class AnnealmentQueryResult:
    def __init__(self, start, length, other_key):
        self.start = start
        self.length = length
        self.other_key = other_key

class Clump:
    def __init__(self):
        self._sequences = tuple()
        self._annealments = tuple()

    def __append_sequence(self, sequence):
        new = Clump()
        new._sequences = self._sequences + (sequence,)
        new._annealments = self._annealments
        return new

    def __replace_annealments(self, replacement):
        new = Clump()
        new._sequences = self._sequences
        new._annealments = replacement
        return new

    def _validate_annealment(self, annealment):
        f0, f1 = (self._sequences[key].fragment(start, annealment.length)
                  for key, start in zip(annealment.keys, annealment.starts))
        if not f0.reverse_complement().has_same_nucleotides(f1):
            raise ValueError("sequences must be reverse complementary over "
                             "annealment region")

    def _has_annealment_conflict(self, annealment):
        return any(_has_overlap(self._sequences[key],
                                (r.start, start),
                                (r.length, annealment.length))
                   for key, start in zip(annealment.keys, annealment.starts)
                   for r in self.query_annealments(key))

    def add_annealment(self, keys, starts, length, *, overwrite=False):
        keys = tuple(keys)
        starts = tuple(starts)
        if not len(keys) == len(starts) == 2:
            raise ValueError("must have two keys and two start points")
        if not length > 0:
            raise ValueError("length must be positive")
        if not all(isinstance(s, int) for s in starts):
            raise TypeError("start points must be integers")
        if not all(s >= 0 for s in starts):
            raise ValueError("start points must be non-negative")
        if not all(isinstance(k, int) for k in keys):
            raise TypeError("keys must be integers")
        if not all(0 <= k < len(self._sequences) for k in keys):
            raise ValueError("invalid keys")
        new = _Annealment(keys, starts, length)
        if overwrite:
            raise NotImplementedError
        if self._has_annealment_conflict(new):
            raise ValueError("refusing to overwrite conflicting annealment(s)")
        self._validate_annealment(new)
        return self.__replace_annealments(self._annealments + (new,))

    def add_sequence(self, new):
        return self.__append_sequence(new)

    def sequences(self):
        return tuple(self._sequences)

    def locate_sequences(self, match):
        return (key for key, seq in enumerate(self._sequences)
                if seq.has_same_nucleotides(match))

    # TODO: this function is nested too deeply
    def query_annealments(self, key):
        for ann in self._annealments:
            for i, o in (0, 1), (1, 0):
                if ann.keys[i] != key:
                    continue
                a = (ann.starts[i], ann.length, ann.keys[o])
                yield AnnealmentQueryResult(*a)
