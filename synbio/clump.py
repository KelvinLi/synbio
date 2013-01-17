"""Raw model of a "clump" of DNA."""

def _has_overlap(starts_ends):
    """
    Returns True if any of the given half-open intervals intersect.

    starts_ends -- (iterable) half-open intervals of the form:
                   [(start0, end0), (start1, end1), ...]
    """
    taken = list()
    for s, e in starts_ends:
        if not s < e:
            raise ValueError
        if any(not (old_e <= s or e <= old_s) for old_s, old_e in taken):
            return True
        taken.append((s, e))
    return False

def _normalized_circular_starts_ends(starts_ends, L):
    """
    Given a set of intervals, normalizes the intervals over modulo arithmetic
    so that intersections can be computed using `_has_overlap'

    starts_ends -- (iterable) half-open intervals, to be normalized modulo `L'
    L           -- (int) modulus
    """
    iter = ((s % L, e % L) for s, e in starts_ends)
    for s, e in iter:
        if s < e:
            yield (s, e)
            continue
        yield (s, L)
        if e:
            yield (0, e)

def _sequence_has_overlap(sequence, starts_ends):
    """
    Returns True if the given half-open intervals overlap on `sequence'.

    sequence -- (BaseSequence) used to extract geometry
    starts_ends -- (iterable) half-open intervals; interpretation
                   depends on geometry
    """
    if sequence.is_circular:
        starts_ends = _normalized_circular_starts_ends(
                      starts_ends, len(sequence))
    return _has_overlap(starts_ends)

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
            if _sequence_has_overlap(common_sequence, common_starts_ends):
                return True
        return False

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
        if any(new.has_overlap(old) for old in self.annealments):
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
