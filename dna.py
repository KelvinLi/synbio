_nucleotides = {"a" : "t",
                "c" : "g",
                "g" : "c",
                "t" : "a"}

_nucleotide_repr = {frozenset("a")    : "A",
                    frozenset("c")    : "C",
                    frozenset("g")    : "G",
                    frozenset("t")    : "T",
                    frozenset("at")   : "W",
                    frozenset("cg")   : "S",
                    frozenset("ac")   : "M",
                    frozenset("gt")   : "K",
                    frozenset("act")  : "H",
                    frozenset("cgt")  : "B",
                    frozenset("acg")  : "V",
                    frozenset("agt")  : "D",
                    frozenset("acgt") : "N"}

def _min_rotation(sequence):
    """
    Returns the amount of left-rotation of `sequence' that minimizes the
    sequence under lexicographical order. The `n'-step left-rotation is:
        sequence[n:] + sequence[:n]
    """
    L = len(sequence)
    depth = 0
    candidates = tuple(range(L))
    while depth < L and len(candidates) > 1:
        vals = tuple(sequence[(c + depth) % L] for c in candidates)
        lightest = min(vals)
        candidates = tuple(c for c, v in zip(candidates, vals)
                           if v == lightest)
        depth += 1
    return candidates[0]

def _circular_fragment(sequence, start, length):
    """
    Returns the subsequence of given length, beginning at index `start',
    treating `sequence' as circular.
    """
    L = len(sequence)
    if not L or length < 0 or length > L:
        raise IndexError
    start %= L
    end = (start + length) % L
    if start < end:
        return sequence[start:end]
    return sequence[start:] + sequence[:end]

def _has_overlap(starts_ends):
    """
    Given an iterable of the form [(start0, end0), (start1, end1), ...],
    returns whether the half-open intervals have any overlap.
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
    Normalizes the input to `_has_overlap', assuming `starts_ends'
    specifies intervals on a circular sequence of size `L'.
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
    Returns True if `starts_ends' half-open intervals overlap on `sequence'.
    """
    return _has_overlap(_normalized_circular_starts_ends(starts_ends,
                                                         len(sequence))
                        if sequence.is_circular()
                        else starts_ends)

class Nucleotide:
    def __init__(self, wildcard):
        w = frozenset(wildcard)
        assert all(base in _nucleotides for base in w)
        self._wildcard = w

    def __str__(self):
        try:
            return _nucleotide_repr[self._wildcard]
        except KeyError as exc:
            raw_message = "No representation for wildcard ({0})"
            message = raw_message.format(self._wildcard)
            raise ValueError(message) from exc

    __repr__ = __str__

    def __eq__(self, other):
        return self is other

    def is_complement(self, other):
        return other._wildcard == \
               frozenset(_nucleotides[base] for base in self._wildcard)

    def sort_key(self):
        return tuple(2 ** i
                     for i, base in enumerate(sorted(_nucleotides))
                     if base in self._wildcard)

class BaseSequence:
    def __init__(self, nucleotides):
        self._nucleotides = nucleotides

    def __len__(self):
        return len(self._nucleotides)

    def __eq__(self, other):
        return self is other

    def dump(self):
        return "".join(str(n) for n in self._nucleotides)

class LinearSequence(BaseSequence):
    def __init__(self, nucleotides):
        super().__init__(tuple(nucleotides))

    def __str__(self):
        return "linear sequence of {0} bases".format(len(self))

    __repr__ = __str__

    def is_reverse_complement(self, other):
        if other.is_circular():
            return False
        return all(s.is_complement(o)
                   for s, o in zip(reversed(self._nucleotides),
                                   other._nucleotides))

    def is_circular(self):
        return False

    # TODO: transfer properties?
    def fragment(self, start, length):
        if length < 0:
            raise IndexError
        end = start + length
        return LinearSequence(self._nucleotides[start:end])

class CircularSequence(BaseSequence):
    def __init__(self, nucleotides):
        nuc = tuple(nucleotides)
        amount = _min_rotation(tuple(n.sort_key() for n in nuc))
        super().__init__(nuc[amount:] + nuc[:amount])

    def __str__(self):
        return "circular sequence of {0} bases".format(len(self))

    __repr__ = __str__

    def is_circular(self):
        return True

    # TODO: transfer properties?
    def fragment(self, start, length):
        frag = _circular_fragment(self._nucleotides, start, length)
        return LinearSequence(frag)

class _Annealment:
    def __init__(self, sequences, starts, length):
        self.sequences = tuple(sequences)
        self.starts = tuple(starts)
        self._length = length
        self._validate()

    def _validate(self):
        f0, f1 = tuple(seq.fragment(start, self._length)
                       for seq, start in zip(self.sequences, self.starts))
        if not f0.is_reverse_complement(f1):
            raise ValueError("sequences must be reverse complementary over " \
                             "annealment region")

    def ends(self):
        return tuple(s + self._length for s in self.starts)

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

class Cluster:
    def __init__(self):
        self._sequences = list()
        self._annealments = list()

    def add_sequence(self, new):
        if new in self._sequences:
            raise ValueError("cluster already contains this sequence")
        self._sequences.append(new)

    def _validate_pre_add_annealment(self, sequences, starts, length):
        if not len(sequences) == len(starts) == 2:
            raise ValueError("must have two sequences and two start points")
        if not length > 0:
            raise ValueError("length must be positive")
        if not all(s >= 0 for s in starts):
            raise ValueError("start points must be non-negative")
        if not all(s in self._sequences for s in sequences):
            raise ValueError("cluster must contain the input sequences")

    def add_annealment(self, sequences, starts, length, *, overwrite=False):
        sequences = tuple(sequences)
        starts = tuple(starts)
        self._validate_pre_add_annealment(sequences, starts, length)
        new = _Annealment(sequences, starts, length)
        if overwrite:
            raise NotImplementedError
        if any(new.has_overlap(old) for old in self._annealments):
            raise ValueError("refusing to overwrite existing annealment")
        self._annealments.append(new)
