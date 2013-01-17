def _min_rotation(sequence):
    """Returns the amount of left-rotation of `sequence' that minimizes the
    sequence under lexicographical order.

    The n-step left-rotation is defined as:
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
    """Returns a subsequence, treating `sequence' as circular.

    start  -- (int) starting index
    length -- (int) number of nucleotides in resulting subsequence
    """
    L = len(sequence)
    if not L or length < 0 or length > L:
        raise IndexError
    start %= L
    end = (start + length) % L
    if start < end:
        return sequence[start:end]
    return sequence[start:] + sequence[:end]

class Nucleotide:
    """A nucleotide is an integer."""
    T, G, C, A = 8, 4, 2, 1

    _complement_table = (
        0, # illegal anyways

              T,
            G  ,
            G|T,
          C    ,
          C|  T,
          C|G  ,
          C|G|T,
        A      ,
        A|    T,
        A|  G  ,
        A|  G|T,
        A|C    ,
        A|C  |T,
        A|C|G  ,
        A|C|G|T,
    )

    @classmethod
    def complement(cls, n):
        assert isinstance(n, int)
        assert 0 < n < 16
        return cls._complement_table[n]

class BaseSequence:
    def __init__(self, nucleotides):
        assert not issubclass(BaseSequence, type(self))
        self.nucleotides = nucleotides

    def __len__(self):
        return len(self.nucleotides)

    def __eq__(self, other):
        return self is other

    def has_same_nucleotides(self, other):
        if self.is_circular ^ other.is_circular:
            raise TypeError
        return self.nucleotides == other.nucleotides

    def reverse_complement(self):
        i = (Nucleotide.complement(n) for n in reversed(self.nucleotides))
        if self.is_circular:
            return CircularSequence(i)
        return LinearSequence(i)

class LinearSequence(BaseSequence):
    is_circular = False

    def __init__(self, nucleotides):
        super().__init__(tuple(nucleotides))

    def fragment(self, start, length):
        if length < 0:
            raise IndexError
        end = start + length
        return LinearSequence(self.nucleotides[start:end])

    def to_circular(self):
        return CircularSequence(self.nucleotides)

class CircularSequence(BaseSequence):
    is_circular = True

    def __init__(self, nucleotides):
        nuc = tuple(nucleotides)

        # exploit the fact that nucleotides are integers
        amount = _min_rotation(nuc)
        super().__init__(nuc[amount:] + nuc[:amount])

    def fragment(self, start, length):
        frag = _circular_fragment(self.nucleotides, start, length)
        return LinearSequence(frag)

    def to_linear(self):
        return self.fragment(0, len(self))
