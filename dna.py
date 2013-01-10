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
    size = len(sequence) # constant
    depth = 0
    candidates = tuple(range(size))
    while depth < size and len(candidates) > 1:
        vals = tuple(sequence[(c + depth) % size] for c in candidates)
        lightest = min(vals)
        candidates = tuple(c for c, v in zip(candidates, vals) if v == lightest)
        depth += 1
    return candidates[0]

# immutable
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

    def complement(self):
        return Nucleotide(_nucleotides[base]
                          for base in self._wildcard)

    def sort_key(self):
        return tuple(2 ** i
                     for i, base in enumerate(sorted(_nucleotides))
                     if base in self._wildcard)

class BaseSequence:
    def __init__(self, nucleotides):
        self._nucleotides = nucleotides

    def __len__(self):
        return len(self._nucleotides)

    def dump(self):
        return "".join(str(n) for n in self._nucleotides)

class LinearSequence(BaseSequence):
    def __init__(self, nucleotides):
        super.__init__(tuple(nucleotides))

    def __str__(self):
        return "linear sequence of {0} bases".format(len(self))

    __repr__ = __str__

class CircularSequence(BaseSequence):
    def __init__(self, nucleotides):
        amount = _min_rotation(n.sort_key() for n in nucleotides)
        super.__init__(tuple(nucleotides[amount:] + nucleotides[:amount]))

    def __str__(self):
        return "circular sequence of {0} bases".format(len(self))

    __repr__ = __str__
