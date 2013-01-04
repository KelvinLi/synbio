_nucleotide_complements = {"a" : "t",
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

# immutable
class Nucleotide:
    def __init__(self, wildcard):
        w = frozenset(wildcard)
        assert(all(base in _nucleotide_complements for base in w))
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
        return Nucleotide(_nucleotide_complements[base]
                          for base in self._wildcard)

# immutable
class Sequence:
    def __init__(self, nucleotides):
        self._nucleotides = tuple(nucleotides)

    def __len__(self):
        return len(self._nucleotides)

    def complement(self):
        return Sequence(nucleotide.complement()
                        for nucleotide in self._nucleotides)

    def reverse(self):
        return Sequence(list(self._nucleotides).reverse())

# immutable
class LinearSequence(Sequence):
    def __str__(self):
        return "".join(str(nucleotide) for nucleotide in self._nucleotides)

    __repr__ = __str__

# immutable
class CircularSequence(Sequence):
    pass
