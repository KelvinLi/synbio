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

    def sort_key(self):
        raise NotImplementedError

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
        return Sequence(reversed(self._nucleotides))

# immutable
class LinearSequence(Sequence):
    def __str__(self):
        return "".join(str(nucleotide) for nucleotide in self._nucleotides)

    __repr__ = __str__

_canonical_rotate(sequence, key=None):
    """
    Finds the rotation of the sequence that lexicographically minimizes the key
    function over the sequence. Returns 2-tuple (rotation, amount), where
    rotation is a list giving the rotated sequence, and amount is an integer
    indicating the amount of rotation to the left.

    One in-place left-rotation of a list 'L' is defined as:
        L.append(L.pop(0))

    If the sequence is periodic with period 'n', then the returned
    rotation amount will be less than n.

    As a corollary, if the key function satisfies the proposition:
        key(x) == key(y) iff x == y
    then the returned rotation is "canonical".
    """
    seq = list(sequence)       # constant
    keys = list(map(key, seq)) # constant
    skip = 0                   # which item we are comparing
    candidates = list(seq)     # indices of seq, equivalent to rotation amount

    while skip < len(seq) and len(candidates) > 1:
        # validate candidates
        candidate_keys = [keys[(c + skip) % len(seq)] for c in candidates]
        lightest = min(candidate_keys)

        # prune bad candidates (the ones that don't have the minimum key)
        # ...
        # look at next item

# immutable
class CircularSequence(Sequence):
    def __init__(self, nucleotides):
        rotated, amount = _canonical_rotate(nucleotides,
                                            lambda n : n.sort_key())
        super.__init__(rotated)
