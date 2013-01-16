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

nucleotides = {"a" : "t",
                "c" : "g",
                "g" : "c",
                "t" : "a"}

nucleotide_repr = {frozenset("a")    : "A",
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

def sequence_has_overlap(sequence, starts_ends):
    """
    Returns True if the given half-open intervals overlap on `sequence'.

    sequence -- (BaseSequence) used to extract geometry
    starts_ends -- (iterable) half-open intervals; interpretation
                   depends on geometry
    """
    return _has_overlap(_normalized_circular_starts_ends(starts_ends,
                                                         len(sequence))
                        if sequence.is_circular
                        else starts_ends)
