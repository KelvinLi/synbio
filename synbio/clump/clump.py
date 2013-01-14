"""Raw model of a "clump" of DNA.
Only the shape layer should be intimately aware of this layer.

The objects in this module assemble together as follows:

    Clump
      |-- Annealment
      |-- BaseSequence
            |-- Nucleotide

BaseSequence is abstract. The concrete subclasses are LinearSequence
and CircularSequence.

The only mutators in this module are in the Clump class. All object
properties are intended to be read-only.
"""

from . import priv

class Nucleotide:
    def __init__(self, wildcard):
        w = frozenset(wildcard)
        if not all(base in priv.nucleotides for base in w):
            raise ValueError
        self.wildcard = w

    def __str__(self):
        try:
            return priv.nucleotide_repr[self.wildcard]
        except KeyError as exc:
            raw_message = "No representation for wildcard ({0})"
            message = raw_message.format(self.wildcard)
            raise ValueError(message) from exc

    __repr__ = __str__

    def __eq__(self, other):
        return self is other

    def is_complement(self, other):
        return other.wildcard == \
               frozenset(priv.nucleotides[base] for base in self.wildcard)

    def sort_key(self):
        return tuple(2 ** i
                     for i, base in enumerate(sorted(priv.nucleotides))
                     if base in self.wildcard)

class BaseSequence:
    def __init__(self, nucleotides):
        self.nucleotides = nucleotides

    def __len__(self):
        return len(self.nucleotides)

    def __eq__(self, other):
        return self is other

    def dump(self):
        return "".join(str(n) for n in self.nucleotides)

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
                   for s, o in zip(reversed(self.nucleotides),
                                   other.nucleotides))

    def is_circular(self):
        return False

    # TODO: transfer properties?
    def fragment(self, start, length):
        if length < 0:
            raise IndexError
        end = start + length
        return LinearSequence(self.nucleotides[start:end])

class CircularSequence(BaseSequence):
    def __init__(self, nucleotides):
        nuc = tuple(nucleotides)
        amount = priv.min_rotation(tuple(n.sort_key() for n in nuc))
        super().__init__(nuc[amount:] + nuc[:amount])

    def __str__(self):
        return "circular sequence of {0} bases".format(len(self))

    __repr__ = __str__

    def is_circular(self):
        return True

    # TODO: transfer properties?
    def fragment(self, start, length):
        frag = priv.circular_fragment(self.nucleotides, start, length)
        return LinearSequence(frag)

class Annealment:
    def __init__(self, sequences, starts, length):
        self.sequences = tuple(sequences)
        self.starts = tuple(starts)
        self.length = length
        self._validate()

    def _validate(self):
        f0, f1 = tuple(seq.fragment(start, self.length)
                       for seq, start in zip(self.sequences, self.starts))
        if not f0.is_reverse_complement(f1):
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
            if priv.sequence_has_overlap(common_sequence, common_starts_ends):
                return True
        return False

class Clump:
    def __init__(self):
        self.sequences = list()
        self.annealments = list()

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
        self.sequences.append(new)

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
        self.annealments.append(new)

    def strip_annealments(self, sequences):
        """
        Strip all annealments referring to any sequence in `sequences'.

        sequences -- (iterable of BaseSequence)
        """
        S = tuple(sequences)
        self.annealments = [ann for ann in self.annealments
                            if not any(seq in ann.sequences for seq in S)]

    def remove_sequence(self, old):
        """
        Remove `old' from this Clump, and strip all annealments
        referring to `old'.

        old -- (BaseSequence) sequence to be removed

        Raises ValueError if `old' is not in this Clump.
        """
        self.strip_annealments((old,))
        self.sequences.remove(old)
        assert old not in self.sequences
