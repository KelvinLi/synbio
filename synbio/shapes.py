from .clump import clump
from . import shape

# = = = =
GENERIC = shape.Shape()
@shape.declare_op(GENERIC, shape.OpNames.VALIDATE)
def _x(clump_obj):
    pass

@shape.declare_op(GENERIC, "count_sequences")
def _x(clump_obj):
    return len(clump_obj.sequences)

# = = = =
LINEAR = shape.Shape((GENERIC,))
@shape.declare_op(LINEAR, shape.OpNames.VALIDATE)
def _x(clump_obj):
    if not all(not s.is_circular() for s in clump_obj.sequences):
        raise shape.ShapeError("all sequences must be linear")

# = = = =
DOUBLE_STRANDED = shape.Shape((GENERIC,))
@shape.declare_op(DOUBLE_STRANDED, shape.OpNames.VALIDATE)
def _x(clump_obj):
    if not len(clump_obj.sequences) == 2:
        raise shape.ShapeError("must have exactly two sequences")
    if not clump_obj.annealments:
        raise shape.ShapeError("must have at least one annealment")
    if not any(all(seq in ann.sequences
                   for seq in clump_obj.sequences)
               for ann in clump_obj.annealments):
        raise shape.ShapeError("the two sequences must be annealed "
                               "to each other")

@shape.declare_op(DOUBLE_STRANDED, "sequence_lengths")
def _x(clump_obj):
    return tuple(len(seq) for seq in clump_obj.sequences)

# = = = =
PCR_TEMPLATE = shape.Shape((LINEAR, DOUBLE_STRANDED))
@shape.declare_op(PCR_TEMPLATE, shape.OpNames.VALIDATE)
def _x(clump_obj):
    pass

@shape.declare_op(PCR_TEMPLATE, shape.OpNames.INIT)
def _x(clump_obj, nucleotide_string):
    nuc = tuple(clump.Nucleotide(c) for c in nucleotide_string)
    revcomp_nuc = tuple(clump.Nucleotide(clump.priv.nucleotides[base]
                                         for base in n.wildcard)
                        for n in reversed(nuc))
    seqs = tuple(clump.LinearSequence(n) for n in (nuc, revcomp_nuc))
    for s in seqs:
        clump_obj.add_sequence(s)
    clump_obj.add_annealment(seqs, (0, 0), len(nuc))
