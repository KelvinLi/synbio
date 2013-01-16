from .clump import clump
from .clump import sequence
from . import shape

# = = = =
GENERIC = shape.ShapeType()
@shape.declare_op(GENERIC, shape.OpNames.VALIDATE)
def _x(clump_obj):
    pass

@shape.declare_op(GENERIC, "count_sequences")
def _x(clump_obj):
    return len(clump_obj.sequences)

# = = = =
LINEAR = shape.ShapeType(GENERIC)
@shape.declare_op(LINEAR, shape.OpNames.VALIDATE)
def _x(clump_obj):
    if not all(not s.is_circular for s in clump_obj.sequences):
        raise shape.ShapeError("all sequences must be linear")

# = = = =
DOUBLE_STRANDED = shape.ShapeType(GENERIC)
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
PCR_TEMPLATE = shape.ShapeType(LINEAR, DOUBLE_STRANDED)
@shape.declare_op(PCR_TEMPLATE, shape.OpNames.VALIDATE)
def _x(clump_obj):
    pass

@shape.declare_op(PCR_TEMPLATE, shape.OpNames.INIT)
def _x(clump_obj, nucleotides):
    seq = sequence.LinearSequence(nucleotides)
    rseq = seq.reverse_complement()
    clump_obj.add_sequence(seq)
    clump_obj.add_sequence(rseq)
    clump_obj.add_annealment((seq, rseq), (0, 0), len(seq))
