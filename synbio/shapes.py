from . import clump, sequence
from .shape import BaseShapeFunc, BaseShapeType, ShapeError

class _GenericShapeFunc(BaseShapeFunc):
    def count_sequences(self):
        return len(self.c.sequences)

class _GenericShapeType(BaseShapeType):
    def accept(self, clump_obj): return _GenericShapeFunc(clump_obj)
    def validate(self, clump_obj): print("generic passed")

GENERIC = _GenericShapeType()

class _LinearShapeFunc(BaseShapeFunc):
    pass

class _LinearShapeType(BaseShapeType):
    def accept(self, clump_obj): return _LinearShapeFunc(clump_obj)
    def validate(self, clump_obj):
        if any(s.is_circular for s in clump_obj.sequences):
            raise ShapeError("all sequences must be linear")
        print("linear passed")

LINEAR = _LinearShapeType(GENERIC)

class _DoubleStrandedShapeFunc(BaseShapeFunc):
    def sequence_lengths(self):
        return tuple(len(seq) for seq in self.c.sequences)

class _DoubleStrandedShapeType(BaseShapeType):
    def accept(self, clump_obj): return _DoubleStrandedShapeFunc(clump_obj)
    def validate(self, clump_obj):
        if not len(clump_obj.sequences) == 2:
            raise shape.ShapeError("must have exactly two sequences")
        if not clump_obj.annealments:
            raise shape.ShapeError("must have at least one annealment")
        if not any(all(seq in ann.sequences for seq in clump_obj.sequences)
                   for ann in clump_obj.annealments):
            raise shape.ShapeError("the two sequences must be annealed to "
                                   "each other")
        print("double passed")

DOUBLE_STRANDED = _DoubleStrandedShapeType(GENERIC)

class _PCRTemplateShapeFunc(BaseShapeFunc):
    pass

class _PCRTemplateShapeType(BaseShapeType):
    def accept(self, clump_obj): return _PCRTemplateShapeFunc(clump_obj)
    def validate(self, clump_obj): print("pcr template passed")

PCR_TEMPLATE = _PCRTemplateShapeType(LINEAR, DOUBLE_STRANDED)

def create_pcr_template(nucleotides):
    seq = sequence.LinearSequence(nucleotides)
    rseq = seq.reverse_complement()
    new = clump.Clump()
    new = new.add_sequence(seq)
    new = new.add_sequence(rseq)
    new = new.add_annealment((seq, rseq), (0, 0), len(seq))
    return new
