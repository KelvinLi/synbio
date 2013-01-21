from .. import clump, sequence
from ..shape import BaseShape, BaseShapeFactory, ShapeError

class _Overhangs:
    def __init__(self, *, top5, top3, bottom5, bottom3):
        if any(oh.is_circular for oh in (top5, top3, bottom5, bottom3)):
            raise ShapeError
        self.top5 = top5
        self.top3 = top3
        self.bottom5 = bottom5
        self.bottom3 = bottom3

class _Shape(BaseShape):
    def __init__(self):
        empty = sequence.LinearSequence(tuple())
        self._top_middle = empty
        self._overhangs = _Overhangs(top5=empty, top3=empty,
                                     bottom5=empty, bottom3=empty)

    def clump(self):
        top = self._top_middle.join((self._overhangs.top5,
                                    self._overhangs.top3))
        bottom = self._top_middle.reverse_complement().join(
                 (self._overhangs.bottom5, self._overhangs.bottom3))
        out = clump.Clump()
        out = out.add_sequence(top).add_sequence(bottom)
        assert len(out.sequences()) == 2
        assert 0 in out.locate_sequences(top)
        assert 1 in out.locate_sequences(bottom)
        return out.add_annealment(
               (0, 1),
               (len(self._overhangs.top5), len(self._overhangs.bottom5)),
               len(self._top_middle))

    def set_top_middle(self, replacement):
        if replacement.is_circular:
            raise ShapeError
        new = _Shape()
        new._top_middle = replacement
        new._overhangs = self._overhangs
        return new

    def set_overhangs(self, *, top5=None, top3=None,
                      bottom5=None, bottom3=None):
        if top5 is None:
            top5 = self._overhangs.top5
        if top3 is None:
            top3 = self._overhangs.top3
        if bottom5 is None:
            bottom5 = self._overhangs.bottom5
        if bottom3 is None:
            bottom3 = self._overhangs.bottom3
        new = _Shape()
        new._top_middle = self._top_middle
        new._overhangs = _Overhangs(top5=top5, top3=top3,
                                    bottom5=bottom5, bottom3=bottom3)
        return new

class _ShapeFactory(BaseShapeFactory):
    """Linear, double stranded, with four overhangs.
    5' -------------- 3'  top strand
           |||||||
    3' -------------- 5'  bottom strand
           ^^^^^^^
           "middle" part of strands
    """
    def analyze_empty(self):
        return _Shape()

    # TODO: ugly
    def analyze(self, clump_obj, dep_shape_objs):
        if len(clump_obj.sequences()) != 2:
            raise ShapeError
        if any(seq.is_circular for seq in clump_obj.sequences()):
            raise ShapeError
        top, bottom = clump_obj.sequences()

        results = list(clump_obj.query_annealments(0))
        if len(results) != 1 or results[0].other_key != 1:
            raise ShapeError
        top_start = results[0].start
        top_length = results[0].length
        top5 = top.fragment(0, top_start)
        top3_start = top_start + top_length
        top3 = top.fragment(top3_start, len(top) - top3_start)

        results = list(clump_obj.query_annealments(1))
        if len(results) != 1 or results[0].other_key != 0:
            raise ShapeError
        bottom5 = bottom.fragment(0, results[0].start)
        bottom3_start = results[0].start + results[0].length
        bottom3 = bottom.fragment(bottom3_start, len(bottom) - bottom3_start)

        return _Shape().set_top_middle(top.fragment(top_start, top_length)) \
                       .set_overhangs(top5=top5, top3=top3,
                                      bottom5=bottom5, bottom3=bottom3)

LINEAR_DOUBLE = _ShapeFactory()
