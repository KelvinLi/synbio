from .. import clump, sequence
from ..shape import BaseShape, BaseShapeFactory, ShapeError

def _get_overhangs(top, bottom, ann):
    top3_start = ann.starts[0] + ann.length
    top3_length = len(top) - top3_start
    bottom3_start = ann.starts[1] + ann.length
    bottom3_length = len(bottom) - bottom3_start
    return _Overhangs(top5=top.fragment(0, ann.starts[0]),
                      bottom5=bottom.fragment(0, ann.starts[1]),
                      top3=top.fragment(top3_start, top3_length),
                      bottom3=bottom.fragment(bottom3_start, bottom3_length))

class _Overhangs:
    def __init__(self, *, top5, top3, bottom5, bottom3):
        if any(oh.is_circular for oh in (top5, top3, bottom5, bottom3)):
            raise ShapeError
        self.top5 = top5
        self.top3 = top3
        self.bottom5 = bottom5
        self.bottom3 = bottom3

class _Shape(BaseShape):
    def __init__(self, top_middle, overhangs):
        self._top_middle = top_middle
        self._overhangs = overhangs

    def clump(self):
        top = self._top_middle.join((self._overhangs.top5,
                                    self._overhangs.top3))
        bottom = self._top_middle.reverse_complement().join(
                 (self._overhangs.bottom5, self._overhangs.bottom3))
        out = clump.Clump()
        out = out.add_sequence(top).add_sequence(bottom)
        return out.add_annealment(
               (top, bottom),
               (len(self._overhangs.top5), len(self._overhangs.bottom5)),
               len(self._top_middle))

    def set_top_middle(self, replacement):
        if replacement.is_circular:
            raise ShapeError
        return _Shape(replacement, self._overhangs)

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
        return _Shape(self._top_middle,
                      _Overhangs(top5=top5, top3=top3,
                                 bottom5=bottom5, bottom3=bottom3))

class _ShapeFactory(BaseShapeFactory):
    """Linear, double stranded, with four overhangs.
    5' -------------- 3'  top strand
           |||||||
    3' -------------- 5'  bottom strand
           ^^^^^^^
           "middle" part of strands
    """
    def analyze_empty(self):
        empty = sequence.LinearSequence(tuple())
        return _Shape(empty, _Overhangs(top5=empty, top3=empty,
                                        bottom5=empty, bottom3=empty))

    def analyze(self, clump_obj, dep_shape_objs):
        if not len(clump_obj.sequences) == 2:
            raise ShapeError
        if any(seq.is_circular for seq in clump_obj.sequences):
            raise ShapeError
        top, bottom = clump_obj.sequences
        query = clump.AnnealmentQuery((top, bottom), (None, None), None)
        results = list(clump_obj.query_annealments(query))
        if len(results) != 1:
            raise ShapeError
        ann = results[0]
        return _Shape(top.fragment(ann.starts[0], ann.length),
                      _get_overhangs(top, bottom, ann))

LINEAR_DOUBLE = _ShapeFactory()
