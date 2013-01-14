from . import shape

class GenericShape(shape.BaseShape):
    def _validate(self):
        pass # allow all clumps

    @shape.shape_method
    def count_sequences(self):
        return len(self._clump.sequences)
