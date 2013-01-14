class ShapeError(Exception):
    pass

class BaseShape:
    def __init__(self, clump):
        self._clump = clump
        self._validate()

def shape_method(func):
    def out(self_shape, *args, **kwargs):
        self_shape._validate()
        return func(self_shape, *args, **kwargs)

    return out
