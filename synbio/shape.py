"""Shape layer -- framework for safely manipulating clumps"""

class ShapeError(Exception):
    pass

class BaseShapeFunc:
    def __init__(self, clump_obj):
        assert not issubclass(BaseShapeFunc, type(self))
        self.c = clump_obj

class BaseShapeType:
    def __init__(self, *generalizations):
        assert not issubclass(BaseShapeType, type(self))
        self._generalizations = generalizations

    def examine(self, clump_obj):
        # TODO: optimize for diamond relationships
        for g in self._generalizations:
            g.examine(clump_obj)
        self.validate(clump_obj)
        return self.accept(clump_obj)
