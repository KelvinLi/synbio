"""Shape layer -- framework for safely manipulating clumps"""

class ShapeError(Exception):
    pass

class BaseShape:
    def __init__(*args, **kwargs):
        raise TypeError

class BaseShapeFactory:
    def __init__(self, *dependencies):
        assert not issubclass(BaseShapeFactory, type(self))
        self._dependencies = dependencies

    def accept(self, clump_obj):
        if not clump_obj.sequences():
            return self.analyze_empty()
        # TODO: optimize for diamond relationships
        dep_shape_objs = tuple(d.accept(clump_obj) for d in self._dependencies)
        return self.analyze(clump_obj, dep_shape_objs)
