"""Shape layer -- framework for extracting features out of clumps"""

from .clump import clump

def declare_op(shape_type, name):
    def out(func):
        shape_type.declare_operator(name, func)
        return func
    return out

class OpNames:
    VALIDATE = "_validate"
    INIT = "_init"

class ShapeError(Exception):
    pass

class ShapeType:
    def __init__(self, *generalizations):
        self._generalizations = tuple(generalizations)
        self._operators = dict()

    def __eq__(self, other):
        return self is other

    def is_specialization(self, shape_type):
        if self is shape_type:
            return True
        return any(g.is_specialization(shape_type)
                   for g in self._generalizations)

    def validate(self, clump_obj):
        # TODO: optimize for diamond relationships
        for g in self._generalizations:
            g.validate(clump_obj)
        self.call_operator(OpNames.VALIDATE, clump_obj, tuple(), dict())

    def call_operator(self, name, clump_obj, args, kwargs):
        return self._operators[name](clump_obj, *args, **kwargs)

    def declare_operator(self, name, func):
        if name in self._operators:
            raise AttributeError
        self._operators[name] = func

class ShapeInstance:
    def __init__(self, shape_type, *args, **kwargs):
        self._shape_type = shape_type
        self._clump_obj = clump.Clump()
        self.operate(OpNames.INIT, *args, **kwargs)

    def cast(self, new=None):
        if new is None:
            self._shape_type.validate(self._clump_obj)
            return
        if not self._shape_type.is_specialization(new):
            new.validate(self._clump_obj)
        self._shape_type = new

    def is_shape(self, shape_type):
        return self._shape_type.is_specialization(shape_type)

    def operate(self, name, *args, **kwargs):
        return self._shape_type.call_operator(
               name, self._clump_obj, args, kwargs)
