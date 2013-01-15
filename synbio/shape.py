"""Shape layer -- framework for extracting features out of clumps"""

from .clump import clump

def declare_op(shape, op_name):
    def out(func):
        shape.declare(op_name, func)
        return func
    return out

class OpNames:
    VALIDATE = "_validate"
    INIT = "_init"

class ShapeError(Exception):
    pass

class Shape:
    def __init__(self, generalizations=tuple()):
        self._generalizations = tuple(generalizations)
        self._ops = dict()

    def __eq__(self, other):
        return self is other

    def is_specialization(self, shape):
        if self is shape:
            return True
        return any(g.is_specialization(shape) for g in self._generalizations)

    def validate(self, clump_obj):
        # TODO: optimize for diamond relationships
        for g in self._generalizations:
            g.validate(clump_obj)
        self.call_operator(OpNames.VALIDATE, clump_obj)

    def call_operator(self, op_name, *args, **kwargs):
        return self._ops[op_name](*args, **kwargs)

    def declare(self, op_name, func):
        if op_name in self._ops:
            raise AttributeError
        self._ops[op_name] = func

class ShapeInstance:
    def __init__(self, shape, *args, **kwargs):
        self._shape = shape
        self._clump_obj = clump.Clump()
        self.operate(OpNames.INIT, *args, **kwargs)

    def cast(self, shape=None):
        if shape is None:
            shape = self._shape
        shape.validate(self._clump_obj)
        self._shape = shape

    def is_shape(self, shape):
        return self._shape.is_specialization(shape)

    def operate(self, op_name, *args, **kwargs):
        return self._shape.call_operator(op_name,
                                         self._clump_obj,
                                         *args, **kwargs)
