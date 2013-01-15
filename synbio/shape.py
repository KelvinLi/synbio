"""Shape layer -- framework for extracting features out of clumps"""

from .clump import clump

def declare_op(shape_obj, op_name):
    def out(func):
        shape_obj.declare(op_name, func)
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

    def is_specialization(self, shape_obj):
        if self is shape_obj:
            return True
        return any(g.is_specialization(shape_obj)
                   for g in self._generalizations)

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
    def __init__(self, shape_obj, *args, **kwargs):
        self._shape_obj = shape_obj
        self._clump_obj = clump.Clump()
        self.operate(OpNames.INIT, *args, **kwargs)

    def cast(self, shape_obj=None):
        if shape_obj is None:
            shape_obj = self._shape_obj
        shape_obj.validate(self._clump_obj)
        self._shape_obj = shape_obj

    def is_shape(self, shape_obj):
        return self._shape_obj.is_specialization(shape_obj)

    def operate(self, op_name, *args, **kwargs):
        return self._shape_obj.call_operator(op_name,
                                             self._clump_obj,
                                             *args, **kwargs)
