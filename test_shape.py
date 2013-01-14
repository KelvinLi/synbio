from synbio import shape, shapes

test_instance = shape.ShapeInstance(shapes.PCR_TEMPLATE, "attacg")
test_instance.cast()
test_instance.cast(shapes.DOUBLE_STRANDED)
