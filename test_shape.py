from synbio import shape, shapes

test_instance = shape.ShapeInstance(shapes.PCR_TEMPLATE, "attacg")
test_instance.cast()

test_instance.cast(shapes.DOUBLE_STRANDED)
print(test_instance.operate("sequence_lengths"))

test_instance.cast(shapes.GENERIC)
print(test_instance.operate("count_sequences"))
