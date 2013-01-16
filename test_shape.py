from synbio import sequence, shape, shapes

def make_nucleotide(s):
    if s == "a":
        return sequence.Nucleotide.A
    elif s == "c":
        return sequence.Nucleotide.C
    elif s == "g":
        return sequence.Nucleotide.G
    elif s == "t":
        return sequence.Nucleotide.T
    else:
        raise ValueError

test_instance = shape.ShapeInstance(shapes.PCR_TEMPLATE, map(make_nucleotide, "attacg"))
test_instance.cast()

test_instance.cast(shapes.DOUBLE_STRANDED)
print(test_instance.operate("sequence_lengths"))

test_instance.cast(shapes.GENERIC)
print(test_instance.operate("count_sequences"))
