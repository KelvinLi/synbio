from synbio import sequence, shapes

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

test_clump = shapes.create_pcr_template(map(make_nucleotide, "attacg"))
shapes.PCR_TEMPLATE.examine(test_clump)
print(shapes.DOUBLE_STRANDED.examine(test_clump).sequence_lengths())
print(shapes.GENERIC.examine(test_clump).count_sequences())
