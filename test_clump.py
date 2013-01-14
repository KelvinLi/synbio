from synbio.clump import clump

def dump_clump(clump):
    print(clump.sequences)
    print(clump.annealments)
    print(clump.sequences[0].dump())
    print(clump.sequences[1].dump())
    print([ann.starts for ann in clump.annealments])
    print([ann.length for ann in clump.annealments])
    print("=========")

test_clump = clump.Clump()

test_clump.add_sequence(clump.LinearSequence(clump.Nucleotide(n) for n in "agctg"))
test_clump.add_sequence(clump.CircularSequence(clump.Nucleotide(n) for n in "gctca"))

test_clump.add_annealment(tuple(test_clump.sequences), (0, 2), 2)
dump_clump(test_clump)

test_clump.add_annealment(tuple(test_clump.sequences), (2, 4), 3)
dump_clump(test_clump)

test_clump.remove_sequence(test_clump.sequences[0])
print(test_clump.sequences)
print(test_clump.annealments)
print("=========")

test_clump.remove_sequence(test_clump.sequences[0])
print(test_clump.sequences)
print(test_clump.annealments)
print("=========")
