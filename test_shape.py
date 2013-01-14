from synbio import shapes
from synbio.clump import clump

test_clump = clump.Clump()
test_clump.add_sequence(clump.LinearSequence(clump.Nucleotide(n)
                                             for n in "attcg"))
generic_shape = shapes.GenericShape(test_clump)

print(generic_shape.count_sequences())
