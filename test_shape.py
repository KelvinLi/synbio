from test_mod import *
from synbio import clump, sequence, shapes

test_shape = shapes.LINEAR_DOUBLE \
    .accept(clump.Clump()) \
    .set_top_middle(sequence.LinearSequence(map(make_nucleotide, "attacg"))) \
    .set_overhangs(top5=sequence.LinearSequence(map(make_nucleotide, "ccgg")))
test_clump = test_shape.clump()
dump_clump(test_clump)
dump_clump(shapes.LINEAR_DOUBLE.accept(test_clump).clump())
