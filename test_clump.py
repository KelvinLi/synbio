from test_mod import *
from synbio import clump
from synbio import sequence

test_clump = clump.Clump()

test_clump = test_clump.add_sequence(
             sequence.LinearSequence(make_nucleotide(n) for n in "agctg"))
test_clump = test_clump.add_sequence(
             sequence.CircularSequence(make_nucleotide(n) for n in "gctca"))

test_clump = test_clump.add_annealment((0, 1), (0, 2), 2)
dump_clump(test_clump)

test_clump = test_clump.add_annealment((0, 1), (2, 4), 3)
dump_clump(test_clump)
