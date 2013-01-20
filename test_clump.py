from test_mod import *
from synbio import clump
from synbio import sequence

test_clump = clump.Clump()

test_clump = test_clump.add_sequence(
             sequence.LinearSequence(make_nucleotide(n) for n in "agctg"))
test_clump = test_clump.add_sequence(
             sequence.CircularSequence(make_nucleotide(n) for n in "gctca"))

test_clump = test_clump.add_annealment(test_clump.sequences, (0, 2), 2)
dump_clump(test_clump)

test_clump = test_clump.add_annealment(tuple(test_clump.sequences), (2, 4), 3)
dump_clump(test_clump)

query_result = list(test_clump.query_annealments(clump.AnnealmentQuery((None, None), (None, None), None)))
for q in query_result:
    dump_query(q)

test_clump = test_clump.remove_sequence(test_clump.sequences[0])
dump_clump(test_clump)

test_clump = test_clump.remove_sequence(test_clump.sequences[0])
dump_clump(test_clump)
