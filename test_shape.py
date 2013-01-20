from test_mod import *
from synbio import sequence, shapes

test_clump = shapes.create_pcr_template(map(make_nucleotide, "attacg"))
shapes.PCR_TEMPLATE.examine(test_clump)
print(shapes.DOUBLE_STRANDED.examine(test_clump).sequence_lengths())
print(shapes.GENERIC.examine(test_clump).count_sequences())
