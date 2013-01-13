from synbio import views
from synbio.clump import clump

test_clump = clump.Clump()
test_clump.add_sequence(clump.LinearSequence(clump.Nucleotide(n)
                                             for n in "attcg"))
generic_view = views.GenericView(test_clump)

print(generic_view.count_sequences())
