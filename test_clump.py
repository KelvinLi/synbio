from synbio import clump
from synbio import sequence

def dump_nuc(n):
    if n == sequence.Nucleotide.A:
        return "A"
    elif n == sequence.Nucleotide.C:
        return "C"
    elif n == sequence.Nucleotide.G:
        return "G"
    elif n == sequence.Nucleotide.T:
        return "T"
    else:
        raise ValueError

def dump_seq(s):
    return "seq: " + "".join(dump_nuc(n) for n in s.nucleotides)

def dump_clump(clump):
    print("CLUMP:")
    print("   ", clump.sequences)
    print("   ", clump.annealments)
    for s in clump.sequences:
        print("   ", dump_seq(s))
    print("   ", [ann.starts for ann in clump.annealments], [ann.length for ann in clump.annealments])

def dump_query(query):
    print("QUERY:")
    print("   ", query.sequences)
    print("   ", query.starts, query.length)

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
