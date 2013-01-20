from synbio import clump, sequence

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
