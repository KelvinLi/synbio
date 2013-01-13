import dna

def dump_cluster(cluster):
    print(cluster._sequences)
    print(cluster._annealments)
    print(cluster._sequences[0].dump())
    print(cluster._sequences[1].dump())
    print([ann.starts for ann in cluster._annealments])
    print([ann._length for ann in cluster._annealments])
    print("=========")

cluster = dna.Cluster()

cluster.add_sequence(dna.LinearSequence(dna.Nucleotide(n) for n in "agctg"))
cluster.add_sequence(dna.CircularSequence(dna.Nucleotide(n) for n in "gctca"))

cluster.add_annealment(tuple(cluster._sequences), (0, 2), 2)
dump_cluster(cluster)

cluster.add_annealment(tuple(cluster._sequences), (2, 4), 3)
dump_cluster(cluster)

cluster.remove_sequence(cluster._sequences[0])
print(cluster._sequences)
print(cluster._annealments)
print("=========")

cluster.remove_sequence(cluster._sequences[0])
print(cluster._sequences)
print(cluster._annealments)
print("=========")
