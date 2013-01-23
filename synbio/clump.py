"""Raw model of a "clump" of DNA."""

import itertools

# TODO: this annealment merge code is an unholy mess
def _sorted_groupby(iterable, key):
    groups = itertools.groupby(sorted(iterable, key=key), key)
    return (list(g) for k, g in groups)

def _group_annealments_by_keys(annealments):
    groups0 = _sorted_groupby(annealments, lambda a: a.keys[0])
    for group0 in groups0:
        groups1 = _sorted_groupby(group0, lambda a: a.keys[1])
        for group1 in groups1:
            yield group1

def _merge_annealment_pair(annA, annB):
    assert annA.keys == annB.keys
    if tuple(s + annA.length for s in annA.starts) != annB.starts:
        return None
    return _Annealment(annA.keys, annA.starts, annA.length + annB.length)

def _merge_annealment_group_once(group):
    modified = False
    i = 0
    while i + 1 < len(group):
        new = _merge_annealment_pair(group[i], group[i + 1])
        if new is None:
            i += 1
            continue
        group.pop(i)
        group.pop(i)
        group.insert(0, new)
        modified = True
        i += 1
    return modified

def _merge_annealment_group(group):
    out = sorted(group, key=lambda a: a.starts[0])
    assert len(out) >= 1
    if len(out) <= 1:
        return out
    while _merge_annealment_group_once(out):
        out = sorted(out, key=lambda a: a.starts[0])
    return out

def _has_overlap(common_seq, starts, lengths):
    diff = starts[0] - starts[1]
    if not common_seq.is_circular:
        return -lengths[0] < diff < lengths[1]
    L = len(common_seq)
    assert all(0 < x <= L for x in lengths)
    return diff % L < lengths[1] or diff % L > L - lengths[0]

class _Annealment:
    def __init__(self, keys, starts, length):
        self.keys = keys
        self.starts = starts
        self.length = length

class AnnealmentQueryResult:
    def __init__(self, start, length, other_key):
        self.start = start
        self.length = length
        self.other_key = other_key

class Clump:
    def __init__(self):
        self._sequences = tuple()
        self._annealments = tuple()

    def __append_sequence(self, sequence):
        new = Clump()
        new._sequences = self._sequences + (sequence,)
        new._annealments = self._annealments
        return new

    def __replace_annealments(self, replacement):
        new = Clump()
        new._sequences = self._sequences
        new._annealments = replacement
        return new

    def _validate_annealment(self, annealment):
        f0, f1 = (self._sequences[key].fragment(start, annealment.length)
                  for key, start in zip(annealment.keys, annealment.starts))
        if not f0.reverse_complement().has_same_nucleotides(f1):
            raise ValueError("sequences must be reverse complementary over "
                             "annealment region")

    def _has_annealment_conflict(self, annealment):
        return any(_has_overlap(self._sequences[key],
                                (r.start, start),
                                (r.length, annealment.length))
                   for key, start in zip(annealment.keys, annealment.starts)
                   for r in self.query_annealments(key))

    def _merge_annealments(self):
        groups = _group_annealments_by_keys(self._annealments)
        merged_groups = (_merge_annealment_group(group) for group in groups)
        replacement = tuple(annealment
                            for merged_group in merged_groups
                            for annealment in merged_group)
        return self.__replace_annealments(replacement)

    def add_annealment(self, keys, starts, length, *, overwrite=False):
        keys = tuple(keys)
        starts = tuple(starts)
        if not len(keys) == len(starts) == 2:
            raise ValueError("must have two keys and two start points")
        if not length > 0:
            raise ValueError("length must be positive")
        if not all(isinstance(s, int) for s in starts):
            raise TypeError("start points must be integers")
        if not all(s >= 0 for s in starts):
            raise ValueError("start points must be non-negative")
        if not all(isinstance(k, int) for k in keys):
            raise TypeError("keys must be integers")
        if not all(0 <= k < len(self._sequences) for k in keys):
            raise ValueError("invalid keys")
        if keys[0] > keys[1]:
            keys = tuple(reversed(keys))
            starts = tuple(reversed(starts))
        new = _Annealment(keys, starts, length)
        if overwrite:
            raise NotImplementedError
        if self._has_annealment_conflict(new):
            raise ValueError("refusing to overwrite conflicting annealment(s)")
        self._validate_annealment(new)
        return self.__replace_annealments(self._annealments + (new,)) \
                   ._merge_annealments()

    def add_sequence(self, new):
        return self.__append_sequence(new)

    def sequences(self):
        return tuple(self._sequences)

    def locate_sequences(self, match):
        return (key for key, seq in enumerate(self._sequences)
                if seq.has_same_nucleotides(match))

    # TODO: this function is nested too deeply
    def query_annealments(self, key):
        for ann in self._annealments:
            for i, o in (0, 1), (1, 0):
                if ann.keys[i] != key:
                    continue
                a = (ann.starts[i], ann.length, ann.keys[o])
                yield AnnealmentQueryResult(*a)
