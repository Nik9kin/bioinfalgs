from itertools import product


def all_nucleotide_seqs(length: int) -> str:
    for seq in product('ATGC', repeat=length):
        yield ''.join(seq)


def complement(seq: str, mod='dna', reverse=True) -> str:
    pairs = {'G': 'C', 'C': 'G', 'g': 'c', 'c': 'g', 'N': 'N', 'n': 'n'}
    if mod == 'dna':
        pairs.update({'A': 'T', 'T': 'A', 'a': 't', 't': 'a'})
    elif mod == 'rna':
        pairs.update({'A': 'U', 'U': 'A', 'a': 'u', 'u': 'a'})
    else:
        raise ValueError(f"mod must be one of 'dna' or 'rna' (got '{mod}')")
    res = [pairs[nucleotide] for nucleotide in seq]
    if reverse:
        res = res[::-1]
    return ''.join(res)


def distance_between_pattern_and_strings(pattern: str, strings: list[str]) -> int:
    k = len(pattern)
    res = 0
    for seq in strings:
        best_dist = k
        for pos in range(len(seq) - k + 1):
            dist = hamming_distance(pattern, seq[pos: pos + k])
            if dist < best_dist:
                best_dist = dist
        res += best_dist
    return res


def hamming_distance(seq1: str, seq2: str) -> int:
    if len(seq1) != len(seq2):
        raise ValueError('input strings must be the same length')
    res = 0
    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            res += 1
    return res


def minimum_skew(dna: str) -> list[int]:
    skew = [0]
    cur = 0
    min_skew = 0
    for nucleotide in dna:
        if nucleotide == "G":
            cur += 1
        if nucleotide == "C":
            cur -= 1
        if cur < min_skew:
            min_skew = cur
        skew.append(cur)
    res = []
    for i in range(len(skew)):
        if skew[i] == min_skew:
            res.append(i)
    return res
