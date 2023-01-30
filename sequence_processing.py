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
