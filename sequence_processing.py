def reverse_complement(seq: str) -> str:
    seq = seq.replace("A", "$").replace("T", "A").replace("$", "T")
    seq = seq.replace("C", "$").replace("G", "C").replace("$", "G")
    return seq[::-1]


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
