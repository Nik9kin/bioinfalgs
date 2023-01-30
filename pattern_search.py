from sequence_processing import complement


def hamming_distance(s1: str, s2: str) -> int:
    if len(s1) != len(s2):
        raise ValueError("input strings must be the same length")
    res = 0
    for i in range(len(s1)):
        if s1[i] != s2[i]:
            res += 1
    return res


def neighbours(dna: str, d: int) -> dict[str: int]:
    if d == 0:
        return {dna: 0}
    if d < 0:
        raise ValueError("d must be non-negative")
    if len(dna) == 1:
        return {"A": 1, "T": 1, "G": 1, "C": 1, dna: 0}

    res = dict()
    suffix_neighbours = neighbours(dna[1:], d)
    diff_nucleotides = set("ATGC").difference(dna[0])
    for neighbour, dist in suffix_neighbours.items():
        res[dna[0] + neighbour] = dist
        if dist < d:
            for nucleotide in diff_nucleotides:
                res[nucleotide + neighbour] = dist + 1
    return res


def pattern_matching(text: str, pattern: str, mismatches: int = 0) -> list[int]:
    res = []
    pattern_len = len(pattern)
    for i in range(len(text) - pattern_len + 1):
        if hamming_distance(text[i:i + pattern_len], pattern) <= mismatches:
            res.append(i)
    return res


def frequent_k_mers(text: str, k: int, freq_threshold=0, mismatches=0,
                    reverse_complementary=False) -> list[tuple[str, int]]:
    freq_dict = {}
    highest_frequency = 0
    for i in range(len(text) - k + 1):
        k_mer = text[i:i + k]
        for neighbor in neighbours(k_mer, mismatches):
            freq_dict[neighbor] = freq_dict.setdefault(neighbor, 0) + 1
            highest_frequency = max(highest_frequency, freq_dict[neighbor])
            if reverse_complementary:
                rc_neighbor = complement(neighbor)
                freq_dict[rc_neighbor] = freq_dict.setdefault(rc_neighbor, 0) + 1
                highest_frequency = max(highest_frequency, freq_dict[rc_neighbor])
    res = []
    for k_mer, freq in freq_dict.items():
        if freq_threshold:
            if freq >= freq_threshold:
                res.append((-freq, k_mer))
        else:
            if freq == highest_frequency:
                res.append((-freq, k_mer))
    res.sort()
    return [(k_mer, -freq) for (freq, k_mer) in res]


def find_clumps(text: str, k: int,
                window_length: int, freq_threshold: int) -> list[str]:
    res = set()
    for i in range(len(text) - window_length + 1):
        window = text[i:i + window_length]
        for k_mer, freq in frequent_k_mers(window, k, freq_threshold):
            res.add(k_mer)
    return sorted(list(res))
