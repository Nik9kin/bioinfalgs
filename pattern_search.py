import sequence_processing as sp


def neighbours(seq: str, radius: int) -> dict[str: int]:
    if radius == 0:
        return {seq: 0}
    if radius < 0:
        raise ValueError('radius must be non-negative')
    if len(seq) == 1:
        return {'A': 1, 'T': 1, 'G': 1, 'C': 1, seq: 0}

    res = dict()
    suffix_neighbours = neighbours(seq[1:], radius)
    diff_nucleotides = set('ATGC').difference(seq[0])
    for neighbour, dist in suffix_neighbours.items():
        res[seq[0] + neighbour] = dist
        if dist < radius:
            for nucleotide in diff_nucleotides:
                res[nucleotide + neighbour] = dist + 1
    return res


def pattern_matching(text: str, pattern: str, mismatches: int = 0) -> list[int]:
    res = []
    pattern_len = len(pattern)
    for i in range(len(text) - pattern_len + 1):
        if sp.hamming_distance(text[i: i + pattern_len], pattern) <= mismatches:
            res.append(i)
    return res


def frequent_k_mers(text: str, k: int, freq_threshold=0, mismatches=0,
                    reverse_complementary=False) -> list[tuple[str, int]]:
    freq_dict = {}
    highest_frequency = 0
    for i in range(len(text) - k + 1):
        k_mer = text[i: i + k]
        for neighbor in neighbours(k_mer, mismatches):
            freq_dict[neighbor] = freq_dict.setdefault(neighbor, 0) + 1
            highest_frequency = max(highest_frequency, freq_dict[neighbor])
            if reverse_complementary:
                rc_neighbor = sp.complement(neighbor)
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
        window = text[i: i + window_length]
        for k_mer, freq in frequent_k_mers(window, k, freq_threshold):
            res.add(k_mer)
    return sorted(list(res))


def motif_search_bruteforce(seqs: list[str], length: int, mismatches: int) -> list[str]:
    motif_candidates = set()
    seq = seqs[0]
    for i in range(len(seq) - length + 1):
        motif_candidates.update(neighbours(seq[i: i + length], mismatches).keys())
    res = list()
    for motif in motif_candidates:
        motif_flag = True
        for seq in seqs:
            seq_flag = False
            for pos in range(len(seq) - length + 1):
                if sp.hamming_distance(seq[pos: pos + length], motif) <= mismatches:
                    seq_flag = True
                    break
            if not seq_flag:
                motif_flag = False
                break
        if motif_flag:
            res.append(motif)
    return res


def median_string(seqs: list[str], length: int) -> (str, int):
    best_dist = length * len(seqs)
    res = 'A' * length
    for pattern in sp.all_nucleotide_seqs(length):
        dist = sp.distance_between_pattern_and_strings(pattern, seqs)
        if dist < best_dist:
            best_dist = dist
            res = pattern
    return res, best_dist
