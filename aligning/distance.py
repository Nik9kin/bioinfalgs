def hamming(seq1, seq2):
    if len(seq1) != len(seq2):
        raise ValueError('Input objects must be the same length')

    return sum(c1 != c2 for c1, c2 in zip(seq1, seq2))


def find_substring(pattern, text):
    pattern_len = len(pattern)
    text_len = len(text)
    if pattern_len > text_len:
        raise ValueError('The pattern length is longer than the text length')

    pos, dist_best = 0, hamming(pattern, text[:pattern_len])
    for i in range(1, text_len - pattern_len + 1):
        dist_cur = hamming(pattern, text[i: i + pattern_len])
        if dist_cur < dist_best:
            pos, dist_best = i, dist_cur
    return pos, text[pos: pos + pattern_len], dist_best


def levenshtein(seq1, seq2):
    n, m = len(seq1), len(seq2)
    if n < m:
        seq1, seq2 = seq2, seq1
        n, m = m, n

    dp_prev, dp_cur = (list(range(m + 1)) for _ in range(2))
    for c1 in seq1:
        dp_cur[0] = dp_prev[0] + 1
        for j, c2 in enumerate(seq2):
            dp_cur[j + 1] = min(dp_cur[j] + 1, dp_prev[j + 1] + 1, dp_prev[j] + (c1 != c2))
        dp_prev, dp_cur = dp_cur, dp_prev
    return dp_prev[-1]


def test_distances(seq1, seq2):
    print(f'Sequence 1: {seq1}')
    print(f'Sequence 2: {seq2}')
    try:
        print(f'Hamming distance: {hamming(seq1, seq2)}')
    except ValueError:
        print('The Hamming distance cannot be calculated: the sequences have different lengths')
    print(f'Levenshtein distance: {levenshtein(seq1, seq2)}')


if __name__ == '__main__':
    from Bio import SeqIO
    gattaca_seqs = [getattr(rec, 'seq') for rec in SeqIO.parse('data/gattaca.fasta', 'fasta')]
    f8_seqs = [getattr(rec, 'seq') for rec in SeqIO.parse('data/f8.fasta', 'fasta')]
    test_distances(*gattaca_seqs)
    test_distances(*f8_seqs)
    for pattern in gattaca_seqs:
        for text in f8_seqs:
            print(find_substring(pattern, text))
