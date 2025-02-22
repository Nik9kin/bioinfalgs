import numpy as np


def matrix_align(seq1, seq2, matrix, gap_score=0):
    m, n = len(seq1), len(seq2)
    dp = -np.inf * np.ones((m + 1, n + 1))
    dp[0, 0] = 0

    for i in range(m + 1):
        for j in range(n + 1):
            if i > 0:
                dp[i, j] = max(dp[i, j], dp[i - 1, j] + gap_score)
            if j > 0:
                dp[i, j] = max(dp[i, j], dp[i, j - 1] + gap_score)
            if i > 0 and j > 0:
                dp[i, j] = max(dp[i, j], dp[i - 1, j - 1] + matrix[(seq1[i - 1], seq2[j - 1])])

    align1, align2 = [], []
    i, j = m, n
    while i + j > 0:
        if i > 0 and dp[i, j] == dp[i - 1, j] + gap_score:
            i -= 1
            align1.append(seq1[i])
            align2.append('_')
        elif j > 0 and dp[i, j] == dp[i, j - 1] + gap_score:
            j -= 1
            align1.append('_')
            align2.append(seq2[j])
        else:
            i -= 1
            j -= 1
            align1.append(seq1[i])
            align2.append(seq2[j])

    return (''.join(align1[::-1]), ''.join(align2[::-1])), dp[m][n]
