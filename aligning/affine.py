import numpy as np


def affine_gap_align(seq1, seq2, matrix, gap_score=-7, gap_continue_score=-2):
    gap_score += gap_continue_score

    m, n = len(seq1), len(seq2)
    dp_cc = -np.inf * np.ones((m + 1, n + 1))
    dp_gc = -np.inf * np.ones((m + 1, n + 1))
    dp_cg = -np.inf * np.ones((m + 1, n + 1))
    dp_cc[0, 0] = 0

    for i in range(m + 1):
        for j in range(n + 1):
            if i > 0:
                dp_cg[i, j] = max(dp_cg[i, j], dp_cc[i - 1, j] + gap_score)
                dp_cg[i, j] = max(dp_cg[i, j], dp_cg[i - 1, j] + gap_continue_score)
                dp_cg[i, j] = max(dp_cg[i, j], dp_gc[i - 1, j] + gap_score)
            if j > 0:
                dp_gc[i, j] = max(dp_gc[i, j], dp_cc[i, j - 1] + gap_score)
                dp_gc[i, j] = max(dp_gc[i, j], dp_cg[i, j - 1] + gap_score)
                dp_gc[i, j] = max(dp_gc[i, j], dp_gc[i, j - 1] + gap_continue_score)
            if i > 0 and j > 0:
                cur_score = matrix[(seq1[i - 1], seq2[j - 1])]
                dp_cc[i, j] = max(dp_cc[i, j], dp_cc[i - 1, j - 1] + cur_score)
                dp_cc[i, j] = max(dp_cc[i, j], dp_cg[i - 1, j - 1] + cur_score)
                dp_cc[i, j] = max(dp_cc[i, j], dp_gc[i - 1, j - 1] + cur_score)

    best_score = max(dp_cc[m, n], dp_cg[m, n], dp_gc[m, n])

    align1, align2 = [], []
    i, j = m, n
    if best_score == dp_cc[m, n]:
        level = "cc"
    elif best_score == dp_cg[m, n]:
        level = "cg"
    else:
        level = "gc"

    while i + j > 0:
        if level == "cc":
            cur_score = matrix[(seq1[i - 1], seq2[j - 1])]
            if dp_cc[i, j] == dp_cg[i - 1, j - 1] + cur_score:
                level = "cg"
            elif dp_cc[i, j] == dp_gc[i - 1, j - 1] + cur_score:
                level = "gc"
            i -= 1
            j -= 1
            align1.append(seq1[i])
            align2.append(seq2[j])
        elif level == "cg":
            if dp_cg[i, j] == dp_cc[i - 1, j] + gap_score:
                level = "cc"
            elif dp_cg[i, j] == dp_gc[i - 1, j] + gap_score:
                level = "gc"
            i -= 1
            align1.append(seq1[i])
            align2.append('_')
        else:
            if dp_gc[i, j] == dp_cc[i, j - 1] + gap_score:
                level = "cc"
            elif dp_gc[i, j] == dp_cg[i, j - 1] + gap_score:
                level = "cg"
            j -= 1
            align1.append('_')
            align2.append(seq2[j])

    return (''.join(align1[::-1]), ''.join(align2[::-1])), best_score
