import numpy as np


class MsaMatrix:
    def __init__(self, seqs: list[str]):
        self._seqs = [seq.upper() for seq in seqs]
        self._length = len(seqs[0])
        for seq in seqs:
            if len(seq) != self._length:
                raise ValueError("Input strings must be the same length")

        self._profile = None
        self._consensus = None
        self._compute_profile()
        self._compute_consensus()

    def __len__(self):
        return self._length

    def __str__(self):
        res = [f"seq # {i:3} | " + seq for i, seq in enumerate(self._seqs)]
        res.append('-' * (self._length + 12))
        res.append("Consensus | " + self._consensus)
        return '\n'.join(res)

    def _compute_profile(self):
        self._profile = [{c: 0 for c in 'ATGC_'} for _ in range(self._length)]
        for seq in self._seqs:
            for pos, char in enumerate(seq):
                self._profile[pos][char] += 1

    def _compute_consensus(self):
        res = []
        for pos in range(self._length):
            res.append(max(self._profile[pos], key=self._profile[pos].get))
        self._consensus = ''.join(res)

    def extend(self, new_consensus):
        new_seqs = []
        for seq in self._seqs:
            new_seq = []
            pos = 0
            for char in new_consensus:
                if char == '.':
                    new_seq.append('_')
                else:
                    new_seq.append(seq[pos])
                    pos += 1
            new_seqs.append(''.join(new_seq))

        self._seqs = new_seqs
        self._length = len(new_consensus)
        self._compute_profile()
        self._compute_consensus()

    def update(self, other):
        if len(other) != self._length:
            raise ValueError("MSA matrices must be the same length")

        self._seqs.extend(other.get_seqs())
        self._compute_profile()
        self._compute_consensus()

    def get_seqs(self):
        return self._seqs

    def get_consensus(self):
        return self._consensus


def align_pair(
        seq1,
        seq2,
        match_score,
        mismatch_penalty,
        gap_penalty,
        return_alignment=False
):
    m, n = len(seq1), len(seq2)
    dp = np.zeros((m + 1, n + 1))
    dp[:, 0] = gap_penalty * np.arange(m + 1)
    dp[0, :] = gap_penalty * np.arange(n + 1)

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if seq1[i - 1] == '_' or seq2[j - 1] == '_':
                score = gap_penalty
            elif seq1[i - 1] == seq2[j - 1]:
                score = match_score
            else:
                score = mismatch_penalty
            dp[i, j] = max(
                dp[i - 1, j] + gap_penalty,
                dp[i, j - 1] + gap_penalty,
                dp[i - 1, j - 1] + score,
            )

    if return_alignment:
        align1, align2 = [], []
        i, j = m, n
        while i + j > 0:
            if i > 0 and dp[i, j] == dp[i - 1, j] + gap_penalty:
                i -= 1
                align1.append(seq1[i])
                align2.append('.')
            elif j > 0 and dp[i, j] == dp[i, j - 1] + gap_penalty:
                j -= 1
                align1.append('.')
                align2.append(seq2[j])
            else:
                i -= 1
                j -= 1
                align1.append(seq1[i])
                align2.append(seq2[j])

        return dp[m][n], (''.join(align1[::-1]), ''.join(align2[::-1]))

    return dp[m][n]


def greedy_msa(seqs: list[str], match_score=1, mismatch_penalty=-1, gap_penalty=-1):
    n = len(seqs)
    msa_parts = [MsaMatrix([seq]) for seq in seqs]
    similarity = np.zeros((n, n), dtype=np.float64)
    np.fill_diagonal(similarity, -np.inf)
    for i in range(n):
        for j in range(i + 1, n):
            similarity[i, j] = align_pair(
                msa_parts[i].get_consensus(),
                msa_parts[j].get_consensus(),
                match_score=match_score,
                mismatch_penalty=mismatch_penalty,
                gap_penalty=gap_penalty,
            )
            similarity[j, i] = similarity[i, j]

    i = 0
    for _ in range(n - 1):
        i, j = np.unravel_index(similarity.argmax(), similarity.shape)
        _, alignment = align_pair(
            msa_parts[i].get_consensus(),
            msa_parts[j].get_consensus(),
            match_score=match_score,
            mismatch_penalty=mismatch_penalty,
            gap_penalty=gap_penalty,
            return_alignment=True
        )
        msa_parts[i].extend(alignment[0])
        msa_parts[j].extend(alignment[1])
        msa_parts[i].update(msa_parts[j])

        similarity[:, j] = -np.inf
        similarity[j, :] = -np.inf

        for k in range(n):
            if similarity[i, k] != -np.inf:
                similarity[i, k] = align_pair(
                    msa_parts[i].get_consensus(),
                    msa_parts[k].get_consensus(),
                    match_score=match_score,
                    mismatch_penalty=mismatch_penalty,
                    gap_penalty=gap_penalty,
                )
                similarity[k, i] = similarity[i, k]

    return msa_parts[i]


if __name__ == '__main__':
    seqs_1 = ["ACT", "ATC", "GCT", "ATCC"]
    seqs_2 = [
        'TCGGGGGTTTTT',
        'CCGGTGACTTAC',
        'ACGGGGATTTTC',
        'TTGGGGACTTTT',
        'AAGGGGACTTCC',
        'TTGGGGACTTCC',
        'TCGGGGATTCAT',
        'TCGGGGATTCCT',
        'TAGGGGAACTAC',
        'TCGGGTATAACC',
    ]
    seqs_3 = [
        "ttgacagctagctcagtcctaggtataatgctagc",
        "ttgacagctagctcagtcctaggtataatgctagc",
        "tttacagctagctcagtcctaggtattatgctagc",
        "ttgacagctagctcagtcctaggtactgtgctagc",
        "ctgatagctagctcagtcctagggattatgctagc",
        "ttgacagctagctcagtcctaggtattgtgctagc",
        "tttacggctagctcagtcctaggtactatgctagc",
        "tttacggctagctcagtcctaggtatagtgctagc",
        "tttacggctagctcagccctaggtattatgctagc",
        "ctgacagctagctcagtcctaggtataatgctagc",
        "tttacagctagctcagtcctagggactgtgctagc",
        "tttacggctagctcagtcctaggtacaatgctagc",
        "ttgacggctagctcagtcctaggtatagtgctagc",
        "ctgatagctagctcagtcctagggattatgctagc",
        "ctgatggctagctcagtcctagggattatgctagc",
        "tttatggctagctcagtcctaggtacaatgctagc",
        "tttatagctagctcagcccttggtacaatgctagc",
        "ttgacagctagctcagtcctagggactatgctagc",
        "ttgacagctagctcagtcctagggattgtgctagc",
        "ttgacggctagctcagtcctaggtattgtgctagc"
    ]

    msa = greedy_msa(seqs_3)
    print(msa)
