from scipy.stats import entropy


class MotifMatrix:
    def __init__(self, motifs: list[str]):
        self._height = len(motifs)
        self._length = len(motifs[0])
        for seq in motifs:
            if len(seq) != self._length:
                raise ValueError('motif strings must be the same length')
        self._motifs = [seq.upper() for seq in motifs]
        self._profile = None
        self._consensus = None
        self._compute_profile()
        self._compute_consensus()

    def _compute_profile(self):
        self._profile = [{c: 0 for c in 'ATGC'} for _ in range(self._length)]
        for seq in self._motifs:
            for pos in range(self._length):
                self._profile[pos][seq[pos]] += 1

    def _compute_consensus(self):
        res = []
        for pos in range(self._length):
            res.append(max(self._profile[pos], key=self._profile[pos].get))
        self._consensus = ''.join(res)

    def extend(self, motif: str):
        if len(motif) != self._length:
            raise ValueError('the motif string must be the same length as the motif matrix')
        self._height += 1
        self._motifs.append(motif.upper())
        new_consensus = []
        for pos in range(self._length):
            profile_col = self._profile[pos]
            profile_col[self._motifs[-1][pos]] += 1
            if profile_col[self._motifs[-1][pos]] > profile_col[self._consensus[pos]]:
                new_consensus.append(self._motifs[-1][pos])
            else:
                new_consensus.append(self._consensus[pos])
        self._consensus = ''.join(new_consensus)

    def likelihood(self, seq: str, smooth=True) -> float:
        if self._length != len(seq):
            return 0.0
        res = 1.0
        for pos in range(self._length):
            if smooth:
                res *= self._profile[pos][seq[pos]] + 1
            else:
                res *= self._profile[pos][seq[pos]]
        if smooth:
            res /= (self._height + 4) ** self._length
        else:
            res /= self._height ** self._length
        return res

    def max_likelihood(self, seq: str, smooth=True) -> (str, int, float):
        best_likelihood = 0.0
        best_start_pos = 0
        for start_pos in range(len(seq) - self._length + 1):
            likelihood = self.likelihood(seq[start_pos: start_pos + self._length], smooth)
            if likelihood > best_likelihood:
                best_likelihood = likelihood
                best_start_pos = start_pos
        return seq[best_start_pos: best_start_pos + self._length], best_start_pos, best_likelihood

    def score(self, mod='count') -> int | float:
        res = 0
        if mod == 'count':
            for pos in range(self._length):
                res += self._height - self._profile[pos][self._consensus[pos]]
        elif mod == 'entropy':
            for pos in range(self._length):
                res += entropy(list(self._profile[pos].values()), base=2)
        else:
            raise ValueError(f"mod must be one of 'count' or 'entropy' (got '{mod}')")
        return res

    def get_motifs(self):
        return self._motifs

    def get_profile(self):
        return self._profile

    def get_consensus(self):
        return self._consensus


def greedy_motif_search(seqs: list[str], length: int, smooth=True) -> MotifMatrix:
    best_motifs = MotifMatrix([seq[:length] for seq in seqs])
    for pos in range(len(seqs[0]) - length + 1):
        cur_motifs = MotifMatrix([seqs[0][pos: pos + length]])
        for i in range(1, len(seqs)):
            motif, *_ = cur_motifs.max_likelihood(seqs[i], smooth)
            cur_motifs.extend(motif)
        if cur_motifs.score() < best_motifs.score():
            best_motifs = cur_motifs
    return best_motifs


if __name__ == '__main__':
    motifs = [
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
    m = MotifMatrix(motifs)
    print(m.get_consensus())
    print(m.get_profile())
    print(m.score(mod='entropy'))
    print(m.likelihood('TCGTGGATTTCC'))
    print(m.likelihood(m.get_consensus()))
