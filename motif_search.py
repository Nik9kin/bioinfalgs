from scipy.stats import entropy


class MotifMatrix:
    def __init__(self, motifs: list[str]):
        self._height = len(motifs)
        self._length = len(motifs[0])
        for seq in motifs:
            if len(seq) != self._length:
                raise ValueError('motif strings must be the same _length')
        self._motifs = [seq.upper() for seq in motifs]
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

    def score(self, mod='count'):
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

    def get_profile(self):
        return self._profile

    def get_consensus(self):
        return self._consensus


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
