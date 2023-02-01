import copy
import random
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

    def __eq__(self, other):
        return self._motifs == other.get_motifs()

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

    def extend(self, motif: str, index: int = None):
        if len(motif) != self._length:
            raise ValueError('the motif string must be the same length as the motif matrix')
        motif = motif.upper()
        if index is None:
            index = self._height
        self._motifs.insert(index, motif)
        self._height += 1
        for pos in range(self._length):
            self._profile[pos][motif[pos]] += 1
        self._compute_consensus()

    def reduce(self, index: int):
        motif = self._motifs.pop(index)
        self._height -= 1
        for pos in range(self._length):
            self._profile[pos][motif[pos]] -= 1
        self._compute_consensus()

    def likelihood(self, seq: str, smooth=True):
        if len(seq) < self._length:
            return 0.0
        likelihoods = []
        for start_pos in range(len(seq) - self._length + 1):
            res = 1.0
            for pos in range(self._length):
                if smooth:
                    res *= self._profile[pos][seq[start_pos + pos]] + 1
                else:
                    res *= self._profile[pos][seq[start_pos + pos]]
            if smooth:
                res /= (self._height + 4) ** self._length
            else:
                res /= self._height ** self._length
            likelihoods.append(res)
        return likelihoods if len(likelihoods) > 1 else likelihoods[0]

    def most_likelihood(self, seq: str, smooth=True) -> (str, int, float):
        if len(seq) <= self._length:
            return seq, 0, float(self.likelihood(seq, smooth=smooth))
        likelihoods = self.likelihood(seq, smooth=smooth)
        best_pos = int(max(range(len(likelihoods)), key=lambda i: likelihoods[i]))
        return seq[best_pos: best_pos + self._length], best_pos, likelihoods[best_pos]

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


def greedy_motif_search(seqs: list[str], length: int, smooth=True, scoring='entropy') -> MotifMatrix:
    best_motifs = MotifMatrix([seq[:length] for seq in seqs])
    for pos in range(len(seqs[0]) - length + 1):
        cur_motifs = MotifMatrix([seqs[0][pos: pos + length]])
        for i in range(1, len(seqs)):
            motif, *_ = cur_motifs.most_likelihood(seqs[i], smooth)
            cur_motifs.extend(motif)
        if cur_motifs.score(mod=scoring) < best_motifs.score(mod=scoring):
            best_motifs = cur_motifs
    return best_motifs


def randomized_motif_search(seqs: list[str],
                            length: int,
                            smooth=True,
                            scoring='entropy',
                            max_steps=200) -> (MotifMatrix, int | float):
    init_pos = [random.randrange(len(seq) - length + 1) for seq in seqs]
    M = MotifMatrix([seq[pos: pos + length] for seq, pos in zip(seqs, init_pos)])
    best_score = float('inf')
    best_M = M
    for _ in range(max_steps):
        new_M = MotifMatrix([M.most_likelihood(seq, smooth)[0] for seq in seqs])
        new_score = new_M.score(mod=scoring)
        if new_score < best_score:
            best_score = new_score
            best_M = M = new_M
        elif new_score == best_score and M == new_M:
            break
        else:
            M = new_M
    return best_M, best_score


def gibbs_motif_search(seqs: list[str],
                       length: int,
                       smooth=True,
                       scoring='entropy',
                       max_steps=2000) -> (MotifMatrix, int | float):
    seqs_num = len(seqs)
    init_pos = [random.randrange(len(seq) - length + 1) for seq in seqs]
    M = MotifMatrix([seq[pos: pos + length] for seq, pos in zip(seqs, init_pos)])
    best_score = float('inf')
    best_M = copy.deepcopy(M)
    for _ in range(max_steps):
        row_for_change = random.randrange(seqs_num)
        M.reduce(row_for_change)
        probs = M.likelihood(seqs[row_for_change], smooth=smooth)
        new_pos = random.choices(range(len(probs)), weights=probs, k=1)[0]
        new_motif = seqs[row_for_change][new_pos: new_pos + length]
        M.extend(new_motif, row_for_change)
        new_score = M.score(mod=scoring)
        if new_score < best_score:
            best_score = new_score
            best_M = copy.deepcopy(M)
    return best_M, best_score


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
