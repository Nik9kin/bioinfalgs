import random
import warnings


class Node:
    def __init__(self, s1: str, s2: str, score, left=None, right=None):
        self.seqs = (s1, s2)
        self.score = score
        self.left = left
        self.right = right
        self.width = self.calc_width()

    def calc_width(self):
        if self.left is None:
            return sum(map(len, self.seqs)) + 3
        else:
            return self.left.width + self.right.width + 1

    def _print_prepare(self, end='\n'):
        cur_node_string = '(' + self.seqs[0] + ',' + self.seqs[1] + ')'
        if self.left is None:
            return [cur_node_string + end]

        left_strings = self.left._print_prepare(end='')
        right_strings = self.right._print_prepare(end=end)
        if len(left_strings) < len(right_strings):
            left_strings.extend(
                [' ' * self.left.width
                 for _ in range(len(right_strings) - len(left_strings))]
            )
        elif len(left_strings) > len(right_strings):
            right_strings.extend(
                [' ' * self.right.width + end
                 for _ in range(len(left_strings) - len(right_strings))]
            )

        left_indent = (self.width - len(cur_node_string)) // 2
        right_indent = self.width - len(cur_node_string) - left_indent
        res = [' ' * left_indent + cur_node_string + ' ' * right_indent + end,
               ' ' * left_indent + '/' + ' ' * (len(cur_node_string) - 2) + '\\' + ' ' * right_indent + end]
        for s1, s2 in zip(left_strings, right_strings):
            res.append(s1 + ' ' + s2)

        return res

    def __str__(self):
        return ''.join(self._print_prepare())


def hirschberg(s1: str,
               s2: str,
               match_score=2,
               mismatch_penalty=-1,
               del_penalty=-2,
               ins_penalty=-2):
    if mismatch_penalty > match_score:
        raise ValueError("The mismatch_penalty is greater than the match_score.")
    if del_penalty + ins_penalty > match_score:
        raise ValueError("The sum of the ins_penalty and del_penalty is greater "
                         "than the match_score.")
    if del_penalty + ins_penalty > mismatch_penalty:
        msg = ("Nonsensical penalties. The sum of the ins_penalty and "
               "del_penalty is greater than the mismatch_penalty. There will "
               "be no mismatches in alignment.")
        warnings.warn(msg)

    seqs = (s1, s2)
    is_swapped = False
    if len(s1) < len(s2):
        s1, s2 = s2, s1
        is_swapped = True
        del_penalty, ins_penalty = ins_penalty, del_penalty

    if not s2:
        return Node(*seqs, score=del_penalty * len(s1))
    if len(s2) == 1:
        if s2 in s1:
            cur_score = match_score
        else:
            cur_score = max(mismatch_penalty, del_penalty + ins_penalty)
        return Node(*seqs, score=del_penalty * (len(s1) - 1) + cur_score)

    n, m = len(s1), len(s2)
    row = [0] * (m + 1)
    for j in range(m):
        row[j + 1] = row[j] + ins_penalty

    for i in range(n // 2):
        prev_val = (0, row[0])
        row[0] += del_penalty
        for j in range(m):
            prev_val = (prev_val[1], row[j + 1])

            if s1[i] == s2[j]:
                cur_score = match_score
            else:
                cur_score = mismatch_penalty

            row[j + 1] = max(row[j] + ins_penalty,
                             row[j + 1] + del_penalty,
                             prev_val[0] + cur_score)

    par = list(range(m + 1))
    for i in range(n // 2, n):
        prev_val = (0, row[0])
        prev_par = (0, par[0])
        row[0] += del_penalty
        for j in range(m):
            prev_val = (prev_val[1], row[j + 1])
            prev_par = (prev_par[1], par[j + 1])

            if s1[i] == s2[j]:
                cur_score = match_score
            else:
                cur_score = mismatch_penalty

            row[j + 1] += del_penalty
            if row[j] + ins_penalty > row[j + 1]:
                row[j + 1] = row[j] + ins_penalty
                par[j + 1] = par[j]
            if prev_val[0] + cur_score > row[j + 1]:
                row[j + 1] = prev_val[0] + cur_score
                par[j + 1] = prev_par[0]

    seqs_left = [s1[: n // 2], s2[:par[m]]]
    seqs_right = [s1[n // 2:], s2[par[m]:]]
    if is_swapped:
        seqs_left = seqs_left[::-1]
        seqs_right = seqs_right[::-1]
    scoring = (match_score, mismatch_penalty, del_penalty, ins_penalty)
    return Node(*seqs,
                score=row[m],
                left=hirschberg(*seqs_left, *scoring),
                right=hirschberg(*seqs_right, *scoring))


if __name__ == '__main__':
    seq1 = "AGTACGCA"
    seq2 = "TATGC"
    root = hirschberg(seq1, seq2)
    print(root)
