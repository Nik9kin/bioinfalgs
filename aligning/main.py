import numpy as np
from Bio.Align import substitution_matrices as matlist

from affine import affine_gap_align
from matrix import matrix_align


matrix = matlist.load('BLOSUM62')

with open("output.txt", 'w') as f_out:
    seq1 = "DVKVDDRQHGRINCPCNSRPKPPLVLLPKWQAKGLFRPFPDPNHRPKDWSFGCFEFIRFRRWNRHTDYAIGSNLMHSYYIHMAWI"
    seq2 = "DVKVDDRQHGRINCAEYHTFCNSRPKPPLVLLPKWQAFLSLFRPFPWSFGCFEFIRFRRWNGSYYIHMAMI"
    res, score = affin_gap_align(seq1, seq2, matrix=matrix, gap_score=-11, gap_continue_score=-1)
    f_out.write(str(score))
    f_out.write('\n')
    f_out.write(res[0])
    f_out.write('\n')
    f_out.write(res[1])
