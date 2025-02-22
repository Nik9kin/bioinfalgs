import numpy as np


def viterbi(pi, A, B, D):
    N, M, T = B.shape[0], B.shape[1], len(D)
    pi, A, B = np.log(pi), np.log(A), np.log(B)

    # delta[t, i] = max_{q_0, ..., q_t-1} log P(q_0, ..., q_t-1, q_t = i, d_0, ..., d_t | pi, A, B)
    delta = np.zeros((T, N))
    psi = np.zeros((T, N), dtype=int)

    delta[0, :] = pi + B[:, D[0]]
    for t in range(1, T):
        delta[t, :] = np.max(delta[t - 1, :] + A.T + B[:, D[t], np.newaxis], axis=1)
        psi[t, :] = np.argmax(delta[t - 1, :] + A.T + B[:, D[t], np.newaxis], axis=1)

    mle_states = [np.argmax(delta[T - 1, :])]
    for t in range(T - 1, 0, -1):
        mle_states.append(psi[t, mle_states[-1]])
    mle_states = mle_states[::-1]

    return mle_states, delta

