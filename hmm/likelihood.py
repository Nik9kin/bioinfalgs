import numpy as np


def forward_backward(pi, A, B, D):
    N, M, T = B.shape[0], B.shape[1], len(D)
    pi, A, B = np.log(pi), np.log(A), np.log(B)

    # alpha[t, i] = log P(d_0, d_1, ..., d_t, q_t = i | pi, A, B)
    # beta[t, i]  = log P(d_t+1, d_t+2, ..., d_T-1 | q_t = i, pi, A, B)
    alpha = np.zeros((T, N))
    beta = np.zeros((T, N))

    # forward pass
    alpha[0, :] = pi + B[:, D[0]]
    for t in range(1, T):
        alpha[t, :] = np.logaddexp.reduce(alpha[t - 1, :] + A.T + B[:, D[t], np.newaxis], axis=1)

    likelihood_1 = np.logaddexp.reduce(alpha[T - 1, :])

    # backward pass
    for t in range(T - 2, -1, -1):
        beta[t, :] = np.logaddexp.reduce(beta[t + 1, :] + A + B[:, D[t + 1]], axis=1)

    likelihood_2 = np.logaddexp.reduce(pi + B[:, D[0]] + beta[0, :])
    assert abs(likelihood_1 - likelihood_2) < 1e-8, (likelihood_1, likelihood_2)

    # gamma[t, i] = log P(q_t = i | D, pi, A, B)
    gamma = alpha + beta
    gamma -= likelihood_1

    return gamma, likelihood_1
