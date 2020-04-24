import numpy as np

def repeat_center(n, repeat):
    return np.zeros((repeat, n))

def fullfact(levels):
    n = len(levels)  # number of factors
    nb_lines = np.prod(levels)  # number of trial conditions
    H = np.zeros((nb_lines, n))

    level_repeat = 1
    range_repeat = np.prod(levels)
    for i in range(n):
        range_repeat //= levels[i]
        lvl = []
        for j in range(levels[i]):
            lvl += [j] * level_repeat
        rng = lvl * range_repeat
        level_repeat *= levels[i]
        H[:, i] = rng
    return H

def ff2n(n):
    return 2 * fullfact([2] * n) - 1

def bbdesign(n, center=None):
    if n == 1:
        return np.array([[-1.0],[1.0],[0.0]])

    # First, compute a factorial DOE with 2 parameters
    H_fact = ff2n(n)
    # Now we populate the real DOE with this DOE

    # We made a factorial design on each pair of dimensions
    # - So, we created a factorial design with two factors
    # - Make two loops
    Index = 0
    nb_lines = (n * (n - 1) / 2) * H_fact.shape[0]
    H = repeat_center(n, int(nb_lines))

    for i in range(n - 1):
        for j in range(i + 1, n):
            Index = Index + 1
            H[max([0, (Index - 1) * H_fact.shape[0]]):Index * H_fact.shape[0], i] = H_fact[:, 0]
            H[max([0, (Index - 1) * H_fact.shape[0]]):Index * H_fact.shape[0], j] = H_fact[:, 1]

    if center is None:
        if n <= 16:
            points = [0, 0, 0, 3, 3, 6, 6, 6, 8, 9, 10, 12, 12, 13, 14, 15, 16]
            center = points[n]
        else:
            center = n

    H = np.c_[H.T, repeat_center(n, center).T].T
    H = np.unique(H,axis=0)
    return H