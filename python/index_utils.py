import numpy as np


def compute_index_ranges(cspan, centeri, centerj):
    I = np.array([-cspan + centeri + k for k in range(2*cspan + 1)])
    J = np.array([-cspan + centerj + k for k in range(2*cspan + 1)])

    return I, J


def index_range_mod(I, J, M, N):
    return np.mod(I, M), np.mod(J, N)


def compute_k_indices(cspan, centeri, centerj, M, N):
    I, J = compute_index_ranges(cspan, centeri, centerj)
    IM, JM = index_range_mod(I, J, M, N)
    
    II = np.tile(IM, reps=2*cspan+1)
    JJ = np.repeat(JM, axis=0, repeats=2*cspan+1)
    KK = II + JJ * M
    
    return II, JJ, KK


if __name__ == '__main__':
    cspan = 2
    M = 64
    N = 64
    centeri, centerj = 2, 0
    
    I, J = compute_index_ranges(cspan, centeri, centerj)
    I = np.tile(I, reps=2*cspan+1)
    J = np.repeat(J, axis=0, repeats=2*cspan+1)
    
    II, JJ, KK = compute_k_indices(cspan, centeri, centerj, M, N)
    
    for idx, (i, j, im, jm, km) in enumerate(zip(I, J, II, JJ, KK)):
        print('{0:4d}: {1:4d}, {2:4d} -> {3:4d}, {4:4d} -> {5:6d}'.format(idx, i, j, im, jm, km))