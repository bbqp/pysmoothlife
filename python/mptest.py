import multiprocessing
import itertools
import numpy as np
import time
import matplotlib.pyplot as plt

def compute_sum_entry(i, j):
    return i + j


def compute_prd_entry(i, j):
    return i * j


if __name__ == '__main__':
    N = 2**np.arange(3, 12)
    num_processes = list(range(0, multiprocessing.cpu_count() + 1))
    runtimes = np.zeros((len(num_processes), len(N)))
    
    for nidx, n in enumerate(N):
        A = np.zeros((n,n))
        B = np.zeros((n,n))
        
        print('Gathering runtimes for initializing {0:4d}x{1:4d} arrays.'.format(n, n))
    
        for pidx, nprocs in enumerate(num_processes):
            if nprocs > 0:
                print('\tGathering runtimes using {0:2d} process(es).'.format(nprocs))
                tstart = time.time()
                with multiprocessing.Pool(processes=nprocs) as p:
                    resultsA = p.starmap_async(compute_sum_entry, itertools.product(range(n), range(n)))
                    resultsB = p.starmap_async(compute_prd_entry, itertools.product(range(n), range(n)))
                    
                    A[:,:] = np.reshape(resultsA.get(), A.shape)
                    B[:,:] = np.reshape(resultsB.get(), B.shape)
                tend = time.time()
            else:
                print('\tGathering runtimes using double for loop.')
                tstart = time.time()
                for i in range(n):
                    for j in range(n):
                        A[i, j] = compute_sum_entry(i, j)
                        B[i, j] = compute_prd_entry(i, j)
                tend = time.time()
            runtimes[pidx, nidx] = tend - tstart
            
            A[:,:] = 0
            B[:,:] = 0
    
    fig, ax = plt.subplots(1, 2)
    fig.set_tight_layout(True)
    colors = ['black', 'purple', 'blue', 'green', 'yellow', 'orange', 'red']
    titles = ['Runtimes for compute_sum_entry()', 'Runtimes for compute_prd_entry()']

    for i, t in enumerate(titles):
        for pidx, nprocs in enumerate(num_processes):
            ax[i].loglog(N, runtimes[pidx,:], marker='x', color=colors[pidx], linestyle='none', label=f'Processes: {nprocs}')
        ax[i].set_xlabel('Dimension')
        ax[i].set_ylabel('Runtime (seconds)')
        ax[i].set_title(t)
        ax[i].legend()
        ax[i].grid(True)
    plt.show()
        