from smoothlifegrid import *
import time
import matplotlib
import matplotlib.pyplot as plt


def benchmark(numcells=[16, 32, 64], ntrials=10, nprocesses=list(range(multiprocessing.cpu_count()+1)), show=True):
    numcells = np.array(numcells)
    nprocesses = np.array(nprocesses)

    # An array for the loop runtimes.
    ltimes = np.zeros((ntrials, len(numcells), len(nprocesses)))

    for pidx, nprocs in enumerate(nprocesses):
        serial = (nprocs == 0)
        
        if serial or nprocs == 1:
            runstring = 'double for loop' if serial else 'multiprocessing pool'
            print('Running serial execution with {0:2d} process(es) and {1:s}'.format(nprocs, runstring))
        else:
            runstring = 'multiprocessing pool'
            print('Running serial execution with {0:2d} process(es) and {1:s}'.format(nprocs, runstring))

        for cidx, ncells in enumerate(numcells):
            print('\tRunning with grid size {0:4d}x{1:4d}'.format(ncells, ncells))

            # Create the grid and set the state.
            slg = SmoothLifeGrid(numcellsw=ncells, numcellsh=ncells)
            slg.set_state_function(np.random.rand(ncells, ncells))

            for trial in range(ntrials):
                print('\t\tTrial {0:3d}'.format(trial))
                lstart = time.time()
                slg.compute_integrals(serial=serial, nprocs=nprocs)
                lend = time.time()
                ltimes[trial, cidx, pidx] = lend - lstart

    ltimes_mean = np.mean(ltimes, axis=0)
    ltimes_std = np.std(ltimes, axis=0)
    ltimes_max = np.max(ltimes, axis=0)

    fig, ax = plt.subplots(nrows=1, ncols=1)
    fig.set_tight_layout(True)

    labels = []
    nsigma = 3
    for k in range(nsigma + 1):
        if nsigma > 0:
            labels.append('$\\mu + {0:d}\\sigma$'.format('' if nsigma == 1 else nsigma))
        else:
            labels.append('$\\mu$ (Mean Runtime)')

    labels.append('Maximum Runtime')

    colors = matplotlib.colormaps['rainbow'](np.linspace(0, 1, len(nprocesses)))

    for pidx, nprocs in enumerate(nprocesses):
        label = 'Serial' if nprocs == 0 else 'Total Processes: {0}'.format(nprocs)
        ax.loglog(numcells, ltimes_mean[:, pidx], color=colors[pidx,:], marker='x', linestyle='none', label=label)

    ax.set_title('Quadrature Loop Runtimes')
    ax.set_xlabel('Dimension')
    ax.set_ylabel('Runtime (seconds)')
    ax.grid(True)
    ax.legend()

    if show:
        plt.show()

    fig.savefig('runtimes.png')

def benchmark_integrals_and_transition(numcells=32, ri_ratio=0.05, b=0.05/4, dtype=np.float32):
    slg = SmoothLifeGrid(numcellsw=numcells, numcellsh=numcells, ri_ratio=ri_ratio, b=b, dtype=dtype)
    slg.compute_integrals()
    
def draw(numcells=32, ri_ratio=0.05, b=0.05/4, dtype=np.float32):
    slg = SmoothLifeGrid(numcellsw=numcells, numcellsh=numcells, ri_ratio=ri_ratio, b=b, dtype=dtype)
    slg.set_state_function(1)
    slg.draw_circular_region()
    slg.draw_annular_region()

if __name__ == '__main__':
    numcells=[32, 64, 128, 256]
    ri_ratios=[0.1, 0.05, 0.025, 0.0125]
    ri_ratios=[ra for ra in ri_ratios]
    b = [ra for ra in ri_ratios]
    #benchmark(numcells=numcells, nprocesses=[0])
    benchmark_integrals_and_transition(numcells=numcells[-2], ri_ratio=ri_ratios[-2], b=b[-2], dtype=np.float64)
    #for k in range(len(numcells)):
    #    draw(numcells=numcells[k], ri_ratio=ri_ratios[k], b=b[k], dtype=np.float16)