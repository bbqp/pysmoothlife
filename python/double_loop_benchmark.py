import numpy as np
import matplotlib.pyplot as plt
import time


unitmap = {
    9 : 'G',
    6 : 'M',
    3 : 'K',
    0 : '',
    -3 : 'm',
    -6 : 'u',
    -9 : 'n',
    -12 : 'p',
    -15 : 'f',
    -18 : 'a'
}

def get_mantissa_and_units(x):
    logx = np.log10(x)

    # Round the exponent to the lowest power of 3 in magnitude.
    if logx < 0:
        exponent = int(np.ceil(logx))
    else:
        exponent = int(np.floor(logx))

    # Now find the nearest power of 3 to the exponent.
    exponent = (exponent // 3) * 3
        
    mantissa = 10**(logx - exponent)
    
    return mantissa, unitmap[exponent] + 's'

def generate_random_subgrid(dim, dtype):
    return np.random.rand(dim, dim).astype(dtype)


def array_sum_test(p0=-1, pn=1, ri_ratio=[0.2, 0.1, 0.05, 0.025], rfactor=3, numcells=2**np.arange(4, 8), dtype=np.float64):
    if (ri_ratio is None) or (len(ri_ratio) != len(numcells)):
        initial_ratio = numcells[0] / 80 if ri_ratio is None else ri_ratio[0]
        ri_ratio = np.array([initial_ratio / 2**k for k in range(len(numcells))])

    ltimes = np.zeros(numcells.shape)
    sumtimes = []
    for k, (ri_rat, nc) in enumerate(zip(ri_ratio, numcells)):
        ri = (pn - p0) * ri_rat;
        ro = rfactor * ri;
        dp = (pn - p0) / nc
        cspani, cspano = int(np.ceil(ri / dp)), int(np.ceil(ro / dp))

        CI = np.zeros((nc, nc), dtype=dtype)
        AI = np.zeros((nc, nc), dtype=dtype)

        print('--- Computing sums for {0:4d}x{0:4d} grid with 2*cspani+1 = {1:2d} and 2*cspano+1 = {2:2d}'.format(nc, 2*cspani+1, 2*cspano+1))
        lstart = time.time()
        for i in range(nc - 1):
            for j in range(nc - 1):
                sstart = time.time()
                cint = np.sum(generate_random_subgrid(2*cspani+1, dtype=dtype))
                #aint = np.sum(generate_random_subgrid(2*cspano+1, dtype=dtype))
                send = time.time()
                
                sumtimes.append(send - sstart)
                
                CI[i, j] = cint
                #AI[i, j] = aint
        lend = time.time()
        ltimes[k] = lend - lstart

    fig, ax = plt.subplots(1)
    ax.loglog(numcells, ltimes, marker='x', linestyle='none')
    ax.set_xlabel('Dimension')
    ax.set_ylabel('Runtime (seconds)')
    
    meantime = np.mean(sumtimes)
    mantissa, units = get_mantissa_and_units(meantime)
    ax.set_title('Double Loop Runtime\n' + 'np.sum mean runtime $\\approx$ {0:.2f} {1:s}'.format(mantissa, units))
    ax.grid(True)
    plt.show()


if __name__ == '__main__':
    array_sum_test(numcells=2**np.arange(4, 8), dtype=np.float32)