import sys
import numpy as np
import time
import matplotlib.pyplot as plt


def get_indices(ci, cj, r, N):
    ci0, ci1 = ci - r, ci + r
    cj0, cj1 = cj - r, cj + r
    
    I = (np.arange(ci0, ci1 + 1) + N) % N
    J = (np.arange(cj0, cj1 + 1) + N) % N

    return I, J
    

def compute_integrals_vectorized(n):
    A = np.random.randn(n, n)
    r = int(np.sqrt(n))
    
    runtimes = {
        'get_indices' : [],
        'slice' : [],
        'sum' : [],
        'loop' : []
    }
    
    I = np.array(range(n))
    J = np.array(range(n))
    
    II, JJ = np.meshgrid(I, J)

    return runtimes


if __name__ == '__main__':
    NMAX = 11
    plot_type = 'loglog'
    nsigma = 3
    
    if len(sys.argv) > 1:
        NMAX = int(sys.argv[1])

    if len(sys.argv) > 2:
        plot_type = int(sys.argv[2])
        
    if len(sys.argv) > 3:
        nsigma = int(sys.argv[3])
    
    N = np.array([2**n for n in range(NMAX)])
    
    runtimes = {}
    
    for n in N:
        print('n = {0}'.format(n))
        A = np.random.randn(n, n)
        r = int(np.sqrt(n))
        
        runtimes[n] = {
            'get_indices' : [],
            'slice' : [],
            'sum' : [],
            'loop' : []
        }
        
        loopstart = time.time()
        for ci in range(n):
            for cj in range(n):
                start = time.time()
                I, J = get_indices(ci, cj, r, n)
                finish = time.time()
                runtimes[n]['get_indices'].append(finish - start)
                
                start = time.time()
                F = A[I, J]
                finish = time.time()
                runtimes[n]['slice'].append(finish - start)

                start = time.time()
                np.sum(F[F > 0])
                finish = time.time()
                runtimes[n]['sum'].append(finish - start)
        loopend = time.time()
        runtimes[n]['loop'].append(loopend - loopstart)
        

    runtimes['get_indices'] = {}
    runtimes['get_indices']['min'] = np.array([np.min(runtimes[n]['get_indices']) for n in N])
    runtimes['get_indices']['max'] = np.array([np.max(runtimes[n]['get_indices']) for n in N])
    runtimes['get_indices']['mean'] = np.array([np.mean(runtimes[n]['get_indices']) for n in N])
    runtimes['get_indices']['std'] = np.array([np.std(runtimes[n]['get_indices']) for n in N])

    runtimes['slice'] = {}
    runtimes['slice']['min'] = np.array([np.min(runtimes[n]['slice']) for n in N])
    runtimes['slice']['max'] = np.array([np.max(runtimes[n]['slice']) for n in N])
    runtimes['slice']['mean'] = np.array([np.mean(runtimes[n]['slice']) for n in N])
    runtimes['slice']['std'] = np.array([np.std(runtimes[n]['slice']) for n in N])

    runtimes['sum'] = {}
    runtimes['sum']['min'] = np.array([np.min(runtimes[n]['sum']) for n in N])
    runtimes['sum']['max'] = np.array([np.max(runtimes[n]['sum']) for n in N])
    runtimes['sum']['mean'] = np.array([np.mean(runtimes[n]['sum']) for n in N])
    runtimes['sum']['std'] = np.array([np.std(runtimes[n]['sum']) for n in N])

    fig, ((ax0, ax1), (ax2, ax3)) = plt.subplots(2, 2)
    
    if plot_type == 'cartesian':
        ax0.plot(N, [runtimes[n]['loop'] for n in N], color='green', marker='o')
    elif plot_type == 'semilogy':
        ax0.semilogy(N, [runtimes[n]['loop'] for n in N], color='green', marker='o')
    elif plot_type == 'semilogx':
        ax0.semilogx(N, [runtimes[n]['loop'] for n in N], color='green', marker='o')
    else:
        ax0.loglog(N, [runtimes[n]['loop'] for n in N], color='green', marker='o')
    
    ax0.set_xlabel('N')
    ax0.set_ylabel('Runtime (seconds)')
    ax0.set_title('Full Loop Runtime')
    ax0.grid(True)
    
    
    rmin = runtimes['get_indices']['min']; rmin[rmin <= 0] = np.nan
    rmean = runtimes['get_indices']['mean']; rmean[rmean <= 0] = np.nan
    rmax = runtimes['get_indices']['max']; rmax[rmax <= 0] = np.nan
    rstd = runtimes['get_indices']['std']; rstd[rstd <= 0] = np.nan

    if plot_type == 'cartesian':
        ax1.plot(N, rmin, color='blue', marker='x', label='min')
        ax1.plot(N, rmean, color='green', marker='o', label='mean')
        ax1.plot(N, rmax, color='red', marker='^', label='max')
    elif plot_type == 'semilogy':
        ax1.semilogy(N, rmin, color='blue', marker='x', label='min')
        ax1.semilogy(N, rmean, color='green', marker='o', label='mean')
        ax1.semilogy(N, rmax, color='red', marker='^', label='max')
    elif plot_type == 'semilogx':
        ax1.loglog(N, rmin, color='blue', marker='x', label='min')
        ax1.loglog(N, rmean, color='green', marker='o', label='mean')
        ax1.loglog(N, rmax, color='red', marker='^', label='max')
    else:
        ax1.semilogx(N, rmin, color='blue', marker='x', label='min')
        ax1.semilogx(N, rmean, color='green', marker='o', label='mean')
        ax1.semilogx(N, rmax, color='red', marker='^', label='max')
        
    alpha = 0.5
    beta = nsigma / (nsigma + 2)
    for k in range(0, nsigma):
        label = '$\mu + ' + str(k) + '\sigma$'
        ax1.fill_between(N, rmean + k*rstd, rmean + (k+1)*rstd, facecolor='green', alpha=alpha * beta**k, label=label)
    
    ax1.legend()
    ax1.set_xlabel('N')
    ax1.set_ylabel('Runtime (seconds)')
    ax1.set_title('get_indices() Runtime')
    
    rmin = runtimes['slice']['min']; rmin[rmin <= 0] = np.nan
    rmean = runtimes['slice']['mean']; rmean[rmean <= 0] = np.nan
    rmax = runtimes['slice']['max']; rmax[rmax <= 0] = np.nan
    rstd = runtimes['get_indices']['std']; rstd[rstd <= 0] = np.nan
    ax2.loglog(N, rmin, color='blue', marker='x', label='min')
    ax2.loglog(N, rmean, color='green', marker='o', label='mean')
    ax2.loglog(N, rmax, color='red', marker='^', label='max')
    
    for k in range(0, nsigma):
        label = '$\mu + ' + str(k) + '\sigma$'
        ax2.fill_between(N, rmean + k*rstd, rmean + (k+1)*rstd, facecolor='green', alpha=alpha * beta**k, label=label)
    
    ax2.legend()
    ax2.set_xlabel('N')
    ax2.set_ylabel('Runtime (seconds)')
    ax2.set_title('slice Runtime')

    rmin = runtimes['sum']['min']; rmin[rmin <= 0] = np.nan
    rmean = runtimes['sum']['mean']; rmean[rmean <= 0] = np.nan
    rmax = runtimes['sum']['max']; rmax[rmax <= 0] = np.nan
    rstd = runtimes['get_indices']['std']; rstd[rstd <= 0] = np.nan
    ax3.loglog(N, rmin, color='blue', marker='x', label='min')
    ax3.loglog(N, rmean, color='green', marker='o', label='mean')
    ax3.loglog(N, rmax, color='red', marker='^', label='max')


    for k in range(0, nsigma):
        label = '$\mu + ' + str(k) + '\sigma$'
        ax3.fill_between(N, rmean + k*rstd, rmean + (k+1)*rstd, facecolor='green', alpha=alpha * beta**k, label=label)
        ax3.loglog(N, rmean + rstd, color='black', linestyle=':', label='$\mu + \sigma$')
    
    ax3.legend()
    ax3.set_xlabel('N')
    ax3.set_ylabel('Runtime (seconds)')
    ax3.set_title('np.sum() Runtime')
    
    fig.savefig('timing_{0:s}_nsigma{1:d}.png'.format(plot_type, nsigma))
    plt.show()
    