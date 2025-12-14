import numpy as np
import time
import matplotlib.pyplot as plt

def run_benchmarking(x0=-1, xn=1, y0=-1, yn=1, npts=None):
    if npts is None:
        npts = 2**np.array(list(range(8)))+1
    
    runtimes = []
    mgridtimes = []
    dsttimes = []

    for k in range(len(npts)):
        N = npts[k]
        
        x = np.linspace(x0, xn, N)
        y = np.linspace(y0, yn, N)
        
        tstart = time.time()
        XX, YY = np.meshgrid(x, y)
        M, N = XX.shape
        XM, YM = XX[:M-1,:N-1], YY[:M-1,:N-1]
        XMT = np.tile(XM, (M*N,1,1))
        YMT = np.tile(YM, (M*N,1,1))
        D = np.zeros(XMT.shape)
        CX = np.tile(XX.flatten()[:, np.newaxis, np.newaxis], (1, M-1, N-1))
        CY = np.tile(YY.flatten()[:, np.newaxis, np.newaxis], (1, M-1, N-1))
        tend = time.time()
        mgridtimes.append(tend - tstart)
        
        print('XX.shape =', XX.shape)
        print('XM.shape =', XM.shape)
        print('XMT.shape =', XMT.shape)
        print('CX.shape =', CX.shape)
        
        print('Running benchmarking for grid of size {0:3d}x{1:3d}'.format(N, N))
        tstart = time.time()
        D, dstart, dend = compute_distances_tensor_prealloc(XMT, YMT, CX, CY, D)
        tend = time.time()
        runtimes.append(tend - tstart)
        dsttimes.append(dend - dstart)

    fig, ax = plt.subplots(1, 3)
    fig.set_tight_layout(True)
    
    ax[0].plot(npts, runtimes, marker='x', linestyle='none', label='Distance Computation Runtimes')
    fitcurve = np.poly1d(np.polyfit(npts, runtimes, 3))
    ax[0].plot(npts, fitcurve(npts), marker='none', linestyle='--', label='Curve of Best Fit')
    ax[0].set_xlabel('Dimension')
    ax[0].set_ylabel('Runtime (seconds)')
    ax[0].grid(True)
    ax[0].legend()
    
    ax[1].plot(npts, mgridtimes, marker='x', linestyle='none', label='Grid Allocation Runtimes')
    # fitcurve = np.poly1d(np.polyfit(npts, runtimes, 3))
    # ax[1].plot(npts, fitcurve(npts), marker='none', linestyle='--', label='Curve of Best Fit')
    ax[1].set_xlabel('Dimension')
    ax[1].set_ylabel('Runtime (seconds)')
    ax[1].grid(True)
    ax[1].legend()

    ax[2].plot(npts, dsttimes, marker='x', linestyle='none', label='Distance Array Operation Runtimes')
    fitcurve = np.poly1d(np.polyfit(npts, dsttimes, 3))
    ax[2].plot(npts, fitcurve(npts), marker='none', linestyle='--', label='Curve of Best Fit')
    ax[2].set_xlabel('Dimension')
    ax[2].set_ylabel('Runtime (seconds)')
    ax[2].grid(True)
    ax[2].legend()
    
    plt.show()
        

def compute_distances_loop(XX, YY):
    M, N = XX.shape
    
    # Compute the midpoints.
    XM, YM = XX[:M-1,:N-1], YY[:M-1,:N-1]
    
    for i in range(M):
        for j in range(N):
            c = XX[i, j], YY[i, j]
            
            

def compute_distances_tensor(XX, YY):
    M, N = XX.shape
    
    # Compute the midpoints.
    XM, YM = XX[:M-1,:N-1], YY[:M-1,:N-1]

    CX = np.tile(XX.flatten()[np.newaxis, np.newaxis, :], (M-1, N-1,1))
    CY = np.tile(YY.flatten()[np.newaxis, np.newaxis, :], (M-1, N-1,1))
    
    dstart = time.time()
    DX = XM[:,:,np.newaxis] - CX
    DY = YM[:,:,np.newaxis] - CY
    D = np.sqrt(DX**2 + DY**2)
    dend = time.time()
    
    amin = np.argmin(D, axis=2)
    amax = np.argmax(D, axis=2)
    
    return D, dstart, dend


def compute_distances_tensor_prealloc(XMT, YMT, CX, CY, D):
    dstart = time.time()
    D[:,:,:] = np.sqrt((XMT - CX)**2 + (YMT - CY)**2)
    dend = time.time()
    
    amin = np.argmin(D, axis=2)
    amax = np.argmax(D, axis=2)
    
    return D, dstart, dend


if __name__ == '__main__':
    npts = np.array([10*k for k in range(1,11)]) + 1
    run_benchmarking(npts=npts)