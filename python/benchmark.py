def benchmark(numcellsx=[16, 32, 64, 128, 256], numcellsy=None, show=True):
    numcellsx = np.array(numcellsx)

    if numcellsy is None:
        numcellsy = numcellsx

    circle = np.zeros((3, len(numcellsx)))
    annulus = np.zeros((3, len(numcellsx)))
    loop = np.zeros((3, len(numcellsx)))

    for counter, (ncx, ncy) in enumerate(zip(numcellsx, numcellsy)):
        print('-I- Running benchmarking for {0:3d}x{1:3d} grid'.format(ncx, ncy))

        # Create the coordinates describing the grid.
        numx = ncx + 1
        numy = ncy + 1
        dx = (xn - x0) / ncx
        dy = (yn - y0) / ncy

        # Crete arrays and grids of coordinates.
        x = np.linspace(x0, xn, numx)
        y = np.linspace(y0, yn, numy)
        XX, YY = np.meshgrid(x, y)

        # Create the state function.
        M, N = ncx, ncy
        F = np.ones((M, N))

        ctimes = []
        atimes = []
        ltimes = []

        lstart = time.time()
        for i in range(numx):
            for j in range(numy):
                c = XX[i, j], YY[i, j]

                cstart = time.time()
                compute_circular_integral(F=F, c=c, r=ri, b=b, width=numx-1, height=numy-1, x0=x0, xn=xn, y0=y0, yn=yn)
                cend = time.time()
                ctimes.append(cend - cstart)

                astart = time.time()
                compute_annular_integral(F=F, c=c, ri=ri, ro=ro, b=b, width=numx-1, height=numy-1, x0=x0, xn=xn, y0=y0, yn=yn)
                aend = time.time()
                atimes.append(aend - astart)

        lend = time.time()
        ltimes.append(lend - lstart)

        circle[0, counter] = np.min(ctimes)
        circle[1, counter] = np.mean(ctimes)
        circle[2, counter] = np.max(ctimes)

        annulus[0, counter] = np.min(atimes)
        annulus[1, counter] = np.mean(atimes)
        annulus[2, counter] = np.max(atimes)

        loop[0, counter] = np.min(ltimes)
        loop[1, counter] = np.mean(ltimes)
        loop[2, counter] = np.max(ltimes)

        counter += 1

    figwidth = 12
    figheight = figwidth
    fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(figwidth,figwidth))
    fig.set_tight_layout(True)
    labels = ['Minimum Runtime', 'Mean Runtime', 'Maximum Runtime', 'Curve of Best Fit for Mean']
    colors = ['blue', 'green', 'red']
    markers = ['s', 'x', '^']

    for k in range(3):
        ax[0,0].semilogy(numcellsx, circle[k,:], color=colors[k], marker=markers[k], linestyle='none', label=labels[k])

    if len(numcellsx) > 1:
        fitcurve = np.poly1d(np.polyfit(numcellsx, circle[1,:], 2))
        ax[0,0].semilogy(numcellsx, fitcurve(numcellsx), color=colors[1], marker='none', linestyle='--', label=labels[-1])

    ax[0,0].set_title('Circular Quadrature Runtimes')
    ax[0,0].set_xlabel('Dimension')
    ax[0,0].set_ylabel('Runtime (seconds)')
    ax[0,0].grid(True)
    ax[0,0].legend()

    for k in range(3):
        ax[0,1].semilogy(numcellsx, annulus[k,:], color=colors[k], marker=markers[k], linestyle='none', label=labels[k])

    if len(numcellsx) > 1:
        fitcurve = np.poly1d(np.polyfit(numcellsx, annulus[1,:], 2))
        ax[0,1].semilogy(numcellsx, fitcurve(numcellsx), color=colors[1], marker='none', linestyle='--', label=labels[-1])

    ax[0,1].set_title('Annualar Quadrature Runtimes')
    ax[0,1].set_xlabel('Dimension')
    ax[0,1].set_ylabel('Runtime (seconds)')
    ax[0,1].grid(True)
    ax[0,1].legend()

    for k in range(3):
        ax[1,0].semilogy(numcellsx, loop[k,:], color=colors[k], marker=markers[k], linestyle='none', label=labels[k])

    if len(numcellsx) > 1:
        fitcurve = np.poly1d(np.polyfit(numcellsx, loop[1,:], 2))
        ax[1,0].semilogy(numcellsx, fitcurve(numcellsx), color=colors[1], marker='none', linestyle='--', label=labels[-1])

    ax[1,0].set_title('Grid Loop Runtimes')
    ax[1,0].set_xlabel('Dimension')
    ax[1,0].set_ylabel('Runtime (seconds)')
    ax[1,0].grid(True)
    ax[1,0].legend()
    
    ax[1,1].set_visible(False)

    if show:
        plt.show()

    fig.savefig('runtimes.png')
    
if __name__ == '__main__':
    