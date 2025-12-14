import time
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

#-------------------------------------------------------
# Methods for finding extreme points on polygonal
# shapes.
#-------------------------------------------------------

def get_extreme_points_on_lines(c, points, where='b'):
    '''
    get_extreme_points_on_line(c, points)

    Get the locations of the minimizer and maximimzer on a line segement
    described by the iterable of xy-coordinates points from the center c.
    Specifically, we return the critical points and the value of

    ||p - c||**2

    for p in points with c being the center of a circle.
    '''

    # Placeholder for dimensions.
    XX, YY = points
    M, N = XX.shape

    # Set X0, X1, Y0, and Y1 depending on the section.
    # XLL, YLL = XX[0:M-1, 0:N-1], YY[0:M-1, 0:N-1]
    # XLR, YLR = XX[1:M  , 0:N-1], YY[1:M  , 0:N-1]
    # XUL, YUL = XX[0:M-1, 1:N  ], YY[0:M-1, 1:N  ]
    # XUR, YUR = XX[1:M  , 1:N  ], YY[1:M  , 1:N  ]

    if where == 'b':
        X0 = XX[0:M-1, 0:N-1]
        X1 = XX[1:M  , 0:N-1]
        Y0 = YY[0:M-1, 0:N-1]
        Y1 = YY[1:M  , 0:N-1]
    elif where == 'r':
        X0 = XX[1:M  , 0:N-1]
        X1 = XX[1:M  , 1:N  ]
        Y0 = YY[1:M  , 0:N-1]
        Y1 = YY[1:M  , 1:N  ]
    elif where == 't':
        X0 = XX[1:M  , 1:N  ]
        X1 = XX[0:M-1, 1:N  ]
        Y0 = YY[1:M  , 1:N  ]
        Y1 = YY[0:M-1, 1:N  ]
    elif where == 'l':
        X0 = XX[0:M-1, 1:N  ]
        X1 = XX[0:M-1, 0:N-1]
        Y0 = YY[0:M-1, 1:N  ]
        Y1 = YY[0:M-1, 0:N-1]

    # Compute the incremental differences in the coordinates.
    DX = X1 - X0
    DY = Y1 - Y0

    # Unpack the center.
    cx, cy = c

    #-------------------------------------------------------------------------
    # Tackle the lower left points.
    #-------------------------------------------------------------------------

    # Function of t to minimize.
    f = lambda tau : (DX * tau + X0 - cx)**2 + (DY * tau + Y0 - cy)**2

    # Times for the critical points.
    TC = ((cx - X0) * DX + (cy - Y0) * DY) / (DX**2 + DY**2)

    # Array for all times.
    T = np.zeros((3, M-1, N-1))
    T[1,:,:] = TC
    T[2,:,:] = 1
    T[1, TC < 0] = 0
    T[1, TC > 1] = 1

    # Evaluate the critical points on this line.
    F = f(T)

    # Get the times of the minimum and maximum values for our small optimization problem.
    amin = np.argmin(F, axis=0)
    amax = np.argmax(F, axis=0)

    # Matrix indices needed for computing the minimum and maximim critical points and values.
    II, JJ = np.meshgrid(range(M-1), range(N-1), indexing='ij')

    # Set the critical points.
    TMIN = T[amin, II, JJ]
    TMAX = T[amax, II, JJ]

    # Set the critical values.
    FMIN = F[amin, II, JJ]
    FMAX = F[amax, II, JJ]

    # Compute the coordinates nearest to and furthest from the center c.
    xc = DX * TMIN + X0
    yc = DY * TMIN + Y0
    xf = DX * TMAX + X0
    yf = DY * TMAX + Y0

    # Return the coordinates with their squared distances from c.
    return (xc, yc, FMIN), (xf, yf, FMAX)


def get_extreme_points(c, points):
    '''
    get_extreme_points_on_shape(c, points)

    Get the locations of the minimizer and maximimzer on a (polygonal)
    shape described by the iterable of xy-coordinates points from the
    center c. Specifically, we return the critical points and the value of

    ||p - c||**2

    for p in points with c being the center of a circle.
    '''

    M, N = points[0].shape

    XC = np.zeros((4, M-1, N-1))
    YC = np.zeros((4, M-1, N-1))
    XF = np.zeros((4, M-1, N-1))
    YF = np.zeros((4, M-1, N-1))
    closest = np.zeros((4, M-1, N-1))
    furthest = np.zeros((4, M-1, N-1))

    sides = ('b', 'r', 't', 'l')

    for i, s in enumerate(sides):
        pc, pf = get_extreme_points_on_lines(c, points, where=s)

        XC[i,:,:], YC[i,:,:], closest[i,:,:] = pc
        XF[i,:,:], YF[i,:,:], furthest[i,:,:] = pf

    # Matrix indices needed for computing the minimum and maximim critical points and values.
    II, JJ = np.meshgrid(range(M-1), range(N-1), indexing='ij')

    amin = np.argmin(closest, axis=0)
    amax = np.argmax(furthest, axis=0)

    return (XC[amin, II, JJ], YC[amin, II, JJ], closest[amin, II, JJ]), \
           (XF[amax, II, JJ], YF[amax, II, JJ], furthest[amax, II, JJ])


#-------------------------------------------------------
# Methods for determining if a polygonal shape is within
# a circle.
#-------------------------------------------------------

def is_shape_in_circle(c, r, b, points):
    rsq = (r - b / 2)**2
    closest, furthest = get_extreme_points(c, points)

    a, b = rsq - closest[-1], furthest[-1] - rsq

    return (a > 0) & (b <= a)


def is_square_in_circular_aliasing_zone(c, r, b, points):
    rmbsq = (r - b / 2)**2
    rpbsq = (r + b / 2)**2
    closest, furthest = get_extreme_points(c, points)

    am, bm = rmbsq - closest[-1], furthest[-1] - rmbsq
    ap, bp = rpbsq - closest[-1], furthest[-1] - rpbsq

    exceeds_inner_radius = (am <= 0) | (bm > am)
    does_not_exceed_outer_radius = (ap > 0) & (bp <= ap)

    return exceeds_inner_radius & does_not_exceed_outer_radius


def is_shape_outside_circle(c, r, b, points):
    rsq = (r + b / 2)**2
    closest, _ = get_extreme_points(c, points)

    a = rsq - closest[-1]

    return a < 0


#-------------------------------------------------------
# Methods for determining if a polygonal shape is within
# an annulus of inner radius ri and outer radius ro.
#-------------------------------------------------------

def is_shape_in_annulus(c, ri, ro, b, points):
    risq = (ri + b / 2)**2
    rosq = (ro - b / 2)**2
    closest, furthest = get_extreme_points(c, points)

    am, bm = risq - closest[-1], furthest[-1] - risq
    ap, bp = rosq - closest[-1], furthest[-1] - rosq

    exceeds_inner_radius = (am <= 0) | (bm > am)
    does_not_exceed_outer_radius = (ap > 0) & (bp <= ap)

    return exceeds_inner_radius & does_not_exceed_outer_radius


def is_shape_in_annulus_aliasing_zone(c, ri, ro, b, points):
    inner_aliasing_region = is_shape_in_annulus(c, ri - b/2, ri + b/2, 0, points)
    outer_aliasing_region = is_shape_in_annulus(c, ro - b/2, ro + b/2, 0, points)

    return inner_aliasing_region, outer_aliasing_region


def is_shape_outside_annulus(c, ri, ro, b, points):
    risq = (ri - b / 2)**2
    rosq = (ro + b / 2)**2
    closest, furthest = get_extreme_points(c, points)

    am, bm = risq - closest[-1], furthest[-1] - risq
    ap, bp = rosq - closest[-1], furthest[-1] - rosq

    does_not_exceed_inner_radius = (am > 0) | (bm < am)
    exceeds_outer_radius = ap < bp

    return does_not_exceed_inner_radius | exceeds_outer_radius


#-------------------------------------------------------
#
#-------------------------------------------------------

def radius_to_cellspan(radius, width, height, xn, x0, yn, y0):
    cpux = width / (xn - x0)
    cpuy = height / (yn - y0)

    return int(np.ceil(radius * np.sqrt(cpux**2 + cpuy**2)))


#-------------------------------------------------------
# Methods for computing integrals on circular and
# annular domains.
#-------------------------------------------------------

def compute_circular_integral(F, c, r, b, width, height, x0, xn, y0, yn):
    # Dimensions of the state function.
    M, N = F.shape

    # Get the ijth coordinate the point that is in the lower left corner
    # of the shape containing the center c.
    i, j = point_to_cell_ij(c, x0, xn, y0, yn, width, height)
    # print('-D- ci = {0}, cj = {1}'.format(i, j))

    # Convert the radial distance from the center c to the number of cells
    # in the grid spanned by the circle, and compute the number of indices
    # over which we need to run.
    cspan = radius_to_cellspan(r + b/2, width, height, x0, xn, y0, yn)
    I = [i - cspan, i + cspan]
    J = [j - cspan, j + cspan]
    num_indices = 2 * cspan + 1
    # print('-D- cspan       =', cspan)
    # print('-D- I           =', *I)
    # print('-D- J           =', *J)
    # print('-D- num_indices =', num_indices)

    # Compute the step size (units per cell).
    dx = (xn - x0) / width
    dy = (xn - x0) / height

    # Compute a locale grid of points.
    xl, xr = x0 + I[0] * dx, x0 + I[-1] * dx
    yl, yr = y0 + J[0] * dy, y0 + J[-1] * dy


    II, JJ = np.meshgrid(range(I[0], I[-1]),
                         range(J[0], J[-1]))
    II = (II + M) % M
    JJ = (JJ + N) % M

    FR = F[II, JJ]

    XX, YY = np.meshgrid(np.linspace(xl, xr, num_indices),
                         np.linspace(yl, yr, num_indices))
    P, Q = XX.shape
    P, Q = P-1, Q-1

    # Mask for the inner region.
    inmask = is_shape_in_circle(c, ri, b, (XX, YY))

    exact_area = np.pi * r**2
    integral = 0
    diffarea = dx * dy

    integral = np.sum(FR[inmask]) * diffarea
    # print('-D- inmask.shape =', inmask.shape)
    # print('-D- np.sum(FR[inmask]) =', np.sum(FR[inmask]))
    # print('-D- diffarea =', diffarea)
    # print('-D- integral =', integral)

    if b > 0:
        # Mask for the aliasing region.
        almask = is_square_in_circular_aliasing_zone(c, ri, b, (XX, YY))

        # Compute the scaling for the
        XM = XX[:P,:Q] + 0.5 * dx
        YM = YY[:P,:Q] + 0.5 * dy
        LD = np.sqrt((XM - c[0])**2 + (YM - c[1])**2)
        scaling = diffarea * ((r + 0.5 * b - LD) / b)

        # Compute the integral.
        integral += np.sum(FR[almask] * scaling[almask])

    return integral, exact_area


def compute_annular_integral(F, c, ri, ro, b, width, height, x0, xn, y0, yn):
    # Dimensions of the state function.
    M, N = F.shape

    # Get the ijth coordinate the point that is in the lower left corner
    # of the shape containing the center c.
    i, j = point_to_cell_ij(c, x0, xn, y0, yn, width, height)
    # print('-D- ci = {0}, cj = {1}'.format(i, j))

    # Convert the radial distance from the center c to the number of cells
    # in the grid spanned by the circle, and compute the number of indices
    # over which we need to run.
    cspan = radius_to_cellspan(ro + b/2, width, height, x0, xn, y0, yn)
    I = [i - cspan, i + cspan]
    J = [j - cspan, j + cspan]
    num_indices = 2 * cspan + 1
    # print('-D- cspan       =', cspan)
    # print('-D- I           =', *I)
    # print('-D- J           =', *J)
    # print('-D- num_indices =', num_indices)

    # Compute the step size (units per cell).
    dx = (xn - x0) / width
    dy = (xn - x0) / height

    # Compute a locale grid of points.
    xl, xr = x0 + I[0] * dx, x0 + I[-1] * dx
    yl, yr = y0 + J[0] * dy, y0 + J[-1] * dy


    II, JJ = np.meshgrid(range(I[0], I[-1]),
                         range(J[0], J[-1]))
    II = np.mod(II + M, M)
    JJ = np.mod(JJ + N, N)

    FR = F[II, JJ]

    XX, YY = np.meshgrid(np.linspace(xl, xr, num_indices),
                         np.linspace(yl, yr, num_indices))
    P, Q = XX.shape
    P, Q = P-1, Q-1

    # Mask for the inner region.
    inmask = is_shape_in_annulus(c, ri, ro, b, (XX, YY))

    exact_area = np.pi * (ro**2 - ri**2)
    diffarea = dx * dy

    integral = np.sum(FR[inmask]) * diffarea
    # print('-D- inmask.shape =', inmask.shape)
    # print('-D- np.sum(FR[inmask]) =', np.sum(FR[inmask]))
    # print('-D- diffarea =', diffarea)
    # print('-D- integral =', integral)

    if b > 0:
        # Mask for the aliasing region.
        iamask, oamask = is_shape_in_annulus_aliasing_zone(c, ri, ro, b, (XX, YY))

        # Compute the midpoints and distances of the midpoint from the center.
        XM = XX[:P,:Q] + 0.5 * dx
        YM = YY[:P,:Q] + 0.5 * dy
        LD = np.sqrt((XM - c[0])**2 + (YM - c[1])**2)

        # Compute the per-cell scaling for the outer aliasing region.
        scaling = (ro + 0.5 * b - LD) / b

        # Add the contribution from the outer aliasing region.
        integral += np.sum(FR[oamask] * scaling[oamask]) * diffarea

        # Compute the per-cell scaling for the inner aliasing region.
        scaling = (LD - (ri - 0.5 * b)) / b

        # Add the contribution from the inner aliasing region.
        integral += np.sum(FR[iamask] * scaling[iamask]) * diffarea

    return integral, exact_area


#-------------------------------------------------------
# Helper methods.
#-------------------------------------------------------

def point_to_cell_ij(p, x0, xn, y0, yn, width, height):
    i = int((p[0] - x0) / (xn - x0) * width)
    j = int((p[1] - y0) / (yn - y0) * height)

    return i, j

def radius_to_cellspan(r, width, height, x0, xn, y0, yn):
    dx = (xn - x0) / width
    dy = (yn - y0) / height

    return int(np.round(max(r / dx, r / dy)))

#
# Methods for testing out drawing and masks.
#

def draw_circular_region(x0=-2, xn=2, numx=33, y0=-2, yn=2, numy=33, c=(0,0), ri=0.5, b=0.5):
    x = np.linspace(x0, xn, numx)
    dx = x[1] - x[0]

    y = np.linspace(y0, yn, numy)
    dy = y[1] - y[0]

    XX, YY = np.meshgrid(x, y)

    # The state function.
    np.random.seed(0)
    #F = np.random.rand(XX.shape[0] - 1, XX.shape[1] - 1)
    F = np.ones((XX.shape[0] - 1, XX.shape[1] - 1))

    fig, ax = plt.subplots()

    # Plot the circular and anti-aliasing curves.
    t = np.linspace(0, 2*np.pi, 100)

    # Draw the grid lines
    for i in range(numx):
        ax.plot([x[i], x[i]], [y0, yn], color='gray', linestyle='-', alpha=0.2)

    for j in range(numy):
        ax.plot([x0, xn], [y[j], y[j]], color='gray', linestyle='-', alpha=0.2)

    # Draw the circular boundaries.
    colors = ['blue', 'green', 'red']
    radii = [ri - b/2, ri, ri + b/2]
    linestyles = ['-', '-', '-']

    for clr, rad, ls in zip(colors, radii, linestyles):
        ax.plot(c[0] + rad * np.cos(t), c[1] + rad * np.sin(t), color=clr, linestyle=ls)

    # Mask for the inner region.
    inmask = is_shape_in_circle(c, ri, b, (XX, YY))

    # Mask for the aliasing region.
    almask = is_square_in_circular_aliasing_zone(c, ri, b, (XX, YY))

    # The approximate integral.
    integral, exact_area = compute_circular_integral(F=F, c=c, r=ri, b=b, width=numx-1, height=numy-1, x0=x0, xn=xn, y0=y0, yn=yn)

    for i in range(inmask.shape[0]):
        for j in range(inmask.shape[1]):
            if inmask[i, j]:
                px = [XX[i,j], XX[i+1,j], XX[i+1,j+1], XX[i,j+1]]
                py = [YY[i,j], YY[i+1,j], YY[i+1,j+1], YY[i,j+1]]

                ax.fill(px, py, color=colors[0], alpha=0.25)

            if almask[i, j]:
                px = [XX[i,j], XX[i+1,j], XX[i+1,j+1], XX[i,j+1]]
                py = [YY[i,j], YY[i+1,j], YY[i+1,j+1], YY[i,j+1]]

                ax.fill(px, py, color=colors[1], alpha=0.25)

    title_str = 'Circular region, $r_i \\approx$ {0:.3f}'.format(ri)
    title_str += '\nApproximate area: {0:.10f}'.format(integral)
    title_str += '\nExact area: {0:.10f}'.format(exact_area)
    title_str += '\nRelative error: {0:.4f}%'.format((integral - exact_area) / exact_area * 100)
    ax.set_title(title_str)
    ax.set_xlim((x0, xn))
    ax.set_ylim((y0, yn))
    plt.show()


def draw_annular_region(x0=-2, xn=2, numx=33, y0=-2, yn=2, numy=33, c=(0,0), ri=0.5, rscale=3, b=0.25):
    x = np.linspace(x0, xn, numx)
    dx = x[1] - x[0]

    y = np.linspace(y0, yn, numy)
    dy = y[1] - y[0]

    XX, YY = np.meshgrid(x, y)

    # The state function.
    np.random.seed(0)
    #F = np.random.rand(XX.shape[0] - 1, XX.shape[1] - 1)
    F = np.ones((XX.shape[0] - 1, XX.shape[1] - 1))

    ro = rscale * ri

    fig, ax = plt.subplots()

    # Plot the circular and anti-aliasing curves.
    t = np.linspace(0, 2*np.pi, 100)

    # Draw the grid lines
    for i in range(numx):
        ax.plot([x[i], x[i]], [y0, yn], color='gray', linestyle='-', alpha=0.2)

    for j in range(numy):
        ax.plot([x0, xn], [y[j], y[j]], color='gray', linestyle='-', alpha=0.2)

    # Draw the circular boundaries.
    colors = ['blue', 'green', 'blue', 'red', 'green', 'red']
    radii = [ri - b/2, ri, ri + b/2, ro - b/2, ro, ro + b/2]
    linestyles = ['-'] * len(radii)

    for clr, rad, ls in zip(colors, radii, linestyles):
        ax.plot(c[0] + rad * np.cos(t), c[1] + rad * np.sin(t), color=clr, linestyle=ls)

    inmask = is_shape_in_annulus(c, ri, ro, b, (XX, YY))
    iamask, oamask = is_shape_in_annulus_aliasing_zone(c, ri, ro, b, (XX, YY))

    # The approximate integral.
    integral, exact_area = compute_annular_integral(F=F, c=c, ri=ri, ro=ro, b=b, width=numx-1, height=numy-1, x0=x0, xn=xn, y0=y0, yn=yn)
    approx_average = integral / exact_area
    exact_average = 1

    for i in range(inmask.shape[0]):
        for j in range(inmask.shape[1]):
            if inmask[i, j]:
                px = [XX[i,j], XX[i+1,j], XX[i+1,j+1], XX[i,j+1]]
                py = [YY[i,j], YY[i+1,j], YY[i+1,j+1], YY[i,j+1]]

                ax.fill(px, py, color='green', alpha=0.25)

            if iamask[i, j]:
                px = [XX[i,j], XX[i+1,j], XX[i+1,j+1], XX[i,j+1]]
                py = [YY[i,j], YY[i+1,j], YY[i+1,j+1], YY[i,j+1]]

                ax.fill(px, py, color='blue', alpha=0.25)

            if oamask[i, j]:
                px = [XX[i,j], XX[i+1,j], XX[i+1,j+1], XX[i,j+1]]
                py = [YY[i,j], YY[i+1,j], YY[i+1,j+1], YY[i,j+1]]

                ax.fill(px, py, color='red', alpha=0.25)

    title_str = 'Annular region, $r_i \\approx$ {0:.3f}, $r_o \\approx$ {1:.3f}'.format(ri, ro)
    title_str += '\nApproximate average: {0:.10f}'.format(approx_average)
    title_str += '\nRelative error: {0:.4f}%'.format((approx_average - 1) * 100)
    ax.set_title(title_str)
    ax.set_xlim((x0, xn))
    ax.set_ylim((y0, yn))
    plt.show()


def benchmark(ri, ro, b, x0, xn, y0, yn, numcellsx=[16, 32, 64, 128, 256], numcellsy=None, show=True):
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

# Main section for debug.
if __name__ == '__main__':
    x0, xn = -1, 1
    y0, yn = -1, 1
    ri = 0.05 * min((xn - x0), (yn - y0))
    rscale = 3
    ro = ri * rscale
    b = 0.5 * (ro - ri)
    numcellsx = list(range(16, 129, 8))

    benchmark(ri, ro, b, x0, xn, y0, yn, numcellsx=numcellsx, show=False)





