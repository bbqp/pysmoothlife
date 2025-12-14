# grid.py
#
# A utility library for finding areas of chapes on a grid.

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import multiprocessing
from transitionfunction import TransitionFunction
import itertools

np.set_printoptions(linewidth=200)


class SmoothLifeGrid:

    def __init__(self, width=640, height=640, x0=-1, xn=1, y0=-1, yn=1, numcellsw=64, numcellsh=64, ri_ratio=0.1, rad_scale=3, b=0.1/4,
                 bint=[0.278, 0.365], dint=[0.267, 0.445], alphan=0.028, alpham=0.147, dtype=np.float64):
        self.dtype = dtype

        self.x0 = self.dtype(x0)
        self.xn = self.dtype(xn)
        self.numx = numcellsw + 1

        self.y0 = self.dtype(y0)
        self.yn = self.dtype(yn)
        self.numy = numcellsh + 1

        self.width = width
        self.height = height

        self.numcellsh = numcellsh
        self.numcellsw = numcellsw

        self.dx = self.dtype(self.xn - self.x0) / self.numcellsw
        self.dy = self.dtype(self.yn - self.y0) / self.numcellsh

        # The state function.
        self.F = np.zeros((self.numcellsw, self.numcellsh))

        # The grid of points.
        self.x = np.linspace(self.x0, self.xn, self.numcellsw + 1)
        self.y = np.linspace(self.y0, self.yn, self.numcellsh + 1)
        self.XX, self.YY = np.meshgrid(self.x, self.y, indexing='ij')

        # The inner and outer radii.
        self.ri_ratio = ri_ratio
        self.ri = ri_ratio * min(xn - x0, yn - y0)
        self.rad_scale = rad_scale
        self.ro = self.ri * rad_scale
        self.b = b

        # Span for the (outer) radius.
        self.cspani = self.radius_to_cellspan(self.ri + self.b / 2)
        self.cspano = self.radius_to_cellspan(self.ro + self.b / 2)

        # Set the birth and death intervals.
        self.bint = np.array(bint, dtype=self.dtype)
        self.dint = np.array(dint, dtype=self.dtype)

        # Define the transition function.
        self.transition_function = TransitionFunction(alpham, alphan, bint, dint, dtype=self.dtype)

        # Set coffiecients for updating the state functions.
        self.alpham = self.dtype(alpham)
        self.alphan = self.dtype(alphan)

        self._initialize_domains_masks_indices()

    def _initialize_domains_masks_indices(self):
        # Create the circular domain.
        xl, xr = -self.dtype(self.cspani + 0.5) * self.dx, self.dtype(self.cspani + 0.5) * self.dx
        yl, yr = -self.dtype(self.cspani + 0.5) * self.dy, self.dtype(self.cspani + 0.5) * self.dy
        self.cdom = np.meshgrid(np.linspace(xl, xr, 2*self.cspani + 2),
                                np.linspace(yl, yr, 2*self.cspani + 2),
                                indexing='ij')

        # Create the annular domain.
        xl, xr = -self.dtype(self.cspano + 0.5) * self.dx, self.dtype(self.cspano + 0.5) * self.dx
        yl, yr = -self.dtype(self.cspano + 0.5) * self.dy, self.dtype(self.cspano + 0.5) * self.dy
        self.adom = np.meshgrid(np.linspace(xl, xr, 2*self.cspano + 2),
                                np.linspace(yl, yr, 2*self.cspano + 2),
                                indexing='ij')

        # Create the midpoint grids for each domain.
        MC, NC = self.cdom[0].shape
        self.cmid = self.cdom[0][:MC-1,:NC-1] + self.dx / 2, \
                    self.cdom[1][:MC-1,:NC-1] + self.dy / 2

        MA, NA = self.adom[0].shape
        self.amid = self.adom[0][:MA-1,:NA-1] + self.dx / 2, \
                    self.adom[1][:MA-1,:NA-1] + self.dy / 2

        # Compute the distances of the midpoints from the center (0,0).
        self.cld = np.sqrt(self.cmid[0]**2 + self.cmid[1]**2)
        self.ald = np.sqrt(self.amid[0]**2 + self.amid[1]**2)

        # Compute the scalings for the anti-aliasing zones. First the circle,
        # and then for the inner and outer annalus anti-aliasing zones.
        self.cascaling = self.dtype(self.ri + 0.5 * self.b - self.cld) / self.b
        self.aiscaling = self.dtype(self.ald - (self.ri - 0.5 * self.b)) / self.b
        self.aoscaling = self.dtype(self.ro + 0.5 * self.b - self.ald) / self.b

        # The origin of the reference domain.
        center = np.array([0, 0], dtype=self.dtype)

        # Create the circular domain masks. Note that the dimension in each direction
        # is one less than the meshgrid of the total domain.
        self.cimask = self.is_shape_in_circle(center, self.ri, self.b, self.cdom)
        self.camask = self.is_shape_in_circular_aliasing_zone(center, self.ri, self.b, self.cdom)
        self.cmask = self.cimask | self.camask
        
        # Compute the circular scaling as one array.
        self.cscaling = np.zeros((MC-1, NC-1), dtype=self.dtype)
        self.cscaling[self.cimask] = self.dtype(1)
        self.cscaling[self.camask] = self.cascaling[self.camask]

        # Create the annular domain masks. Note that the dimension in each direction
        # is one less than the meshgrid of the total domain.
        self.aimask = self.is_shape_in_annulus(center, self.ri, self.ro, self.b, self.adom)
        self.aamask_inner, self.aamask_outer = self.is_shape_in_annulus_aliasing_zone(center, self.ri, self.ro, self.b, self.adom)
        self.amask = self.aimask | self.aamask_inner | self.aamask_outer
        
        self.ascaling = np.zeros((MA-1, NA-1), dtype=self.dtype)
        self.ascaling[self.aimask] = self.dtype(1)
        self.ascaling[self.aamask_inner] = self.aiscaling[self.aamask_inner]
        self.ascaling[self.aamask_outer] = self.aoscaling[self.aamask_outer]

        # Arrays that hold the integrals for the circles and annuli in each cell.
        self.CI = np.zeros((self.numcellsw, self.numcellsh), dtype=self.dtype)
        self.AI = np.zeros((self.numcellsw, self.numcellsh), dtype=self.dtype)
        
        # Index grids for the circular and annular domains.
        circ_range = range(-self.cspani, self.cspani + 1)
        self.CII, self.CJJ = np.meshgrid(circ_range, circ_range, indexing='ij')

        annulus_range = range(-self.cspano, self.cspano + 1)
        self.AII, self.AJJ = np.meshgrid(annulus_range, annulus_range, indexing='ij')

    def set_state_function(self, state):
        self.F[:,:] = self.dtype(state)

    def radius_to_cellspan(self, radius):
        return int(np.round(self.dtype(radius) / min(self.dx, self.dy)))

    def update(self):
        self.compute_integrals()
        self.F[:,:] = self.transition_function.__call__(self.AI, self.CI)

    #-------------------------------------------------------
    # Methods for finding extreme points on polygonal
    # shapes.
    #-------------------------------------------------------

    def get_extreme_points_on_line(self, c, points, where='b'):
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
        T = np.zeros((3, M-1, N-1), dtype=self.dtype)
        T[1,:,:] = TC
        T[2,:,:] = self.dtype(1)
        T[1, TC < 0] = self.dtype(0)
        T[1, TC > 1] = self.dtype(1)

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

    def get_extreme_points(self, c, points):
        '''
        get_extreme_points_on_shape(c, points)

        Get the locations of the minimizer and maximimzer on a (polygonal)
        shape described by the iterable of xy-coordinates points from the
        center c. Specifically, we return the critical points and the value of

        ||p - c||**2

        for p in points with c being the center of a circle.
        '''

        M, N = points[0].shape

        XC = np.zeros((4, M-1, N-1), dtype=self.dtype)
        YC = np.zeros((4, M-1, N-1), dtype=self.dtype)
        XF = np.zeros((4, M-1, N-1), dtype=self.dtype)
        YF = np.zeros((4, M-1, N-1), dtype=self.dtype)
        closest = np.zeros((4, M-1, N-1), dtype=self.dtype)
        furthest = np.zeros((4, M-1, N-1), dtype=self.dtype)

        sides = ('b', 'r', 't', 'l')

        for i, s in enumerate(sides):
            pc, pf = self.get_extreme_points_on_line(c, points, where=s)

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

    def is_shape_in_circle(self, c, r, b, points):
        rsq = self.dtype(r - b / 2)**2
        closest, furthest = self.get_extreme_points(c, points)

        a, b = rsq - closest[-1], furthest[-1] - rsq

        return (a > 0) & (b <= a)

    def is_shape_in_circular_aliasing_zone(self, c, r, b, points):
        rmbsq = self.dtype(r - b / 2)**2
        rpbsq = self.dtype(r + b / 2)**2
        closest, furthest = self.get_extreme_points(c, points)

        am, bm = rmbsq - closest[-1], furthest[-1] - rmbsq
        ap, bp = rpbsq - closest[-1], furthest[-1] - rpbsq

        exceeds_inner_radius = (am <= 0) | (bm > am)
        does_not_exceed_outer_radius = (ap > 0) & (bp <= ap)

        return exceeds_inner_radius & does_not_exceed_outer_radius

    def is_shape_outside_circle(self, c, r, b, points):
        rsq = self.dtype(r + b / 2)**2
        closest, _ = self.get_extreme_points(c, points)

        a = rsq - closest[-1]

        return a < 0


    #-------------------------------------------------------
    # Methods for determining if a polygonal shape is within
    # an annulus of inner radius ri and outer radius ro.
    #-------------------------------------------------------

    def is_shape_in_annulus(self, c, ri, ro, b, points):
        risq = self.dtype(ri + b / 2)**2
        rosq = self.dtype(ro - b / 2)**2
        closest, furthest = self.get_extreme_points(c, points)

        am, bm = risq - closest[-1], furthest[-1] - risq
        ap, bp = rosq - closest[-1], furthest[-1] - rosq

        exceeds_inner_radius = (am <= 0) | (bm > am)
        does_not_exceed_outer_radius = (ap > 0) & (bp <= ap)

        return exceeds_inner_radius & does_not_exceed_outer_radius

    def is_shape_in_annulus_aliasing_zone(self, c, ri, ro, b, points):
        inner_aliasing_region = self.is_shape_in_annulus(c, ri - b/2, ri + b/2, 0, points)
        outer_aliasing_region = self.is_shape_in_annulus(c, ro - b/2, ro + b/2, 0, points)

        return inner_aliasing_region, outer_aliasing_region

    def is_shape_outside_annulus(self, c, ri, ro, b, points):
        risq = self.dtype(ri - b / 2)**2
        rosq = self.dtype(ro + b / 2)**2
        closest, furthest = self.get_extreme_points(c, points)

        am, bm = risq - closest[-1], furthest[-1] - risq
        ap, bp = rosq - closest[-1], furthest[-1] - rosq

        does_not_exceed_inner_radius = (am > 0) | (bm < am)
        exceeds_outer_radius = ap < bp

        return does_not_exceed_inner_radius | exceeds_outer_radius


    #-------------------------------------------------------
    # Methods for computing integrals on circular and
    # annular domains.
    #-------------------------------------------------------

    def compute_circular_integral(self, i, j):
        # Dimensions of the state function.
        M, N = self.F.shape

        # Get the value of the state function restricted to the circular subgrid.
        II = np.mod(self.CII + i, M)
        JJ = np.mod(self.CJJ + j, N)
        FR = self.F[II, JJ]

        exact_area = self.dtype(np.pi) * self.ri**2
        diffarea = self.dx * self.dy

        # integral = np.sum(FR[self.cimask]) * diffarea

        # # Add the contributions from the aliasing region.
        # if self.b > 0:
            # integral += np.sum(FR[self.camask] * self.cascaling[self.camask]) * diffarea

        integral = np.sum(FR[self.cmask] * self.cscaling[self.cmask]) * diffarea

        return integral / exact_area


    def compute_annular_integral(self, i, j):
        # Dimensions of the state function.
        M, N = self.F.shape
        
        II = np.mod(self.AII + i, M)
        JJ = np.mod(self.AJJ + j, N)

        FR = self.F[II, JJ]

        exact_area = self.dtype(np.pi) * (self.ro**2 - self.ri**2)
        diffarea = self.dx * self.dy
        integral = np.sum(FR[self.amask] * self.ascaling[self.amask]) * diffarea

        return integral / exact_area

    def compute_integrals(self, serial=True, nprocs=multiprocessing.cpu_count()):
        M, N = self.numcellsw, self.numcellsh

        if serial:
            for i in range(M):
                for j in range(N):
                    self.CI[i, j] = self.compute_circular_integral(i, j)
                    self.AI[i, j] = self.compute_annular_integral(i, j)
        else:
            ci_cartprod = itertools.product(range(M), range(N))
            ai_cartprod = itertools.product(range(M), range(N))

            #with multiprocessing.Pool(processes=nprocs) as p:
            with multiprocessing.ThreadPool(processes=nprocs) as p:
                ci_results = p.starmap_async(self.compute_circular_integral, ci_cartprod)
                ai_results = p.starmap_async(self.compute_annular_integral, ai_cartprod)

                self.CI[:,:] = np.reshape(ci_results.get(), (M, N))
                self.AI[:,:] = np.reshape(ai_results.get(), (M, N))

    #-------------------------------------------------------
    # Methods for testing out drawing and masks.
    #-------------------------------------------------------

    def draw_circular_region(self):
        M, N = self.XX.shape
        XX, YY = self.XX[:M-1,:N-1] + self.dx/2, self.YY[:M-1,:N-1] + self.dy/2
        #XX, YY = self.cdom
        x = XX[:, 0]
        y = YY[0, :]
        x0, xn = x[0], x[-1]
        y0, yn = y[0], y[-1]
        numx = len(x)
        numy = len(y)

        # Center of the circle.
        c = (0, 0)

        # Inner radius for the circle.
        ri = self.ri

        # The figure for the circular region.
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
        radii = [self.ri - self.b/2, self.ri, self.ri + self.b/2]
        linestyles = ['-', '-', '-']

        for clr, rad, ls in zip(colors, radii, linestyles):
            ax.plot(c[0] + rad * np.cos(t), c[1] + rad * np.sin(t), color=clr, linestyle=ls)

        inmask = self.is_shape_in_circle(c, self.ri, self.b, (XX, YY))
        camask = self.is_shape_in_circular_aliasing_zone(c, self.ri, self.b, (XX, YY))

        for i in range(inmask.shape[0]):
            for j in range(inmask.shape[1]):
                if inmask[i, j]:
                    px = [XX[i,j], XX[i+1,j], XX[i+1,j+1], XX[i,j+1]]
                    py = [YY[i,j], YY[i+1,j], YY[i+1,j+1], YY[i,j+1]]

                    ax.fill(px, py, color=colors[0], alpha=0.25)

                if camask[i, j]:
                    px = [XX[i,j], XX[i+1,j], XX[i+1,j+1], XX[i,j+1]]
                    py = [YY[i,j], YY[i+1,j], YY[i+1,j+1], YY[i,j+1]]

                    ax.fill(px, py, color=colors[1], alpha=0.25)

        # The approximate integral.
        i = int(np.round((c[0] - x0) / self.dx))
        j = int(np.round((c[1] - y0) / self.dy))
        approx_average = self.compute_circular_integral(i, j)

        title_str = 'Circular region, $r_i \\approx$ {0:.3f}'.format(ri)
        title_str += '\nApproximate average: {0:.10f}'.format(approx_average)
        title_str += '\nRelative error: {0:.4f}%'.format((approx_average - 1) * 100)
        ax.set_title(title_str)
        ax.set_xlim((x0, xn))
        ax.set_ylim((y0, yn))
        plt.show()


    def draw_annular_region(self):
        M, N = self.XX.shape
        XX, YY = self.XX[:M-1,:N-1] + self.dx/2, self.YY[:M-1,:N-1] + self.dy/2
        x = XX[:, 0]
        y = YY[0, :]
        x0, xn = x[0], x[-1]
        y0, yn = y[0], y[-1]
        numx = len(x)
        numy = len(y)

        # Center of the circle.
        c = (0, 0)

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
        radii = [self.ri - self.b/2, self.ri, self.ri + self.b/2,
                 self.ro - self.b/2, self.ro, self.ro + self.b/2]
        linestyles = ['-'] * len(radii)

        for clr, rad, ls in zip(colors, radii, linestyles):
            ax.plot(c[0] + rad * np.cos(t), c[1] + rad * np.sin(t), color=clr, linestyle=ls)

        inmask = self.is_shape_in_annulus(c, self.ri, self.ro, self.b, (XX, YY))
        iamask, oamask = self.is_shape_in_annulus_aliasing_zone(c, self.ri, self.ro, self.b, (XX, YY))

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

        # The approximate integral.
        i = int(np.round((c[0] - x0) / self.dx))
        j = int(np.round((c[1] - y0) / self.dy))
        approx_average = self.compute_annular_integral(i, j)

        title_str = 'Annular region, $r_i \\approx$ {0:.3f}, $r_o \\approx$ {1:.3f}'.format(self.ri, self.ro)
        title_str += '\nApproximate average: {0:.10f}'.format(approx_average)
        title_str += '\nRelative error: {0:.4f}%'.format((approx_average - 1) * 100)
        ax.set_title(title_str)
        ax.set_xlim((x0, xn))
        ax.set_ylim((y0, yn))
        plt.show()



