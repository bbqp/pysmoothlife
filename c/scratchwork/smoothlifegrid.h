/*
 * smoothlifegrid.c
 *
 * Module for speeding things up since python double loops are slow.
 */

#include <concept>
#include <vector>
#include <algorithm>

#include "interval.h"
#include "interval.h"
#include "meshgrid.h"
#import "utils.h"

template class SmoothLifeGrid<std::floating_point T>
{
	private:

	T x0, xn;
	T y0, yn;
	int numcellsw, numcellsh;
	int numx, numy;
	T dx, dy;
	T ri_ratio;
	T rad_scale;
	T b;
	T bint[2];
	T dint[2];
	T alpham, alpham;
	Interval<T> x, y;
	Meshgrid<T> domain;
	Meshgrid<T> cdom, adom;
	Meshgrid<T> cmid, amid;
		
	public:
	
	SmoothLifeGrid(T x0=-1, T xn=1, T y0=-1, T yn=1, int numcellsw=64, int numcellsh=64,
				   T ri_ratio=0.1, T rad_scale=3, T b=0.1/4,
                   bint=std::tuple(0.278, 0.365), std::tuple<T>(0.267, 0.445),
				   T alphan=0.028, T alpham=0.147)
	{
        this->x0 = x0;
        this->xn = xn;
        this->numx = numcellsw + 1;

        this->y0 = y0;
        this->yn = yn;
        this->numy = numcellsh + 1;

        this->width = width;
        this->height = height;

        this->numcellsh = numcellsh;
        this->numcellsw = numcellsw;

        this->dx = (this->xn - this->x0) / this->numcellsw;
        this->dy = (this->yn - this->y0) / this->numcellsh;

        // The state function.
        this->F = std::vector<T>(this->numcellsw * this->numcellsh);

        // The grid of points.
        this->x = Interval<T>(this->x0, this->xn, this->numcellsw + 1);
        this->y = Interval<T>((this->y0, this->yn, this->numcellsh + 1);
        this->domain = Meshgrid<T>(this->x, this->y);

        // The inner and outer radii.
        this->ri_ratio = ri_ratio;
        this->ri = ri_ratio * std::min(xn - x0, yn - y0);
        this->rad_scale = rad_scale;
        this->ro = this->ri * rad_scale;
        this->b = b;

        // Span for the (outer) radius.
        this->cspani = this->radius_to_cellspan(this->ri + this->b / 2);
        this->cspano = this->radius_to_cellspan(this->ro + this->b / 2);

        // Set the birth and death intervals.
        this->bint = std::make_tuple(bint, dtype=this->dtype);
        this->dint = np.array(dint, dtype=this->dtype);

        // Define the transition function.
        this->transition_function = TransitionFunction(alpham, alphan, bint, dint, dtype=this->dtype);

        // Set coffiecients for updating the state functions.
        this->alpham = alpham;
        this->alphan = alphan;

        this->initialize_domains_masks_indices();
	}
	
	~SmoothLifeGrid() {}
	
	void initialize_domains_masks_indices():
		// Create the circular domain.
		T xl = -(this->cspani + 0.5) * this->dx
		T xr =  (this->cspani + 0.5) * this->dx
		xl, xr = -(this->cspani + 0.5) * this->dx, (this->cspani + 0.5) * this->dx
		yl, yr = -(this->cspani + 0.5) * this->dy, (this->cspani + 0.5) * this->dy
		this->cdom = Meshgrid<T>(xl, xr, 2*this->cspani + 2, yl, yr, 2*this->cspani + 2);

		// Create the annular domain.
		xl = -(this->cspano + 0.5) * this->dx
		xr =  (this->cspano + 0.5) * this->dx
		yl = -(this->cspano + 0.5) * this->dy
		yr =  (this->cspano + 0.5) * this->dy
		this->adom = Meshgrid<T>(xl, xr, 2*this->cspano + 2, yl, yr, 2*this->cspano + 2);

		// Create the midpoint grids for each domain.
		int MC = this->cdom.numx;
		int NC = this->cdom.numy;
		this->cmid = this->cdom.slice_copy_and_offset(0, MC - 1, 0, NC - 1, self->dx / 2, self.dy / 2);

		int MA = this->adom.numx;
		int NA = this->adom.numy;
		this->amid = this->adom.slice_copy_and_offset(0, MA - 1, 0, NA - 1, self->dx / 2, self.dy / 2);

		// Compute the distances of the midpoints from the center (0,0).
		this->cld = distance_from<T>(this->cmid.XX, this->cmid.YY, 0, 0);
		this->ald = distance_from<T>(this->amid.XX, this->amid.YY, 0, 0);

		# Compute the scalings for the anti-aliasing zones. First the circle,
		# and then for the inner and outer annalus anti-aliasing zones.
		this->cascaling = (this->ri + 0.5 * this->b - this->cld) / this->b
		this->aiscaling = (this->ald - (this->ri - 0.5 * this->b)) / this->b
		this->aoscaling = (this->ro + 0.5 * this->b - this->ald) / this->b

		# The origin of the reference domain.
		center = std::pair<T>(0, 0)

		# Create the circular domain masks. Note that the dimension in each direction
		# is one less than the meshgrid of the total domain.
		this->cimask = this->is_shape_in_circle(center, this->ri, this->b, this->cdom)
		this->camask = this->is_shape_in_circular_aliasing_zone(center, this->ri, this->b, this->cdom)

		# Create the annular domain masks. Note that the dimension in each direction
		# is one less than the meshgrid of the total domain.
		this->aimask = this->is_shape_in_annulus(center, this->ri, this->ro, this->b, this->adom)
		this->aamask_inner, this->aamask_outer = this->is_shape_in_annulus_aliasing_zone(center, this->ri, this->ro, this->b, this->adom)

		# Arrays that hold the integrals for the circles and annuli in each cell.
		this->CI = std::vector(this->numcellsw * this->numcellsh, 0);
		this->AI = std::vector(this->numcellsw * this->numcellsh, 0);
		
		# Index grids for the circular and annular domains.
		int idx0 = -this->cspani, idxn = this->cspani + 1
		this->cind = Meshgrid<int>(circ_range, circ_range, indexing='ij')

		annulus_range = range(-this->cspano, this->cspano + 1)
		this->AII, this->AJJ = np.meshgrid(annulus_range, annulus_range, indexing='ij')
};




    #-------------------------------------------------------
    # Methods for finding extreme points on polygonal
    # shapes.
    #-------------------------------------------------------

    auto get_extreme_points_on_line(const std::pair<T> &center, const Meshgrid<T> &subdomain, char where='b'):
        /*
        get_extreme_points_on_line(c, points)

        Get the locations of the minimizer and maximimzer on a line segement
        described by the iterable of xy-coordinates points from the center c.
        Specifically, we return the critical points and the value of

        ||p - c||**2

        for p in points with c being the center of a circle.
        */

        // Placeholder for dimensions.
        int M = subdomain.numx - 1, N = subdomain.numy - 1;
		Meshgrid<T> D0, D1;

		switch (where) {
			case 'r':
				D0 = subdomain.slice_copy(1, M  , 0, N-1);
				D1 = subdomain.slice_copy(1, M  , 1, N  );
				
				break;
				
			case 't':
				D0 = subdomain.slice_copy(1, M  , 1, N  );
				D1 = subdomain.slice_copy(0, M-1, 1, N  );
				
				break;
				
			case 'l':
				D0 = subdomain.slice_copy(0, M-1, 1, N  );
				D1 = subdomain.slice_copy(0, M-1, 0, N-1);

			case 'b':
			default:
				D0 = subdomain.slice_copy(0, M-1, 0, N-1);
				D1 = subdomain.slice_copy(1, M  , 0, N-1);
		}
		
		T dx = D0.dx, dy = D0.dy;

        // Unpack the center.
        T cx = center.first, cy = center.second.

        // Function of t to minimize.
        f = lambda tau : (DX * tau + X0 - cx)**2 + (DY * tau + Y0 - cy)**2

        // Times for the critical points.
		int stride = 3;
		std::vector<T> FC(stride * D0.numx * D0.numy, 0);
		std::vector<T> FMIN(D0.numx * D0.numy, 0);
		std::vector<T> FMAX(D0.numx * D0.numy, 0);
		std::vector<T> TC(stride * D0.numx * D0.numy, 0);
		std::vector<T> deltax(D0.numx * D0.numy, 0);
		std::vector<T> deltay(D0.numx * D0.numy, 0);
		int N = stride * D0.numx * D0.numy;

		for(int k = 0; k < N; k += stride) {
			T tcritical = ((cx - D0.XX[k]) * dx + (cy - D0.YY[k]) * dy) / (dx**2 + dy**2)
			
			TC[k] = 0;
			TC[k + 1] = tcritical < 0 ? 0 : (tcricital > 1 ? 1 : tcritical);
			TC[k + 2] = 1;
			
			FMIN[k / stride] = std::min({F[k], F[k+1], F[k+2]});
			FMAX[k / stride] = std::min({F[k], F[k+1], F[k+2]});
		}
		
        // Evaluate the critical points on this line.
        fcrit = f(T)

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

        XC = np.zeros((4, M-1, N-1))
        YC = np.zeros((4, M-1, N-1))
        XF = np.zeros((4, M-1, N-1))
        YF = np.zeros((4, M-1, N-1))
        closest = np.zeros((4, M-1, N-1))
        furthest = np.zeros((4, M-1, N-1))

        sides = ('b', 'r', 't', 'l')

        for i, s in enumerate(sides):
            pc, pf = this->get_extreme_points_on_line(c, points, where=s)

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
        rsq = (r - b / 2)**2
        closest, furthest = this->get_extreme_points(c, points)

        a, b = rsq - closest[-1], furthest[-1] - rsq

        return (a > 0) & (b <= a)

    def is_shape_in_circular_aliasing_zone(self, c, r, b, points):
        rmbsq = (r - b / 2)**2
        rpbsq = (r + b / 2)**2
        closest, furthest = this->get_extreme_points(c, points)

        am, bm = rmbsq - closest[-1], furthest[-1] - rmbsq
        ap, bp = rpbsq - closest[-1], furthest[-1] - rpbsq

        exceeds_inner_radius = (am <= 0) | (bm > am)
        does_not_exceed_outer_radius = (ap > 0) & (bp <= ap)

        return exceeds_inner_radius & does_not_exceed_outer_radius

    def is_shape_outside_circle(self, c, r, b, points):
        rsq = (r + b / 2)**2
        closest, _ = this->get_extreme_points(c, points)

        a = rsq - closest[-1]

        return a < 0


    #-------------------------------------------------------
    # Methods for determining if a polygonal shape is within
    # an annulus of inner radius ri and outer radius ro.
    #-------------------------------------------------------

    def is_shape_in_annulus(self, c, ri, ro, b, points):
        risq = (ri + b / 2)**2
        rosq = (ro - b / 2)**2
        closest, furthest = this->get_extreme_points(c, points)

        am, bm = risq - closest[-1], furthest[-1] - risq
        ap, bp = rosq - closest[-1], furthest[-1] - rosq

        exceeds_inner_radius = (am <= 0) | (bm > am)
        does_not_exceed_outer_radius = (ap > 0) & (bp <= ap)

        return exceeds_inner_radius & does_not_exceed_outer_radius

    def is_shape_in_annulus_aliasing_zone(self, c, ri, ro, b, points):
        inner_aliasing_region = this->is_shape_in_annulus(c, ri - b/2, ri + b/2, 0, points)
        outer_aliasing_region = this->is_shape_in_annulus(c, ro - b/2, ro + b/2, 0, points)

        return inner_aliasing_region, outer_aliasing_region

    def is_shape_outside_annulus(self, c, ri, ro, b, points):
        risq = (ri - b / 2)**2
        rosq = (ro + b / 2)**2
        closest, furthest = this->get_extreme_points(c, points)

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
        M, N = this->F.shape

        # Get the value of the state function restricted to the circular subgrid.
        FR = this->F[np.mod(this->CII + i, M), np.mod(this->CJJ + j, N)]

        exact_area = np.pi * this->ri**2
        integral = 0
        diffarea = this->dx * this->dy

        integral = np.sum(FR[this->cimask]) * diffarea

        # Add the contributions from the aliasing region.
        if this->b > 0:
            integral += np.sum(FR[this->camask] * this->cascaling[this->camask]) * diffarea

        return integral / exact_area


    def compute_annular_integral(self, i, j):
        # Dimensions of the state function.
        M, N = this->F.shape

        # Get the value of the state function restricted to the circular subgrid.
        FR = this->F[np.mod(this->AII + i, M), np.mod(this->AJJ + j, N)]

        exact_area = np.pi * (this->ro**2 - this->ri**2)
        diffarea = this->dx * this->dy
        integral = np.sum(FR[this->aimask]) * diffarea

        if this->b > 0:
            # Add the contribution from the inner aliasing region.
            integral += np.sum(FR[this->aamask_inner] * this->aiscaling[this->aamask_inner]) * diffarea

            # Add the contribution from the outer aliasing region.
            integral += np.sum(FR[this->aamask_outer] * this->aoscaling[this->aamask_outer]) * diffarea

        return integral / exact_area

    def compute_integrals(self, serial=True, nprocs=multiprocessing.cpu_count()):
        M, N = this->numcellsw, this->numcellsh

        if serial:
            for i in range(M):
                for j in range(N):
                    this->CI[i, j] = this->compute_circular_integral(i, j)
                    this->AI[i, j] = this->compute_annular_integral(i, j)
        else:
            ci_cartprod = itertools.product(range(M), range(N))
            ai_cartprod = itertools.product(range(M), range(N))

            with multiprocessing.Pool(processes=nprocs) as p:
                ci_results = p.starmap_async(this->compute_circular_integral, ci_cartprod)
                ai_results = p.starmap_async(this->compute_annular_integral, ai_cartprod)

                this->CI[:,:] = np.reshape(ci_results.get(), (M, N))
                this->AI[:,:] = np.reshape(ai_results.get(), (M, N))

    #-------------------------------------------------------
    # Methods for testing out drawing and masks.
    #-------------------------------------------------------

    def draw_circular_region(self):
        M, N = this->XX.shape
        XX, YY = this->XX[:M-1,:N-1] + this->dx/2, this->YY[:M-1,:N-1] + this->dy/2
        #XX, YY = this->cdom
        x = XX[:, 0]
        y = YY[0, :]
        x0, xn = x[0], x[-1]
        y0, yn = y[0], y[-1]
        numx = len(x)
        numy = len(y)

        # Center of the circle.
        c = (0, 0)

        # Inner radius for the circle.
        ri = this->ri

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
        radii = [this->ri - this->b/2, this->ri, this->ri + this->b/2]
        linestyles = ['-', '-', '-']

        for clr, rad, ls in zip(colors, radii, linestyles):
            ax.plot(c[0] + rad * np.cos(t), c[1] + rad * np.sin(t), color=clr, linestyle=ls)


        inmask = this->is_shape_in_circle(c, this->ri, this->b, (XX, YY))
        camask = this->is_shape_in_circular_aliasing_zone(c, this->ri, this->b, (XX, YY))

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
        i = int(np.round((c[0] - x0) / this->dx))
        j = int(np.round((c[1] - y0) / this->dy))
        approx_average = this->compute_circular_integral(i, j)

        title_str = 'Circular region, $r_i \\approx$ {0:.3f}'.format(ri)
        title_str += '\nApproximate average: {0:.10f}'.format(approx_average)
        title_str += '\nRelative error: {0:.4f}%'.format((approx_average - 1) * 100)
        ax.set_title(title_str)
        ax.set_xlim((x0, xn))
        ax.set_ylim((y0, yn))
        plt.show()


    def draw_annular_region(self):
        M, N = this->XX.shape
        XX, YY = this->XX[:M-1,:N-1] + this->dx/2, this->YY[:M-1,:N-1] + this->dy/2
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
        radii = [this->ri - this->b/2, this->ri, this->ri + this->b/2,
                 this->ro - this->b/2, this->ro, this->ro + this->b/2]
        linestyles = ['-'] * len(radii)

        for clr, rad, ls in zip(colors, radii, linestyles):
            ax.plot(c[0] + rad * np.cos(t), c[1] + rad * np.sin(t), color=clr, linestyle=ls)

        inmask = this->is_shape_in_annulus(c, this->ri, this->ro, this->b, (XX, YY))
        iamask, oamask = this->is_shape_in_annulus_aliasing_zone(c, this->ri, this->ro, this->b, (XX, YY))

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
        i = int(np.round((c[0] - x0) / this->dx))
        j = int(np.round((c[1] - y0) / this->dy))
        approx_average = this->compute_annular_integral(i, j)

        title_str = 'Annular region, $r_i \\approx$ {0:.3f}, $r_o \\approx$ {1:.3f}'.format(this->ri, this->ro)
        title_str += '\nApproximate average: {0:.10f}'.format(approx_average)
        title_str += '\nRelative error: {0:.4f}%'.format((approx_average - 1) * 100)
        ax.set_title(title_str)
        ax.set_xlim((x0, xn))
        ax.set_ylim((y0, yn))
        plt.show()



