#include "grid.h"
#include "transitionfunction.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#ifdef USE_INTRINSICS
#include "intrinsics_utils.h"
#include <immintrin.h>
#endif

#define FORMAT_BUFFER_SIZE 256
#define NUM_REGIONS 5


int grid_init(struct grid_s *grid, float x0, float xn, int numx, float y0, float yn, int numy, float ri_ratio, float b)
{
	float deltax, deltay, mindelta, mindifferential;
	int init_status = 1;
	
	grid->x0 = x0;
	grid->xn = xn;
	
	grid->y0 = y0;
	grid->yn = yn;
	
	deltax = grid->xn - grid->x0;
	deltay = grid->yn - grid->y0;

	grid->numx = numx;
	grid->numy = numy;
	
	grid->dx = deltax / (grid->numx - 1);
	grid->dy = deltay / (grid->numy - 1);
	grid->dA = grid->dx * grid->dy;

	mindelta = deltax < deltay ? deltax : deltay;
	mindifferential = grid->dx < grid->dy ? grid->dx : grid->dy;
	
	grid->ri = ri_ratio * mindelta;
	grid->ri_ratio = ri_ratio;
	grid->ro = 3 * grid->ri;
	grid->b = b;
	
	grid->cspani = (int)ceil(grid->ri / mindifferential);
	grid->cspano = (int)ceil(grid->ro / mindifferential);
	
	grid->cdiami = 2 * grid->cspani + 1;
	grid->cdiamo = 2 * grid->cspano + 1;
	
	grid->cstart = 0;
	grid->cend = grid->cdiami * grid->cdiami;
	grid->astart = grid->cend;
	grid->aend = grid->astart + grid->cdiamo * grid->cdiamo;
	grid->fstart = grid->aend;
	grid->fend = grid->fstart + (grid->numx - 1) * (grid->numy - 1);

	grid->circle_area = M_PI * grid->ri * grid->ri;
	grid->annulus_area = M_PI * (grid->ro * grid->ro - grid->ri * grid->ri);

#ifdef DEBUG
	printf("DEBUG: Allocating memory\n");
#endif

	grid->data = calloc(grid->fend, sizeof(float));
	grid->distsq_closest = calloc(grid->fstart, sizeof(float));
	grid->distsq_furthest = calloc(grid->fstart, sizeof(float));
	grid->dist_midpoint = calloc(grid->fstart, sizeof(float));
	grid->weights = calloc(grid->fstart, sizeof(float));
	grid->domain_indices = calloc(2 * (grid->cdiami + grid->cdiamo), sizeof(int));
	
	if(grid->data == NULL || grid->distsq_closest == NULL || grid->distsq_furthest == NULL ||
	   grid->weights == NULL || grid->dist_midpoint == NULL || grid->domain_indices == NULL) {

		init_status = 0;

		free(grid->data);
		free(grid->distsq_closest);
		free(grid->distsq_furthest);
		free(grid->weights);
		free(grid->domain_indices);
		
		grid->data = NULL;
		grid->distsq_closest = NULL;
		grid->distsq_furthest = NULL;
		grid->weights = NULL;
		grid->domain_indices = NULL;
	} else {
#ifdef DEBUG
		printf("DEBUG: Allocation successful.\n");
		grid_print(grid);
#endif

#ifdef DEBUG
		printf("DEBUG: Setting indices.\n");
#endif

		grid_set_indices(grid, 0, 0);

#ifdef DEBUG
		printf("DEBUG: Computing midpoint distances.\n");
#endif

		grid_compute_midpoint_distances(grid, 0, 0);

#ifdef DEBUG
		printf("DEBUG: Computing extreme points.\n");
#endif

		grid_get_extreme_points(grid);

#ifdef DEBUG
		printf("DEBUG: Computing weights.\n");
#endif

		grid_compute_weights(grid);
		
#ifdef DEBUG
		printf("DEBUG: Setting transition function state.\n");
#endif

		tfunc_state_init(&grid->tfs);
	}
	
	return init_status;
}


void grid_set_indices(struct grid_s *grid, int center_index_i, int center_index_j)
{
	int *indicesi = grid->domain_indices;
	int *indicesj = grid->domain_indices + grid->cdiami;
	int k;
	
#ifdef USE_INTRINSICS
	__m256i istartreg = _mm256_set1_epi32(-grid->cspani + center_index_i);
	__m256i jstartreg = _mm256_set1_epi32(-grid->cspani + center_index_j);
	__m256i seqreg = _mm256_setr_epi32(0, 1, 2, 3, 4, 5, 6, 7);
	__m256i ireg, jreg, kreg;
	__m256i allmask = _mm256_set_mask_epi32(0);
	__m256i cutoffmask;
	int cutoff = grid->cdiami % INT32_PER_M256_REG;
	
	if (cutoff > 0) {
		cutoffmask = _mm256_set_mask_epi32(cutoff);
		
		kreg = _mm256_set1_epi32(k);
		kreg = _mm256_add_epi32(kreg, seqreg);
		ireg = _mm256_add_epi32(istartreg, kreg);
		jreg = _mm256_add_epi32(jstartreg, kreg);
		
		_mm256_maskstore_epi32(indicesi, cutoffmask, ireg);
		_mm256_maskstore_epi32(indicesj, cutoffmask, jreg);
	}

	for (k = cutoff; k < grid->cdiami; k += INT32_PER_M128_REG) {
		kreg = _mm256_set1_epi32(k);
		kreg = _mm256_add_epi32(kreg, seqreg);
		ireg = _mm256_add_epi32(istartreg, kreg);
		jreg = _mm256_add_epi32(jstartreg, kreg);

		_mm256_maskstore_epi32(indicesi + k, allmask, ireg);
		_mm256_maskstore_epi32(indicesj + k, allmask, jreg);
	}
	
	indicesi = grid->domain_indices + 2*grid->cdiami;
	indicesj = indicesi + grid->cdiamo;
	
	cutoff = grid->cdiami % INT32_PER_M256_REG;
	
	if (cutoff > 0) {
		cutoffmask = _mm256_set_mask_epi32(cutoff);
		
		kreg = _mm256_set1_epi32(k);
		kreg = _mm256_add_epi32(kreg, seqreg);
		ireg = _mm256_add_epi32(istartreg, kreg);
		jreg = _mm256_add_epi32(jstartreg, kreg);
		
		_mm256_maskstore_epi32(indicesi, cutoffmask, ireg);
		_mm256_maskstore_epi32(indicesj, cutoffmask, jreg);
	}
	
	for (k = cutoff; k < grid->cdiamo; k += INT32_PER_M128_REG) {
		kreg = _mm256_set1_epi32(k);
		kreg = _mm256_add_epi32(kreg, seqreg);
		ireg = _mm256_add_epi32(istartreg, kreg);
		jreg = _mm256_add_epi32(jstartreg, kreg);

		_mm256_maskstore_epi32(indicesi + k, allmask, ireg);
		_mm256_maskstore_epi32(indicesj + k, allmask, jreg);
	}

#else

	for (k = 0; k < grid->cdiami; k++) {
		indicesi[k] = -grid->cspani + k + center_index_i;
		indicesj[k] = -grid->cspani + k + center_index_j;
	}
	
	indicesi = grid->domain_indices + 2 * grid->cdiami;
	indicesj = indicesi + grid->cdiamo;

	for (k = 0; k < grid->cdiamo; k++) {
		indicesi[k] = -grid->cspano + k + center_index_i;
		indicesj[k] = -grid->cspano + k + center_index_j;
	}
	
#endif
}

void grid_compute_grid_index_modulus(struct grid_s *grid)
{
	int *indicesi = grid->domain_indices;
	int *indicesj = grid->domain_indices + grid->cdiami;
	int M = grid->numx - 1;
	int N = grid->numy - 1;
	int k;
	
	// Set the indices for the grid.
	for (k = 0; k < grid->cdiami; k++) {
		indicesi[k] = indicesi[k] % M;
		indicesj[k] = indicesj[k] % N;
		
		if (indicesi[k] < 0) {
			indicesi[k] += M;
		}
		
		if (indicesj[k] < 0) {
			indicesj[k] += N;
		}
	}
	
	indicesi = grid->domain_indices + 2*grid->cdiami;
	indicesj = indicesi + grid->cdiamo;

	for (k = 0; k < grid->cdiamo; k++) {
		indicesi[k] = indicesi[k] % M;
		indicesj[k] = indicesj[k] % N;
		
		if (indicesi[k] < 0) {
			indicesi[k] += M;
		}
		
		if (indicesj[k] < 0) {
			indicesj[k] += N;
		}
	}
}



void grid_compute_midpoint_distances(struct grid_s *grid, float cx, float cy)
{
	float *dist = grid->dist_midpoint + grid->cstart;
	int *indicesi = grid->domain_indices;
	int *indicesj = grid->domain_indices + grid->cdiami;
	int M, N;
	int i, j, jind, k = 0;
	int istart, iend, jstart, jend;
	float xm, ym;
	float distx2, disty2;
	float x0, y0;
	
	for (j = 0 ; j < grid->cdiami; j++) {
		// ym = grid->y0 + 0.5 * grid->dy * indicesj[j];
		ym = indicesj[j] * grid->dy;
		disty2 = (ym - cy) * (ym - cy);
		jind = j * grid->cdiami;

		for (i = 0; i < grid->cdiami; i++) {
			// xm = grid->x0 + 0.5 * grid->dx * indicesi[i];
			xm = indicesi[i] * grid->dx;
			distx2 = (xm - cx) * (xm - cx);

			dist[jind + i] = sqrt(distx2 + disty2);
			
// #ifdef DEBUG
			// printf("DEBUG: Circular domain -> indicesi[%2d] = %+2d, indicesj[%2d] = %+2d, xm = %+f, ym = %+f, dist[%3d, %3d] = dist[%3d] = %f\n", i, indicesi[i], j, indicesj[j], xm, ym, i, j, jind + i, dist[jind + i]);
// #endif
		}
	}
	
// #ifdef DEBUG
	// printf("\n");
// #endif
	
	dist = grid->dist_midpoint + grid->astart;
	indicesi = grid->domain_indices + 2 * grid->cdiami;
	indicesj = indicesi + grid->cdiamo;
	
	for (j = 0 ; j < grid->cdiamo; j++) {
		// ym = grid->y0 + 0.5 * grid->dy * indicesj[j];
		ym = indicesj[j] * grid->dy;
		disty2 = (ym - cy) * (ym - cy);
		jind = j * grid->cdiamo;

		for (i = 0; i < grid->cdiamo; i++) {
			// xm = grid->x0 + 0.5 * grid->dx * indicesi[i];
			xm = indicesi[i] * grid->dx;
			distx2 = (xm - cx) * (xm - cx);

			dist[jind + i] = sqrt(distx2 + disty2);
// #ifdef DEBUG
			// printf("DEBUG: Annular domain -> indicesi[%2d] = %+2d, indicesj[%2d] = %+2d, xm = %+f, ym = %+f, dist[%3d, %3d] = dist[%3d] = %f\n", i, indicesi[i], j, indicesj[j], xm, ym, i, j, jind + i, dist[jind + i]);
// #endif
		}
	}
}


void grid_get_extreme_points_on_line(struct grid_s *grid, float cx, float cy, char where, int i, int j, float *xcrit, float *ycrit, float *dstsq)
{
	/*
	get_extreme_points_on_line(c, points)

	Get the locations of the minimizer and maximimzer on a line segement
	described by the iterable of xy-coordinates points from the center c.
	Specifically, we return the critical points and the value of

	||p - c||**2

	for p in points with c being the center of a circle.
	*/

	// Compute the incremental differences in the coordinates.
	float dx = grid->dx;
	float dy = grid->dy;
	
	float tc[3];
	float distvals[3];
	
	int k;
	int minidx, maxidx;
	int M = grid->numx - 1;
	int N = grid->numy - 1;
	
	float x0, x1, y0, y1;
	float differential_sum;
	float xdist, ydist;
	
	if (where == 'b') {
		dx = grid->dx;
		dy = 0;
		
		// x0 = grid->x0 + i * grid->dx;
		// y0 = grid->y0 + j * grid->dy;
		
		x0 = i * grid->dx - 0.5 * grid->dx;
		y0 = j * grid->dy - 0.5 * grid->dy;
	} else if (where == 't') {
		dx = -grid->dx;
		dy = 0;
		
		// x0 = grid->x0 + (i + 1) * grid->dx;
		// y0 = grid->y0 + (j + 1) * grid->dy;
		
		x0 = (i + 1) * grid->dx - 0.5 * grid->dx;
		y0 = (j + 1) * grid->dy - 0.5 * grid->dy;
	} else if (where == 'r') {
		dx = 0;
		dy = grid->dy;
		
		// x0 = grid->x0 + (i + 1) * grid->dx;
		// y0 = grid->y0 + j * grid->dy;
		
		x0 = (i + 1) * grid->dx - 0.5 * grid->dx;
		y0 = j * grid->dy - 0.5 * grid->dy;
	} else if (where == 'l') {
		dx = 0;
		dy = -grid->dy;
		
		// x0 = grid->x0 + i * grid->dx;
		// y0 = grid->y0 + (j + 1) * grid->dy;
		
		x0 = i * grid->dx - 0.5 * grid->dx;
		y0 = (j + 1) * grid->dy - 0.5 * grid->dy;
	}
	
	differential_sum = dx * dx + dy * dy;

	// The critical values for time are at the endpoints.
	tc[0] = 0;
	tc[1] = ((cx - x0) * dx + (cy - y0) * dy) / differential_sum;
	tc[2] = 1;
	
	if (tc[1] < 0) {
		tc[1] = 0;
	} else if (tc[1] > 1) {
		tc[1] = 1;
	}
	
	// Compute the squared distances.
	for (k = 0; k < 3; k++) {
		xdist = (dx * tc[k] + x0 - cx);
		ydist = (dy * tc[k] + y0 - cy);

		distvals[k] = xdist * xdist + ydist * ydist;
	}

	grid_get_extreme_indices(distvals, 3, &minidx, &maxidx);
	
	xcrit[0] = dx * tc[minidx] + x0;
	ycrit[0] = dy * tc[minidx] + y0;
	dstsq[0] = distvals[minidx];

	xcrit[1] = dx * tc[maxidx] + x0;
	ycrit[1] = dy * tc[maxidx] + y0;
	dstsq[1] = distvals[maxidx];
}

void grid_get_extreme_indices(float *values, int length, int *minindex, int *maxindex)
{
	int minidx, maxidx;
	float minval, maxval;
	
	minidx = maxidx = 0;
	minval = maxval = values[0];
	
	for (int i = 1; i < length; i++) {
		if (minval > values[i]) {
			minidx = i;
			minval = values[i];
		}
		
		if (maxval <  values[i]) {
			maxidx = i;
			maxval = values[i];
		}
	}
	
	*minindex = minidx;
	*maxindex = maxidx;
}

void grid_get_extreme_points(struct grid_s *grid)
{
	/*
	get_extreme_points_on_shape(c, points)

	Get the locations of the minimizer and maximimzer on a (polygonal)
	shape described by the iterable of xy-coordinates points from the
	center c. Specifically, we return the critical points and the value of

	||p - c||**2

	for p in points with c being the center of a circle.
	*/
	
	int i, j, k = 0;
	int di, dj;
	float cx = 0, cy = 0;
	float xcrit[8];
	float ycrit[8];
	float dstsq[8];
	int minidx, maxidx;
	int *indicesi = grid->domain_indices;
	int *indicesj = grid->domain_indices + grid->cdiami;

	for (j = 0; j < grid->cdiami; j++) {
		dj = indicesj[j];
		
		for (i = 0; i < grid->cdiami; i++) {
			di = indicesi[i];

			grid_get_extreme_points_on_line(grid, cx, cy, 'b', di, dj, xcrit, ycrit, dstsq);
			grid_get_extreme_points_on_line(grid, cx, cy, 't', di, dj, xcrit + 2, ycrit + 2, dstsq + 2);
			grid_get_extreme_points_on_line(grid, cx, cy, 'l', di, dj, xcrit + 4, ycrit + 4, dstsq + 4);
			grid_get_extreme_points_on_line(grid, cx, cy, 'r', di, dj, xcrit + 6, ycrit + 6, dstsq + 6);

			
			grid_get_extreme_indices(dstsq, 8, &minidx, &maxidx);
			grid->distsq_closest[k] = dstsq[minidx];
			grid->distsq_furthest[k] = dstsq[maxidx];

			k++;
		}
	}

	k = grid->astart;
	indicesi = grid->domain_indices + 2*grid->cdiami;
	indicesj = indicesi + grid->cdiamo;

	for (j = 0; j < grid->cdiamo; j++) {
		dj = indicesj[j];
		
		for (i = 0; i < grid->cdiamo; i++) {
			di = indicesi[i];

			grid_get_extreme_points_on_line(grid, cx, cy, 'b', di, dj, xcrit, ycrit, dstsq);
			grid_get_extreme_points_on_line(grid, cx, cy, 't', di, dj, xcrit + 2, ycrit + 2, dstsq + 2);
			grid_get_extreme_points_on_line(grid, cx, cy, 'l', di, dj, xcrit + 4, ycrit + 4, dstsq + 4);
			grid_get_extreme_points_on_line(grid, cx, cy, 'r', di, dj, xcrit + 6, ycrit + 6, dstsq + 6);
			
			grid_get_extreme_indices(dstsq, 8, &minidx, &maxidx);
			grid->distsq_closest[k] = dstsq[minidx];
			grid->distsq_furthest[k] = dstsq[maxidx];
			
			k++;
		}
	}
}

void grid_compute_weights(struct grid_s *grid)
{
	int k;
	float *weights = grid->weights;
	float *mdist = grid->dist_midpoint;
	float ri = grid->ri;
	float ro = grid->ro;
	float b = grid->b;
	
	memset(weights, 0, grid->aend * sizeof(float));

	for (k = grid->cstart; k < grid->cend; k++) {
		if (grid_is_cell_inside_circle(grid, k, ri, b)) {
			weights[k] = 1;
		} else if (grid_is_cell_inside_circle_aliasing_zone(grid, k, ri, b)) {
			weights[k] = (ri + 0.5 * b - mdist[k]) / b;
		}
	}
	
	for (k = grid->astart; k < grid->aend; k++) {
		if (grid_is_cell_in_annulus(grid, k, ri, ro, b)) {
			weights[k] = 1;
		} else if (grid_is_shape_in_inner_annulus_aliasing_zone(grid, k, ri, ro, b)) {
			weights[k] = (mdist[k] - (ri - 0.5 * b)) / b;
		} else if (grid_is_shape_in_outer_annulus_aliasing_zone(grid, k, ri, ro, b)) {
			weights[k] = (ro + 0.5 * b - mdist[k]) / b;
		}
	}
}

void grid_compute_inner_circle_weights_avx2(struct grid_s *grid)
{
	int k;
	float *weights = grid->weights;
	float *mdist = grid->dist_midpoint;
	float *cdist = grid->distsq_closest;
	float *fdist = grid->distsq_furthest;
	float ri = grid->ri;
	float ro = grid->ro;
	float b = grid->b;
	
	__m256 hbreg = _mm256_set1_ps(0.5f * b);
	__m256 ones = _mm256_set1_ps(1.0f);
	__m256 zero = _mm256_set1_ps(0.0f);
	
	__m256 areg, breg;
	__m256 cdreg;
	__m256 fdreg;
	__m256 mdreg;
	
	// Weight register.
	__m256 wreg;
	
	__m256 rmask[3];  // Region masks
	
	__m256 rsq = _mm256_set1_ps(ri);
	rsq = _mm256_sub_ps(rsq, hbreg);
	rsq = _mm256_mul_ps(rsq, rsq);
	
	if (cutoff > 0) {
		lsmask = _mm256_set_mask_epi32(cutoff);
		
		// Load the closest, furthest, and midpoint distances.
		cdreg = _mm256_maskload_ps(cdist, lsmask);
		fdreg = _mm256_maskload_ps(fdist, lsmask);
		mdreg = _mm256_maskload_ps(mdist, lsmask); // 8
		
		// Set the masks.
		areg = _mm256_sub_ps(rsq, cdist);
		breg = _mm256_sub_ps(fdist, rsq); // 4
		rmask[0] = _mm256_cmp_ps(areg, zero, _CMP_GT_OS);
		rmask[1] = _mm256_cmp_ps(breg, areg, _CMP_LT_OS); // 4
		
		// Compute the region mask.
		rmask[2] = _mm256_and_ps(rmask[0], rmask[1]); // 1
		
		// Set the weights.
		// (1.0 / b) * (ri - mdist) + 0.5
		wreg = _mm256_and_ps(rmask[2], ones); // 1
		
		// Store the weights.
		_mm256_store_ps(weights, wreg); // 5
		
		// Store the weights.
		_mm256_maskstore_ps(weights, lsmask, rmask[0]); // 5
	}

	for (k = grid->cstart; k < grid->cend; k += FLOAT_PER_M256_REG) {
		// Load the closest, furthest, and midpoint distances.
		cdreg = _mm256_load_ps(cdist + k);
		fdreg = _mm256_load_ps(fdist + k);
		mdreg = _mm256_load_ps(mdist + k); // 6

		// Set the masks.
		areg = _mm256_sub_ps(rsq, cdist);
		breg = _mm256_sub_ps(fdist, rsq); // 4
		rmask[0] = _mm256_cmp_ps(areg, zero, _CMP_GT_OS);
		rmask[1] = _mm256_cmp_ps(breg, areg, _CMP_LT_OS); // 4
		
		// Compute the region mask.
		rmask[2] = _mm256_and_ps(rmask[0], rmask[1]); // 1
		
		// Set the weights.
		// (1.0 / b) * (ri - mdist) + 0.5
		wreg = _mm256_and_ps(rmask[2], ones); // 1
		
		// Store the weights.
		_mm256_store_ps(weights, wreg); // 5
	}
}

void grid_compute_aliasing_circle_weights_avx2(struct grid_s *grid)
{
	int k;
	float *weights = grid->weights;
	float *mdist = grid->dist_midpoint;
	float *cdist = grid->distsq_closest;
	float *fdist = grid->distsq_furthest;
	float ri = grid->ri;
	float ro = grid->ro;
	float b = grid->b;
	int cutoff = (grid->cend - grid->cstart) % FLOAT_PER_M256_REG;
	
	__m256 hreg = _mm256_set1_ps(0.5f);
	__m256 hbreg = _mm256_set1_ps(0.5f * b);
	__m256 brreg = _mm256_set1_ps(1.0f / b);
	__m256 zero = _mm256_set1_ps(0.0f);
	
	__m256 areg, breg;
	__m256 cdreg;
	__m256 fdreg;
	__m256 mdreg;
	
	__m256 lsmask;  // Load-store mask.
	__m256 rmask[2];  // Region masks
	
	__m256 rireg = _mm256_set1_ps(ri);
	__m256 roreg = _mm256_set1_ps(ro);
	__m256 wreg;
	__m256 rsq[2];
	
	rsq[0] = _mm256_sub_ps(rireg, hbreg);
	rsq[1] = _mm256_add_ps(rireg, hbreg);
	
	rsq[0] = _mm256_mul_ps(rsq[0], rsq[0]);
	rsq[1] = _mm256_mul_ps(rsq[1], rsq[1]);
	
	
	// (ri + b/2 - mdist) / b = (1/b)*ri + 0.5 - (1/b)*mdist
	wreg = _mm256_fmadd_ps(brreg, rireg, hreg);
	
	if (cutoff > 0) {
		lsmask = _mm256_set_mask_epi32(cutoff);
		
		// Load the closest, furthest, and midpoint distances.
		cdreg = _mm256_maskload_ps(cdist, lsmask);
		fdreg = _mm256_maskload_ps(fdist, lsmask);
		mdreg = _mm256_maskload_ps(mdist, lsmask); // 8

		// Set the masks for the anti-aliasing region.
		areg = _mm256_sub_ps(rsq[0], cdreg);
		breg = _mm256_sub_ps(fdreg, rsq[0]); // 4
		rmask[0] = _mm256_cmp_ps(areg, zero, _CMP_LE_OS);
		rmask[1] = _mm256_cmp_ps(breg, areg, _CMP_GT_OS); // 4
		rmask[0] = _mm256_or_ps(rmask[0], rmask[1]); // 1
		rmask[0] = _mm256_and_ps(rmask[0], _mm256_castsi256_ps(lsmask)); // 1
		
		areg = _mm256_sub_ps(rsq[1], cdreg);
		breg = _mm256_sub_ps(fdreg, rsq[1]); // 4
		rmask[0] = _mm256_cmp_ps(areg, zero, _CMP_GT_OS);
		rmask[1] = _mm256_cmp_ps(breg, areg, _CMP_LE_OS); // 4
		rmask[1] = _mm256_and_ps(rmask[0], rmask[1]); // 1
		rmask[1] = _mm256_and_ps(rmask[1], _mm256_castsi256_ps(lsmask)); // 1
		
		// The anti-aliasing mask proper.
		rmask[1] = _mm256_and_ps(rmask[0], rmask[1]); // 1
		
		// Set the weights.
		// (1.0 / b) * (ri - mdist) + 0.5
		rmask[0] = _mm256_mul_ps(brreg, mdist); // 4
		rmask[0] = _mm256_sub_ps(wreg, rmask[0]); // 4
		rmask[0] = _mm256_and_ps(rmask[0], rmask[1]); // 1
		
		// Store the weights.
		_mm256_maskstore_ps(weights, lsmask, rmask[0]); // 5
	}
	
	for (k = grid->cstart; k < grid->cend; k += FLOAT_PER_M256_REG) {
		// Load the closest, furthest, and midpoint distances.
		cdreg = _mm256_load_ps(cdist + k);
		fdreg = _mm256_load_ps(fdist + k);
		mdreg = _mm256_load_ps(mdist + k); // 6

		// Set the masks for the anti-aliasing region.
		areg = _mm256_sub_ps(rsq[0], cdist);
		breg = _mm256_sub_ps(fdist, rsq[0]); // 4
		rmask[0] = _mm256_cmp_ps(areg, zero, _CMP_LE_OS);
		rmask[1] = _mm256_cmp_ps(breg, areg, _CMP_GT_OS); // 4
		rmask[0] = _mm256_or_ps(rmask[0], rmask[1]); // 1
		
		areg = _mm256_sub_ps(rsq[1], cdist);
		breg = _mm256_sub_ps(fdist, rsq[1]); // 4
		rmask[0] = _mm256_cmp_ps(areg, zero, _CMP_GT_OS);
		rmask[1] = _mm256_cmp_ps(breg, areg, _CMP_LE_OS); // 4
		rmask[1] = _mm256_and_ps(rmask[0], rmask[1]); // 1
		
		// The anti-aliasing mask proper.
		rmask[1] = _mm256_and_ps(rmask[0], rmask[1]); // 1
		
		// Set the weights.
		// (1.0 / b) * (ri - mdist) + 0.5
		rmask[0] = _mm256_mul_ps(brreg, mdist); // 4
		rmask[0] = _mm256_sub_ps(wreg, rmask[0]); // 4
		rmask[0] = _mm256_and_ps(rmask[0], rmask[1]); // 1
		
		// Store the weights.
		_mm256_store_ps(weights + k, wreg); // 5
	}
}

int grid_is_cell_inside_circle(const struct grid_s *grid, int k, float ri, float b)
{
	float rsq = (ri - b / 2) * (ri - b / 2);
	float am = rsq - grid->distsq_closest[k];
	float bm = grid->distsq_furthest[k] - rsq;

	return (am > 0) && (bm <= am);
}

int grid_is_cell_inside_circle_aliasing_zone(const struct grid_s *grid, int k, float ri, float b)
{
	return grid_is_cell_in_annulus(grid, k, ri - b/2, ri + b/2, 0);
}


int grid_is_cell_in_annulus(const struct grid_s *grid, int k, float ri, float ro, float b) {
	float risq = (ri + b / 2) * (ri + b / 2);
	float rosq = (ro - b / 2) * (ro - b / 2);

	float am = risq - grid->distsq_closest[k];
	float bm = grid->distsq_furthest[k] - risq;
	float ap = rosq - grid->distsq_closest[k];
	float bp = grid->distsq_furthest[k] - rosq;

	int exceeds_inner_radius = (am <= 0) || (bm > am);
	int does_not_exceed_outer_radius = (ap > 0) && (bp <= ap);

	return exceeds_inner_radius && does_not_exceed_outer_radius;
}

int grid_is_shape_in_inner_annulus_aliasing_zone(const struct grid_s *grid, int k, float ri, float ro, float b)
{
	return grid_is_cell_in_annulus(grid, k, ri - b/2, ri + b/2, 0);
}

int grid_is_shape_in_outer_annulus_aliasing_zone(const struct grid_s *grid, int k, float ri, float ro, float b)
{
	return grid_is_cell_in_annulus(grid, k, ro - b/2, ro + b/2, 0);
}


void grid_clear(struct grid_s *grid)
{
	free(grid->data);
	free(grid->distsq_closest);
	free(grid->distsq_furthest);
	free(grid->dist_midpoint);
	free(grid->weights);
	free(grid->domain_indices);
	
	grid->data = NULL;
	grid->distsq_closest = NULL;
	grid->distsq_furthest = NULL;
	grid->weights = NULL;
	grid->domain_indices = NULL;
}

void grid_set_random(struct grid_s *grid)
{
	int k;
	float *data = grid->data;

	for (k = grid->fstart; k < grid->fend; k++) {
		data[k] = (float)rand() / RAND_MAX;
	}
}


void grid_set_value(struct grid_s *grid, float value)
{
	int k;
	float *data = grid->data;
	
#ifdef USE_INTRINSICS
	k = grid->fend - grid->fstart;
	fset_avx2(data + grid->fstart, k, value);
#else
	for (k = grid->fstart; k < grid->fend; k++) {
		data[k] = value;
	}
#endif
}

void grid_compute_integral(struct grid_s *grid, char which, float *integral, float *average)
{
	int i, n;
	float sum;
	float *data, *weights;

#ifdef USE_INTRINSICS
	if (which == 'i') {
		data = grid->data + grid->cstart;
		weights = grid->weights + grid->cstart;
		n = grid->cend - grid->cstart;
	} else {
		data = grid->data + grid->astart;
		weights = grid->weights + grid->astart;
		n = grid->aend - grid->astart;
	}
		
	sum = fdot_avx2(data, weights, n);
#else
	sum = 0;
	*data = grid->data;
	*weights = grid->weights;

	if (which == 'i') {
		for (i = grid->cstart; i < grid->cend; i++) {
			sum += data[i] * weights[i];
		}
	} else {
		for (i = grid->astart; i < grid->aend; i++) {
			sum += data[i] * weights[i];
		}
	}
#endif

	*integral = sum * grid->dA;
	*average = *integral / ((which == 'i') ? grid->circle_area : grid->annulus_area);
}

void grid_center_state_data_at(struct grid_s *grid, int icenter, int jcenter)
{
	int *indicesi = grid->domain_indices;
	int *indicesj = grid->domain_indices + grid->cdiami;
	float *circular_data = grid->data + grid->cstart;
	float *annular_data = grid->data + grid->astart;
	float *state_data = grid->data + grid->fstart;
	int M = grid->numx - 1;
	int N = grid->numy - 1;
	int i, j;
	int jind, cind, aind;
	
	grid_set_indices(grid, icenter, jcenter);
	grid_compute_grid_index_modulus(grid);
	
	for (j = 0; j < grid->cdiami; j++) {
		jind = indicesj[j] * M;
		cind = j * grid->cdiami;
		
		for (i = 0; i < grid->cdiami; i++) {
			circular_data[cind + i] = state_data[jind + indicesi[i]];
		}
	}

	indicesi = grid->domain_indices + 2 * grid->cdiami;
	indicesj = indicesi + grid->cdiamo;
	
	for (j = 0; j < grid->cdiamo; j++) {
		jind = indicesj[j] * M;
		aind = j * grid->cdiamo;
		
		for (i = 0; i < grid->cdiamo; i++) {
			annular_data[aind + i] = state_data[jind + indicesi[i]];
		}
	}
}

void grid_update_state_function(struct grid_s *grid)
{
	float cint, cmean;
	float aint, amean;
	int i, j, jind;
	int M = grid->numy - 1;
	int N = grid->numx - 1;
	float *data = grid->data + grid->fstart;

	for (j = 0; j < N; j++) {
		jind = j * M;

		for (i = 0; i < M; i++) {
			grid_center_state_data_at(grid, i, j);
			
			grid_compute_integral(grid, 'i', &cint, &cmean);
			grid_compute_integral(grid, 'o', &aint, &amean);
			
			data[jind + i] = tfunc_apply(&grid->tfs, amean, cmean);
		}
	}
}


void grid_print(const struct grid_s *grid)
{
	printf("data            = %X\n", grid->data);
	printf("distsq_closest  = %X\n", grid->distsq_closest);
	printf("distsq_furthest = %X\n", grid->distsq_furthest);
	printf("dist_midpoint   = %X\n", grid->dist_midpoint);
	printf("weights         = %X\n", grid->weights);
	printf("domain_indices  = %X\n\n", grid->domain_indices);
	
	printf("x0              = %f\n", grid->x0);
	printf("xn              = %f\n", grid->xn);
	printf("numx            = %d\n", grid->numx);
	printf("dx              = %f\n\n", grid->dx);
	
	printf("y0              = %f\n", grid->y0);
	printf("yn              = %f\n", grid->yn);
	printf("numy            = %d\n", grid->numy);
	printf("dy              = %f\n\n", grid->dy);
	
	printf("dA              = %f\n", grid->dA);
	printf("ri              = %f\n", grid->ri);
	printf("ri_ratio        = %f\n", grid->ri_ratio);
	printf("ro              = %f\n", grid->ro);
	printf("b               = %f\n\n", grid->b);
	
	printf("cspani          = %d\n", grid->cspani);
	printf("cspano          = %d\n", grid->cspano);
	printf("cdiami          = %d\n", grid->cdiami);
	printf("cdiamo          = %d\n\n", grid->cdiamo);
	
	printf("cstart          = %d\n", grid->cstart);
	printf("cend            = %d\n", grid->cend);
	printf("astart          = %d\n", grid->astart);
	printf("aend            = %d\n", grid->aend);
	printf("fstart          = %d\n", grid->fstart);
	printf("fend            = %d\n\n", grid->fend);
	
	printf("circle_area     = %f\n", grid->circle_area);
	printf("annulus_area    = %f\n", grid->annulus_area);
}

void grid_print_domain(const struct grid_s *grid, char domain)
{
	int k, newline_modulus;
	int start, end;
	float weight;
	
	if (domain == 'i') {
		printf("Circular Domain\n");
		start = grid->cstart;
		end = grid->cend;
		newline_modulus = grid->cdiami;
	} else {
		printf("Annular Domain\n");
		start = grid->astart;
		end = grid->aend;
		newline_modulus = grid->cdiamo;
	}
	
	for (k = start; k < end; k++) {
		weight = grid->weights[k];

		if (weight == 0) {
			printf(" 0 ");
		} else if (weight == 1) {
			printf(" 1 ");
		} else {
			printf(" A ");
		}

		if ((k - start + 1) % newline_modulus == 0) {
			printf("\n");
		}
	}
	
	printf("\n");
}

static void copy_format_specifier(char *format_buffer, const char *format_specifier)
{
	char *default_format = "%.3e";
	int default_len = strlen(default_format);
	int format_specifier_len;

	memset(format_buffer, '\0', FORMAT_BUFFER_SIZE * sizeof(char));

	if (format_specifier == NULL) {
		strncpy(format_buffer, default_format, default_len);
		format_buffer[default_len] = ' ';
	} else {
		format_specifier_len = strlen(format_specifier);

		if (format_specifier_len < FORMAT_BUFFER_SIZE) {
			strncpy(format_buffer, format_specifier, format_specifier_len);
		} else {
			strncpy(format_buffer, default_format, default_len);
		}
	}
}

void grid_print_weights(const struct grid_s *grid, char domain, const char *format_specifier)
{
	int k, newline_modulus;
	int start, end;
	float weight;
	char format_buffer[FORMAT_BUFFER_SIZE];
	
	copy_format_specifier(format_buffer, format_specifier);
	
	if (domain == 'i') {
		printf("Circular Domain Weights\n");
		start = grid->cstart;
		end = grid->cend;
		newline_modulus = grid->cdiami;
	} else {
		printf("Annular Domain Weights\n");
		start = grid->astart;
		end = grid->aend;
		newline_modulus = grid->cdiamo;
	}
	
	for (k = start; k < end; k++) {
		weight = grid->weights[k];

		printf(format_buffer, weight);

		if ((k - start + 1) % newline_modulus == 0) {
			printf("\n");
		}
	}
	
	printf("\n");
}

void grid_print_state(const struct grid_s *grid, char domain, const char *format_specifier)
{
	int k, newline_modulus;
	int start, end;
	float state;
	char format_buffer[FORMAT_BUFFER_SIZE];
	
	copy_format_specifier(format_buffer, format_specifier);
	
	if (domain == 'i') {
		printf("Circular Domain State\n");
		start = grid->cstart;
		end = grid->cend;
		newline_modulus = grid->cdiami;
	} else {
		printf("Annular Domain State\n");
		start = grid->astart;
		end = grid->aend;
		newline_modulus = grid->cdiamo;
	}
	
	for (k = start; k < end; k++) {
		state = grid->data[k];

		printf(format_buffer, state);

		if ((k - start + 1) % newline_modulus == 0) {
			printf("\n");
		}
	}
	
	printf("\n");
}

void grid_print_indices(const struct grid_s *grid)
{
	int i, j, k;
	int *indicesi, *indicesj;
	
	indicesi = grid->domain_indices;
	indicesj = indicesi + grid->cdiami;
	
	printf("Circular Domain Indices\n");
	for (k = 0; k < grid->cdiami; k++) {
		i = indicesi[k];
		j = indicesj[k];
		
		printf("i = %4d, j = %4d\n", i, j);
	}

	indicesi = grid->domain_indices + 2 * grid->cdiami;
	indicesj = indicesi + grid->cdiamo;
	
	printf("\nAnnular Domain Indices\n");
	for (k = 0; k < grid->cdiamo; k++) {
		i = indicesi[k];
		j = indicesj[k];
		
		printf("i = %4d, j = %4d\n", i, j);
	}
}








