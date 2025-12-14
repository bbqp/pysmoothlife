#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define NANOSEC_TO_SEC 1e-9

enum grid_data_type {
	GDT_NONE = -1,
	GDT_FLOAT,
	GDT_DOUBLE
};

typedef union ftype
{
	float f32;
	double f64;
} fltype;

struct grid_s
{
	void *data;
	void *cdist;
	int *idx;
	int cistart, ciend, aistart, aiend;
	int cjstart, cjend, ajstart, ajend;
	int cstart, cend;
	int astart, aend;
	int fstart, fend;
	fltype x0, xn;
	fltype y0, yn;
	fltype dx, dy;
	fltype ri, ri_ratio, ro;
	int cspani, cspano;
	int numx;
	int numy;
	enum grid_data_type type;
};


int grid_init(struct grid_s *grid, fltype x0, fltype xn, int numx, fltype y0, fltype yn, int numy, fltype ri_ratio, enum grid_data_type type)
{
	fltype deltax, deltay, mindelta, mindvar;
	
	grid->type = type;

	if (type == GDT_FLOAT) {
		grid->x0.f32 = x0.f32;
		grid->xn.f32 = xn.f32;
		
		grid->y0.f32 = y0.f32;
		grid->yn.f32 = yn.f32;
		
		deltax.f32 = grid->xn.f32 - grid->x0.f32;
		deltay.f32 = grid->yn.f32 - grid->y0.f32;

		grid->numx = numx;
		grid->numy = numy;
		
		grid->dx.f32 = deltax.f32 / (grid->numx - 1);
		grid->dy.f32 = deltay.f32 / (grid->numy - 1);

		mindelta.f32 = deltax.f32 < deltay.f32 ? deltax.f32 : deltay.f32;
		mindvar.f32 = grid->dx.f32 < grid->dy.f32 ? grid->dx.f32 : grid->dy.f32;
		
		grid->ri.f32 = ri_ratio.f32 * mindelta.f32;
		grid->ri_ratio.f32 = ri_ratio.f32;
		grid->ro.f32 = 3 * grid->ri.f32;
		
		grid->cspani = 2 * (int)ceil(grid->ri.f32 / mindvar.f32) + 2;
		grid->cspano = 2 * (int)ceil(grid->ro.f32 / mindvar.f32) + 2;
		
		grid->cstart = 0;
		grid->cend = grid->cspani * grid->cspani - 1;
		grid->astart = grid->cend + 1;
		grid->aend = grid->astart + grid->cspano * grid->cspano;
		grid->fstart = grid->aend + 1;
		grid->fend = grid->fstart + (grid->numx - 1) * (grid->numy - 1);
		
		
		grid->data = calloc(grid->fend + 1, sizeof(float));
		
		// printf("x0 = %f\n", grid->x0.f32);
		// printf("xn = %f\n", grid->xn.f32);
		// printf("y0 = %f\n", grid->y0.f32);
		// printf("yn = %f\n", grid->yn.f32);

		// printf("deltax = %f\n", deltax.f32);
		// printf("deltay = %f\n", deltay.f32);
		// printf("mindvar = %f\n", mindvar.f32);

		// printf("ri = %f\n", grid->ri.f32);
		// printf("ri_ratio = %f\n", grid->ri_ratio.f32);
		// printf("ro = %f\n", grid->ro.f32);
		
		// printf("cspani = %d\n", grid->cspani);
		// printf("cspano = %d\n", grid->cspano);
		
		// printf("numx = %d\n", grid->numx);
		// printf("numy = %d\n", grid->numy);
		
		// printf("deltax = %f\n", grid->dx.f32);
		// printf("deltay = %f\n", grid->dy.f32);
	} else {
		grid->x0.f64 = x0.f64;
		grid->xn.f64 = xn.f64;
		
		grid->y0.f64 = y0.f64;
		grid->yn.f64 = yn.f64;
		
		deltax.f64 = grid->xn.f64 - grid->x0.f64;
		deltay.f64 = grid->yn.f64 - grid->y0.f64;
		
		grid->numx = numx;
		grid->numy = numy;
		
		grid->dx.f64 = deltax.f64 / (grid->numx - 1);
		grid->dy.f64 = deltay.f64 / (grid->numy - 1);
		
		mindelta.f64 = deltax.f64 < deltay.f64 ? deltax.f64 : deltay.f64;
		mindvar.f64 = grid->dx.f64 < grid->dy.f64 ? grid->dx.f64 : grid->dy.f64;
		
		grid->ri.f64 = ri_ratio.f64 * mindelta.f64;
		grid->ri_ratio.f64 = ri_ratio.f64;
		grid->ro.f64 = 3 * grid->ri.f64;
		
		grid->cspani = 2 * (int)ceil(grid->ri.f64 / mindvar.f64) + 2;
		grid->cspano = 2 * (int)ceil(grid->ro.f64 / mindvar.f64) + 2;
		
		grid->cstart = 0;
		grid->cend = grid->cspani * grid->cspani - 1;
		grid->astart = grid->cend + 1;
		grid->aend = grid->astart + grid->cspano * grid->cspano;
		grid->fstart = grid->aend + 1;
		grid->fend = grid->fstart + (grid->numx - 1) * (grid->numy - 1);
		
		grid->data = calloc(grid->fend + 1, sizeof(double));
		
		
		grid->cistart = 0;
		grid->ciend = grid->cspani - 1;
		grid->cjstart = grid->ciend + 1;
		grid->cjend = grid->cjstart + grid->cspani - 1;
		
		grid->aistart = grid->cjend;
		grid->aiend = grid->aistart + grid->cspano - 1;
		grid->ajstart = grid->aiend + 1;
		grid->ajend + grid->ajstart + grid->cspano - 1;
		
		grid->idx = calloc(grid->ajend + 1, sizeof(int));
		
		// printf("x0 = %lf\n", grid->x0.f64);
		// printf("xn = %lf\n", grid->xn.f64);
		// printf("y0 = %lf\n", grid->y0.f64);
		// printf("yn = %lf\n", grid->yn.f64);

		// printf("deltax = %lf\n", deltax.f64);
		// printf("deltay = %lf\n", deltay.f64);
		// printf("mindvar = %lf\n", mindvar.f64);

		// printf("ri = %lf\n", grid->ri.f64);
		// printf("ri_ratio = %lf\n", grid->ri_ratio.f64);
		// printf("ro = %lf\n", grid->ro.f64);
		
		// printf("cspani = %d\n", grid->cspani);
		// printf("cspano = %d\n", grid->cspano);
		
		// printf("numx = %d\n", grid->numx);
		// printf("numy = %d\n", grid->numy);
		
		// printf("deltax = %lf\n", grid->dx.f64);
		// printf("deltay = %lf\n", grid->dy.f64);
	}
	
	if(grid->data == NULL || grid->idx) {
		free(grid->data);
		free(grid->idx);
		
		grid->data = NULL;
		grid->idx = NULL;
	}
	
	return grid->data != NULL && grid->idx != NULL;
}

void grid_clear(struct grid_s *grid)
{
	free(grid->data);
	grid->data = NULL;
	grid->numx = grid->numy = 0;
	grid->type = GDT_NONE;
}

void grid_set_index_span(struct grid_s *grid, int i, int j)
{
	int NC = grid->cspani;
	int *Cind = grid->idx;
	int *Aind = grid->idx + NC * NC;
	int M = grid->numx - 1;
	int N = grid->numy - 1;
	
	int ileft = i - (grid->cspani - 2) / 2;
	int irght = i + (grid->cspani - 2) / 2;
	int jleft = j - (grid->cspani - 2) / 2;
	int jrght = j + (grid->cspani - 2) / 2;
	int k = 0;
	
	for (int j = jleft; j <= jrght; j++) {
		jind = ((j + N) % N) * M;

		for (int i = ileft; i <= irght; i++) {
			Cind[k++] = jind + ((i + M) % M);
		}
	}
	
	k = 0;
	ileft = i - (grid->cspano - 2) / 2;
	irght = i + (grid->cspano - 2) / 2;
	jleft = j - (grid->cspano - 2) / 2;
	jrght = j + (grid->cspano - 2) / 2;

	for (int j = jleft; j <= jrght; j++) {
		jind = ((j + N) % N) * M;

		for (int i = ileft; i <= irght; i++) {
			Aind[k++] = jind + ((i + M) % M);
		}
	}
	
}

void set_random(struct grid_s *grid)
{
	int k;
	int NC = grid->cspani;
	int NA = grid->cspano;
	int start = 0, end = NC * NC + NA * NA;
	float *dataf;
	double *datad;
	
	//printf("In set_random(): NC = %d, NA = %d\n", NC, NA);
	
	if (grid->type == GDT_FLOAT) {
		//printf("grid->type = %d\n", grid->type);
		dataf = (float *)(grid->data);

		for (k = start; k < end; k++) {
			//printf("\tSetting grid->data[%]\n", k);
			dataf[k] = 2 * ((float)rand() / RAND_MAX) - 1;
		}
	} else {
		//printf("grid->type = %d\n", grid->type);
		datad = (double *)(grid->data);

		for (k = start; k < end; k++) {
			datad[k] = 2 * ((double)rand() / RAND_MAX) - 1;
		}
	}
}

fltype compute_sum(struct grid_s *grid, char which)
{
	int k;
	int NC = grid->cspani;
	int NA = grid->cspano;
	int start, end;
	fltype sum;
	float sumf = 0, *dataf;
	double sumd = 0, *datad;
	
	if (which == 'i') {
		start = 0;
		end = NC * NC;
	} else {
		start = NC * NC;
		end = start + NA * NA;
	}

	if (grid->type == GDT_FLOAT) {
		// printf("Function compute_sum(): In GDT_FLOAT clause\n");
		dataf = (float *)(grid->data);
		
		// printf("\tstart = %d\n", start);
		// printf("\tend = %d\n", end);
		// printf("\tdataf = %X\n", dataf);

		for (k = start; k < end; k++) {
			sumf += dataf[k];
		}

		sum.f32 = sumf;
	} else {
		// printf("Function compute_sum(): In GDT_FLOAT clause\n");
		datad = (double *)(grid->data);

		// printf("\tstart = %d\n", start);
		// printf("\tend = %d\n", end);
		// printf("\tdatad = %X\n", datad);

		for (k = start; k < end; k++) {
			sumd += datad[k];
		}

		sum.f64 = sumd;
	}
		
	return sum;
}

fltype compute_sum_at(struct grid_s *grid, int i, int j, char which)
{
	int i, j, jind;
	fltype sum;
	float sumf = 0, *dataf;
	double sumd = 0, *datad;
	int *CI = grid->idx;
	int *CJ = grid->idx + grid->cjstart;
	int *AI = grid->idx + grid->aistart;
	int *AJ = grid->idx + grid->ajstart;
	int M = grid->numx - 1;
	int N = grid->numy - 1;
	
	grid_set_index_span(grid, i, j);

	if (grid->type == GDT_FLOAT) {
		// printf("Function compute_sum(): In GDT_FLOAT clause\n");
		dataf = (float *)(grid->data);
		
		// printf("\tstart = %d\n", start);
		// printf("\tend = %d\n", end);
		// printf("\tdataf = %X\n", dataf);

		if (which == 'i') {
			for (j = 0; j < grid->cspani; j++) {
				jind = CJ[j] * M;
				
				for (i = 0; i < grid->cspani; i++) {
					sumf += dataf[jind + CI[i]];
				}
			}
		} else {
			for (j = 0; j < grid->cspano; j++) {
				jind = AJ[j] * M;
				
				for (i = 0; i < grid->cspano; i++) {
					sumf += dataf[jind + AI[i]];
				}
			}
		}

		sum.f32 = sumf;
	} else {
		// printf("Function compute_sum(): In GDT_FLOAT clause\n");
		datad = (double *)(grid->data);

		// printf("\tstart = %d\n", start);
		// printf("\tend = %d\n", end);
		// printf("\tdatad = %X\n", datad);

		if (which == 'i') {
			for (j = 0; j < grid->cspani; j++) {
				jind = CJ[j] * M;
				
				for (i = 0; i < grid->cspani; i++) {
					sumd += datad[jind + CI[i]];
				}
			}
		} else {
			for (j = 0; j < grid->cspano; j++) {
				jind = AJ[j] * M;
				
				for (i = 0; i < grid->cspano; i++) {
					sumd += datad[jind + AI[i]];
				}
			}
		}

		sum.f64 = sumd;
	}
		
	return sum;
}

double samplemin(const double *x, int n)
{
	double currmin = x[0];
	
	for (int k = 1; k < n; k++) {
		if (x[k] < currmin) {
			currmin = x[k];
		}
	}
	
	return currmin;
}

double samplemax(const double *x, int n)
{
	double currmax = x[0];
	
	for (int k = 1; k < n; k++) {
		if (x[k] > currmax) {
			currmax = x[k];
		}
	}
	
	return currmax;
}

double sampleave(const double *x, int n)
{
	double sum = 0;
	
	for (int k = 0; k < n; k++) {
		sum += x[k];
	}
	
	return sum / n;
}

double samplevar(const double *x, int n, double mean)
{
	double sse = 0;
	double xk;
	
	for (int k = 0; k < n; k++) {
		xk = x[k];

		sse += (xk - mean) * (xk - mean);
	}
	
	return sse / n;
}

double samplestd(const double *x, int n, double mean)
{
	return sqrt(samplevar(x, n, mean));
}

double compute_elapsed_time(struct timespec *elapsed, const struct timespec *end, const struct timespec *start)
{
	elapsed->tv_sec = end->tv_sec - start->tv_sec;

	if (start->tv_nsec > end->tv_nsec) {
		elapsed->tv_sec -= 1;
		elapsed->tv_nsec = start->tv_nsec - end->tv_nsec;
	} else {
		elapsed->tv_nsec = end->tv_nsec - start->tv_nsec;
	}
	
	return (double)elapsed->tv_sec + ((double)elapsed->tv_nsec) * NANOSEC_TO_SEC;
}

void benchmark(int sizeMin, int sizeMax, int ntrials, double ri_factor, double b, double load_scale, int seed)
{
	int ndims = sizeMax - sizeMin + 1;
	int N = ndims * ntrials;
	int dimension, size, idx, trial;
	double *runtimes = calloc(2 * N, sizeof(double));
	double *runstats = calloc(4 * (2 * ndims), sizeof(double));
	double *mintimes;
	double *maxtimes;
	double *avetimes;
	double *stdtimes;
	double *ftimes;
	double *dtimes;
	struct grid_s fgrid, dgrid;
	fltype x0, xn, y0, yn;
	fltype ri_ratio;
	struct timespec start, end, elapsed;
	int alloc_status;
	
	if (runtimes == NULL || runstats == NULL) {
		free(runtimes);
		free(runstats);
		
		fprintf(stderr, "Unable to allocate %dx%d array for run computations.\n", ntrials, sizeMax - sizeMin + 1);
		
		return;
	}
	
	mintimes = runstats;
	maxtimes = runstats + 2 * ndims;
	avetimes = runstats + 4 * ndims;
	stdtimes = runstats + 6 * ndims;
	ftimes = runtimes;
	dtimes = runtimes + N;
	
	// Seed the random number generator.
	srand(seed);
	
	for (size = sizeMin; size <= sizeMax; size++) {
		// Index for the runtime information.
		idx = (size - sizeMin) * ntrials;
		
		dimension = (1 << size) + 1;
		printf("Setting dimension to %d\n", dimension);
		
		//printf("Initializing single-precision grid.\n", dimension);
		x0.f32 = y0.f32 = -1;
		xn.f32 = yn.f32 = 1;
		ri_ratio.f32 = ((float)(ri_factor + b/2)) * powf(load_scale, size - sizeMin);
		alloc_status = grid_init(&fgrid, x0, xn, dimension, y0, yn, dimension, ri_ratio, GDT_FLOAT);
		
		if (!alloc_status) {
			printf("Single-precision grid not allocated, returning.\n");
			
			return;
		}

		//printf("Initializing double-precision grid.\n", dimension);
		x0.f64 = y0.f64 = -1;
		xn.f64 = yn.f64 = 1;
		ri_ratio.f64 = (ri_factor + b/2) * pow(load_scale, size - sizeMin);
		alloc_status = grid_init(&dgrid, x0, xn, dimension, y0, yn, dimension, ri_ratio, GDT_DOUBLE);
		
		if (!alloc_status) {
			printf("Double-precision grid not allocated, returning.\n");
			grid_clear(&fgrid);
			
			return;
		}
		
		// printf("Begin trial runs.\n", dimension);
		for (trial = 0; trial < ntrials; trial++) {
			//printf("Trial %d\n", trial);
			// Compute the elapsed times for the float sum.
			//timespec_get(&start, TIME_UTC);
			// printf("\tSetting random values on single-precision grid.\n");
			set_random(&fgrid);
			// printf("\tComputing inner sum on single-precision grid.\n");
			timespec_get(&start, TIME_UTC);
			compute_sum(&fgrid, 'i');
			// printf("\tComputing outer sum on single-precision grid.\n");
			compute_sum(&fgrid, 'o');
			timespec_get(&end, TIME_UTC);
			ftimes[idx + trial] = compute_elapsed_time(&elapsed, &end, &start);

			// Compute the elapsed times for the double sum.
			//timespec_get(&start, TIME_UTC);
			// printf("\tSetting random values on double-precision grid.\n");
			set_random(&dgrid);
			// printf("\tComputing inner sum on double-precision grid.\n");
			timespec_get(&start, TIME_UTC);
			compute_sum(&dgrid, 'i');
			// printf("\tComputing outer sum on double-precision grid.\n");
			compute_sum(&dgrid, 'o');
			timespec_get(&end, TIME_UTC);
			dtimes[idx + trial] = compute_elapsed_time(&elapsed, &end, &start);
		}
		
		// Runtime statistics for single precision loop.
		mintimes[size - sizeMin] = samplemin(ftimes + idx, ntrials);
		maxtimes[size - sizeMin] = samplemax(ftimes + idx, ntrials);
		avetimes[size - sizeMin] = sampleave(ftimes + idx, ntrials);
		stdtimes[size - sizeMin] = samplestd(ftimes + idx, ntrials, avetimes[size - sizeMin]);
		
		// Runtime statistics for double precision loop.
		mintimes[ndims + size - sizeMin] = samplemin(dtimes + idx, ntrials);
		maxtimes[ndims + size - sizeMin] = samplemax(dtimes + idx, ntrials);
		avetimes[ndims + size - sizeMin] = sampleave(dtimes + idx, ntrials);
		stdtimes[ndims + size - sizeMin] = samplestd(dtimes + idx, ntrials, avetimes[ndims + size - sizeMin]);

		// Print out the runtime information.
		printf("\nRuntime Statistics for grid of dimension %3dx%3d with %3d trial runs\n", dimension, dimension, ntrials);
		printf("--------------------------------------------------------------------\n");
		
		printf("\t++Single Precision Runtime Statistics++\n\n");
		printf("\tCircular Domain Dimension:            %d\n", fgrid.cspani);
		printf("\tAnnular Domain Dimension:             %d\n\n", fgrid.cspano);
		printf("\tMinimum Runtime (Seconds):            %e\n", mintimes[size - sizeMin]);
		printf("\tMaximum Runtime (Seconds):            %e\n", maxtimes[size - sizeMin]);
		printf("\tAverage Runtime (Seconds):            %e\n", avetimes[size - sizeMin]);
		printf("\tRuntime Standard Deviation (Seconds): %e\n\n", stdtimes[size - sizeMin]);
		
		printf("\t++Double Precision Runtime Statistics++\n\n");
		printf("\tCircular Domain Dimension:            %d\n", dgrid.cspani);
		printf("\tAnnular Domain Dimension:             %d\n\n", dgrid.cspano);
		printf("\tMinimum Runtime (Seconds):            %e\n", mintimes[ndims + size - sizeMin]);
		printf("\tMaximum Runtime (Seconds):            %e\n", maxtimes[ndims + size - sizeMin]);
		printf("\tAverage Runtime (Seconds):            %e\n", avetimes[ndims + size - sizeMin]);
		printf("\tRuntime Standard Deviation (Seconds): %e\n\n", stdtimes[ndims + size - sizeMin]);

		
		// Clear the grids when finished.
		printf("\tClearing the grid objects.\n", trial);
		printf("--------------------------------------------------------------------\n");
		grid_clear(&fgrid);
		grid_clear(&dgrid);
	}
	
	free(runtimes);
	free(runstats);
}


int main(int argc, char *argv[])
{
	int sizeMin = 5, sizeMax = 11;
	int ntrials = 100;
	double load_scale = 1.0;
	double ri_factor = 0.05;
	double b = 0;
	int seed = 0;
	
	benchmark(sizeMin, sizeMax, ntrials, ri_factor, b, load_scale, seed);

	return 0;
}












