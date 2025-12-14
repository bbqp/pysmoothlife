#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define NANOSEC_TO_SEC 1e-9

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
		x0 = y0 = -1;
		xn = yn = 1;
		ri_ratio = ((float)(ri_factor + b/2)) * powf(load_scale, size - sizeMin);
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

