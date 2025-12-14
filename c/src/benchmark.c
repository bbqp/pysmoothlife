#include "grid.h"
#include "transitionfunction.h"
#include "benchmark_utils.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <omp.h>

void benchmark(int ntrials, int size_min_index, int size_max_index, int step, int base, int offset, enum sequence_type type, double ri_factor, double bfactor, double load_scale, int seed)
{
	int dimension, size_index, trial;
	struct grid_s grid;
	float res_ratio = 1080 / 1920.0;
	float x0 = -1, xn = 1, y0 = -1, yn = 1;
	int numx, numy;
	float ri_ratio;
	float ri, ro;
	float b;
	float cint, cmean, aint, amean;
	struct timespec start, end, elapsed;
	double elapsed_time;
	struct runtime_info_s rinfo;
	int i, j, jind;
	int alloc_status = rinfo_init(&rinfo, ntrials, size_min_index, size_max_index, step, base, offset, type);
	
	if (!alloc_status) {
		fprintf(stderr, "Unable to allocate structure for run computations, returning.\n");
		
		return;
	}
	
	// Seed the random number generator.
	srand(seed);
	
	for (size_index = rinfo.size_min_index; size_index <= rinfo.size_max_index; size_index++) {
		dimension = rinfo_get_dimension(&rinfo, size_index);
		numx = dimension;
		numy = dimension;
		printf("Setting dimension to %dx%d\n", numx, numy);
		
		ri_ratio = ri_factor * powf(load_scale, size_index - rinfo.size_min_index);
		ri = ri_ratio * (xn - x0);
		ro = 3 * ri;
		b = bfactor * (ro - ri);
		alloc_status = grid_init(&grid, x0, xn, numx, y0, yn, numy, ri_ratio, b);
		
		if (!alloc_status) {
			printf("Grid not allocated, returning.\n");
			rinfo_clear(&rinfo);
			
			return;
		}
		
		// printf("Begin trial runs.\n", dimension);
		for (trial = 0; trial < ntrials; trial++) {
			if ((trial + 1) % 10 == 0) {
				printf("Trial %d\n", trial + 1);
			}

			grid_set_random(&grid);
			
			// Compute the elapsed times for the state function update.
			timespec_get(&start, TIME_UTC);
			grid_update_state_function(&grid);
			timespec_get(&end, TIME_UTC);
			elapsed_time = compute_elapsed_time(&elapsed, &end, &start);

			rinfo_set_runtime(&rinfo, size_index, trial, elapsed_time);
		}
		
		// Runtime statistics for the integral loop.
		rinfo_set_stat(&rinfo, size_index, STAT_MIN);
		rinfo_set_stat(&rinfo, size_index, STAT_MAX);
		rinfo_set_stat(&rinfo, size_index, STAT_MEAN);
		rinfo_set_stat(&rinfo, size_index, STAT_STD);

		// Print out the runtime information.
		printf("\nRuntime Statistics for grid of dimension %3dx%3d with %3d trial runs\n", dimension, dimension, ntrials);
		printf("--------------------------------------------------------------------\n");
		
		printf("\t++Runtime Statistics++\n\n");
		printf("\tCircular Domain Dimension:            %d\n", grid.cspani);
		printf("\tAnnular Domain Dimension:             %d\n\n", grid.cspano);
		printf("\tMinimum Runtime (Seconds):            %e\n", rinfo_get_stat(&rinfo, size_index, STAT_MIN));
		printf("\tMaximum Runtime (Seconds):            %e\n", rinfo_get_stat(&rinfo, size_index, STAT_MAX));
		printf("\tAverage Runtime (Seconds):            %e\n", rinfo_get_stat(&rinfo, size_index, STAT_MEAN));
		printf("\tRuntime Standard Deviation (Seconds): %e\n\n", rinfo_get_stat(&rinfo, size_index, STAT_STD));
		
		// Clear the grids when finished.
		printf("\tClearing the grid objects.\n", trial);
		printf("--------------------------------------------------------------------\n");
		grid_clear(&grid);
	}
	
	rinfo_write_to_file(&rinfo, "benchmark_data.dat");
	rinfo_clear(&rinfo);
}


int main(int argc, char *argv[])
{
	int ntrials = 50;
	int size_min_index = 8, size_max_index = 8;
	int step = 1;
	int base = 2;
	int offset = 1;
	enum sequence_type type = SEQT_EXPONENTIAL;
	double ri_factor = 0.0125;
	double bfactor = 0.25;
	double load_scale = 0.5;
	int seed = 0;
	
	benchmark(ntrials, size_min_index, size_max_index, step, base, offset, type, ri_factor, bfactor, load_scale, seed);

	return 0;
}

