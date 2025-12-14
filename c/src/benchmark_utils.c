#include "benchmark_utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

int rinfo_init(struct runtime_info_s *rinfo, int ntrials, int size_min_index, int size_max_index, int step, int base, int offset, enum sequence_type type)
{
	rinfo->ntrials = ntrials;
	rinfo->ndims = size_max_index - size_min_index + 1;
	rinfo->size_min_index = size_min_index;
	rinfo->size_max_index = size_max_index;
	rinfo->step = step;
	rinfo->base = base;
	rinfo->offset = offset;
	rinfo->type = type;

	rinfo->runtimes = calloc(rinfo->ndims * rinfo->ntrials, sizeof(double));
	rinfo->runstats = calloc(NSTATS * rinfo->ndims, sizeof(double));
	
	if (rinfo->runtimes == NULL || rinfo->runstats == NULL) {
		free(rinfo->runtimes);
		free(rinfo->runstats);

		rinfo->mintimes = rinfo->maxtimes = rinfo->meantimes = rinfo->stdtimes = NULL;
		rinfo->runtimes = rinfo->runstats = NULL;
		rinfo->ntrials = 0;
		rinfo->ndims = 0;
		rinfo->size_min_index = 0;
		rinfo->size_max_index = 0;
		rinfo->step = 0;
		rinfo->base = 0;
		rinfo->offset = 0;
		rinfo->type = SEQT_NONE;
	} else {
		rinfo->mintimes = rinfo->runstats + STAT_MIN * rinfo->ndims;
		rinfo->maxtimes = rinfo->runstats + STAT_MAX * rinfo->ndims;
		rinfo->meantimes = rinfo->runstats + STAT_MEAN * rinfo->ndims;
		rinfo->stdtimes = rinfo->runstats + STAT_STD * rinfo->ndims;
	}
	
	return rinfo->runtimes != NULL && rinfo->runstats != NULL;
}

float rinfo_get_runtime(const struct runtime_info_s *rinfo, int size_index, int trial)
{
	int idx = (size_index - rinfo->size_min_index) * rinfo->ntrials + trial;

	return rinfo->runtimes[idx];
}

float rinfo_set_runtime(struct runtime_info_s *rinfo, int size_index, int trial, double runtime)
{
	int idx = (size_index - rinfo->size_min_index) * rinfo->ntrials + trial;

	return rinfo->runtimes[idx] = runtime;
}

float rinfo_get_stat(const struct runtime_info_s *rinfo, int size_index, enum statistic_type stype)
{
	int idx = size_index - rinfo->size_min_index;
	double runtime_stat = 0;
	
	switch (stype)
	{
		case STAT_MIN:
			runtime_stat = rinfo->mintimes[idx];
			break;
		
		case STAT_MAX:
			runtime_stat = rinfo->maxtimes[idx];
			break;
		
		case STAT_MEAN:
			runtime_stat = rinfo->meantimes[idx];
			break;
		
		case STAT_STD:
			runtime_stat = rinfo->stdtimes[idx];
	}
	
	return runtime_stat;
}

void rinfo_set_stat(struct runtime_info_s *rinfo, int size_index, enum statistic_type stype)
{
	int idx = size_index - rinfo->size_min_index;
	int trial_idx = idx * rinfo->ntrials;
	double mean;
	
	switch (stype)
	{
		case STAT_MIN:
			rinfo->mintimes[idx] = sample_min(rinfo->runtimes + trial_idx, rinfo->ntrials);
			break;
		
		case STAT_MAX:
			rinfo->maxtimes[idx] = sample_max(rinfo->runtimes + trial_idx, rinfo->ntrials);
			break;
		
		case STAT_MEAN:
			rinfo->meantimes[idx] = sample_mean(rinfo->runtimes + trial_idx, rinfo->ntrials);
			break;
		
		case STAT_STD:
			mean = sample_mean(rinfo->runtimes + trial_idx, rinfo->ntrials);
			rinfo->stdtimes[idx] = sample_std(rinfo->runtimes + trial_idx, rinfo->ntrials, mean);
	}
}

int rinfo_get_dimension(const struct runtime_info_s *rinfo, int size_index)
{
	int dimension = 1, factor = 1;
	
	if (rinfo->type == SEQT_LINEAR) {
		dimension = size_index * rinfo->step + rinfo->offset;
	} else if(rinfo->type == SEQT_EXPONENTIAL) {
		for (int i = 0; i < rinfo->step; i++) {
			factor *= rinfo->base;
		}

		for (int i = 0; i < size_index; i++) {
			dimension *= factor;
		}
		
		dimension += rinfo->offset;
	}
	
	return dimension;
}

void rinfo_clear(struct runtime_info_s *rinfo)
{
	free(rinfo->runtimes);
	free(rinfo->runstats);

	rinfo->mintimes = rinfo->maxtimes = rinfo->meantimes = rinfo->stdtimes = NULL;
	rinfo->runtimes = rinfo->runstats = NULL;
	rinfo->ntrials = 0;
	rinfo->ndims = 0;
	rinfo->size_min_index = 0;
	rinfo->size_max_index = 0;
	rinfo->step = 0;
	rinfo->offset = 0;
	rinfo->base = 0;
	rinfo->type = SEQT_NONE;
}


int rinfo_write_to_file(const struct runtime_info_s *rinfo, const char *filename)
{
	FILE *outfile = NULL;
	char *type = (rinfo->type == SEQT_LINEAR) ? "linear" : "exponential";
	int typelen = strlen(type);
	int nstats = NSTATS;
	int open_status = 0;

	outfile = fopen(filename, "wb");
	open_status = (outfile != NULL);

	if (open_status) {
		fwrite(&rinfo->ntrials, sizeof(int), 1, outfile);
		fwrite(&rinfo->ndims, sizeof(int), 1, outfile);
		fwrite(&rinfo->size_min_index, sizeof(int), 1, outfile);
		fwrite(&rinfo->size_max_index, sizeof(int), 1, outfile);
		fwrite(&rinfo->step, sizeof(int), 1, outfile);
		fwrite(&rinfo->offset, sizeof(int), 1, outfile);
		fwrite(&rinfo->base, sizeof(int), 1, outfile);
		fwrite(&nstats, sizeof(int), 1, outfile);
		fwrite(&typelen, sizeof(int), 1, outfile);
		fwrite(type, sizeof(char), strlen(type), outfile);
		fwrite(rinfo->runtimes, sizeof(double), rinfo->ndims * rinfo->ntrials, outfile);
		fwrite(rinfo->runstats, sizeof(double), nstats * rinfo->ndims, outfile);

		fclose(outfile);
	} else {
		printf("ERROR: Unable to write output data to file.\n");
	}
	
	return open_status;
}

//----------------------------------------------------------------------------
// Helper functions for computing statistics.
//----------------------------------------------------------------------------

double sample_min(const double *x, int n)
{
	double currmin = x[0];
	int k;

	for (k = 1; k < n; k++) {
		if (x[k] < currmin) {
			currmin = x[k];
		}
	}

	return currmin;
}

double sample_max(const double *x, int n)
{
	int maxidx = 0;
	double currmax = x[0];
	
	for (int k = 1; k < n; k++) {
		if (x[k] > currmax) {
			maxidx = k;
			currmax = x[k];
		}
	}
	
	return currmax;
}

double sample_mean(const double *x, int n)
{
	double sum = 0;
	
	for (int k = 0; k < n; k++) {
		sum += x[k];
	}
	
	return sum / n;
}

double sample_var(const double *x, int n, double mean)
{
	double sse = 0;
	double xk;
	
	for (int k = 0; k < n; k++) {
		xk = x[k];

		sse += (xk - mean) * (xk - mean);
	}
	
	return sse / n;
}

double sample_std(const double *x, int n, double mean)
{
	return sqrt(sample_var(x, n, mean));
}

double compute_elapsed_time(struct timespec *elapsed, const struct timespec *end, const struct timespec *start)
{
	double elapsed_sec;
	double elapsed_nsec;
	// elapsed->tv_sec = end->tv_sec - start->tv_sec;

	// if (start->tv_nsec > end->tv_nsec) {
		// elapsed->tv_sec -= 1;
		// elapsed->tv_nsec = start->tv_nsec - end->tv_nsec;
	// } else {
		// elapsed->tv_nsec = end->tv_nsec - start->tv_nsec;
	// }
	
	elapsed_sec = end->tv_sec - start->tv_sec;
	elapsed_nsec = (end->tv_nsec * NANOSEC_TO_SEC) - (start->tv_nsec * NANOSEC_TO_SEC);
	
	return elapsed_sec + elapsed_nsec;
}
