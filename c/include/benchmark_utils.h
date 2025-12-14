#ifndef BENCHMARK_UTILS_H
#define BENCHMARK_UTILS_H

#include <time.h>

#define NANOSEC_TO_SEC 1e-9
#define NSTATS 4

//----------------------------------------------------------------------------
// Enumerated types for sequence types and statistics.
//----------------------------------------------------------------------------

enum sequence_type {
	SEQT_NONE = -1,
	SEQT_LINEAR,
	SEQT_EXPONENTIAL
};

enum statistic_type {
	STAT_NONE = -1,
	STAT_MIN,
	STAT_MAX,
	STAT_MEAN,
	STAT_STD
};

//----------------------------------------------------------------------------
// Struct to hold state of runtime information.
//----------------------------------------------------------------------------

struct runtime_info_s
{
	double *runtimes;
	double *runstats;

	double *mintimes;
	double *maxtimes;
	double *meantimes;
	double *stdtimes;

	int ndims;
	int ntrials;
	int size_min_index;
	int size_max_index;
	int step;
	int base;
	int offset;
	enum sequence_type type;
};

//----------------------------------------------------------------------------
// Forward declarations of functions.
//----------------------------------------------------------------------------

int rinfo_init(struct runtime_info_s *, int, int, int, int, int, int, enum sequence_type);
float rinfo_get_runtime(const struct runtime_info_s *, int, int);
float rinfo_set_runtime(struct runtime_info_s *, int, int, double);
float rinfo_get_stat(const struct runtime_info_s *, int, enum statistic_type);
void rinfo_set_stat(struct runtime_info_s *, int, enum statistic_type);
int rinfo_get_dimension(const struct runtime_info_s *, int);
void rinfo_clear(struct runtime_info_s *);
int rinfo_write_to_file(const struct runtime_info_s *, const char *);

double sample_min(const double *, int);
double sample_max(const double *, int);
double sample_mean(const double *, int);
double sample_var(const double *, int, double);
double sample_std(const double *, int, double);

double compute_elapsed_time(struct timespec *, const struct timespec *, const struct timespec *);
void benchmark(int, int, int, int, int, int, enum sequence_type, double, double, double, int);

#endif