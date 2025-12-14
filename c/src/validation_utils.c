#include "validation_utils.h"
#include <time.h>

#ifdef USE_INTRINSICS
#include "instrinsics_utils.h"
#include <immintrin.h>
#endif


double sample_min(const double *x, int n)
{
	double currmin = x[0];
	int k;

#ifdef USE_INTRINSICS
	__m256d xreg;
	__m256d yreg;
	__m256i loadmask, endmask;
	__m256d ireg = _mm256_set1_ps(DBL_MAX);
	int cutoff = n % FLOAT_PER_M256_REG;
	
	if (cutoff > 0) {
		loadmask = _mm256_set_mask_epi64(cutoff);
		endmask = _mm256_andnot_si256(loadmask, _mm256_set1_epi32(ALL_BITS));
	} else {
		xreg = _mm256_load_ps(x);
	}
	
	for (k = cutoff; k < n; k++) {
		yreg = _mm256_load_ps(x + k);
		xreg = _mm256_min_ps(xreg, yreg);
	}
	
	currmin = _mm256_comp_regmin_ps(xreg);
#else
	for (k = 1; k < n; k++) {
		if (x[k] < currmin) {
			currmin = x[k];
		}
	}
#endif

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

#ifdef USE_INTRINSICS
	int cutoff = n % DOUBLE_PER_M256_REG;
	__m256i mask;
	__m256d sreg = _mm256_set1_pd(0);
	__m256d xreg;
	
	if (cutoff > 0) {
		mask = _mm256_set_mask_epi64(cutoff);
		xreg = _mm256_maskload_pd(x, mask);
		sreg = _mm256_add_pd(sreg, xreg);
	}
	
	for (int k = cutoff; k < n; k += DOUBLE_PER_M256_REG) {
		xreg = _mm256_load_pd(x + k);
		sreg = _mm256_add_pd(sreg, xreg);
	}
	
	sum = _mm256_register_sum_pd(sreg);
#else
	for (int k = 0; k < n; k++) {
		sum += x[k];
	}
#endif

	return sum / n;
}

double sample_var(const double *x, int n, double mean)
{
	double sse = 0;

#ifdef USE_INTRINSICS
	int cutoff = n % DOUBLE_PER_M256_REG;
	__m256i mask;
	__m256d sreg = _mm256_set1_pd(0);
	__m256d mreg = _mm256_set1_pd(mean);
	__m256d dreg;
	__m256d preg;
	__m256d xreg;
	
	if (cutoff > 0) {
		mask = _mm256_set_mask_epi64(cutoff);
		
		xreg = _mm256_maskload_pd(x, mask);
		dreg = _mm256_sub_pd(xreg, mreg);
		preg = _mm256_mul_pd(dreg, dreg);
		sreg = _mm256_add_pd(sreg, preg);
	}
	
	for (int k = cutoff; k < n; k += DOUBLE_PER_M256_REG) {
		xreg = _mm256_load_pd(x + k);
		dreg = _mm256_sub_pd(xreg, mreg);
		preg = _mm256_mul_pd(dreg, dreg);
		sreg = _mm256_add_pd(sreg, preg);
	}
	
	sse = _mm256_register_sum_pd(sreg);
#else
	double xk;

	for (int k = 0; k < n; k++) {
		xk = x[k];

		sse += (xk - mean) * (xk - mean);
	}
#endif
	
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
	
	elapsed_sec = (double)end->tv_sec - (double)start->tv_sec;
	elapsed_nsec = (end->tv_nsec * NANOSEC_TO_SEC) - (start->tv_nsec * NANOSEC_TO_SEC);
	
	return elapsed_sec + elapsed_nsec;
}