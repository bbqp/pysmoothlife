#include <immintrin.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>


#define FLOAT_PER_M256_REG 8
#define INT32_PER_M256_REG FLOAT_PER_M256_REG

#define FLOAT_PER_M128_REG 4
#define INT32_PER_M128_REG FLOAT_PER_M128_REG

#define INT32_ALLBITS 0xFFFFFFFF
#define INT32_HIGHBIT 0x80000000

#define CAPACITY(M, regsize) ((int)ceilf((float)M / regsize))
#define SETBIT(bit) (1 << (bit))
#define SETINTMSB(bit) SETBIT(8)
#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#define MAX(a, b) (((a) > (b)) ? (a) : (b))

void _mm256_modulus(int *a, int b)
{
	// Load a into a register.
	__m256i allmask = _mm_set1_epi32(ALLBITS);
	
	// The divisor in a register.
	__m256i areg = _mm256_maskload_epi32(a, allmask);
	__m256i breg = _mm256_set1_epi32(b);
	
	// Registers for comparing and masking data.
	__m256i zreg = _mm256_set1_epi32(0);

	// Masks for computing elements >= b.
	__m256i gtbmask = _mm256_cmpgt_epi32(areg, breg);
	__m256i eqbmask = _mm256_cmpeq_epi32(areg, breg);
	__m256i gebmask = _mm256_or_si256(gtbmask, eqbmask);
	__m256i ltbmask = _mm256_andnot_si256(gebmask, eqbmask);
	
	// Masks for computing elements < 0.
	__m256i gtzmask = _mm256_cmpgt_epi32(areg, breg);
	__m256i eqzmask = _mm256_cmpeq_epi32(areg, breg);
	__m256i gezmask = _mm256_or_si256(gtzmask, gezmask);
	__m256i ltzmask = _mm256_andnot_si256(gezmask, allmask);

	__m256i outmask = _mm_or_si128(ltzmask, gebmask);
	int some_outside_range = _mm256_testz_si256(outmask, allmask);
	
	if (some_outside_range) {
		// Tackle negative and >= b values in separate registers.
		__m256i anreg = _mm256_and_si256(areg, ltzmask);
		__m256i bnreg = _mm256_and_si256(breg, ltzmask);

		__m256i apreg = _mm256_and_si256(areg, gebmask);
		__m256i bpreg = _mm256_and_si256(breg, gebmask);
		
		while (some_outside_range) {
			anreg = _mm256_add_epi32(anreg, bnreg);
			apreg = _mm256_sub_epi32(apreg, bpreg);

			gtzmask = _mm256_cmpgt_epi32(anreg, breg);
			eqzmask = _mm256_cmpeq_epi32(anreg, breg);
			gezmask = _mm256_or_si256(gtzmask, eqzmask);
			ltzmask = _mm256_andnot_si256(gezmask, allmask);
			
			gtbmask = _mm256_cmpgt_epi32(apreg, breg);
			eqbmask = _mm256_cmpeq_epi32(apreg, breg);
			gebmask = _mm256_or_si256(gtbmask, eqbmask);
			
			some_outside_range = _mm256_testz_si256(ltzmask, gebmask);
		}
		7
		ltzmask = _mm_cmplt_epi32(areg, zreg);
		ltbmask = _mm_cmplt_epi32(areg, breg);
		gebmask = _mm_andnot_si128(ltbmask, allmask);
		
		_mm256_maskstore_epi32(a, ltzmask, anreg);
		_mm256_maskstore_epi32(a, gebmask, apreg);
	}
}


void _mm_modulus(int *a, int b)
{
	// Load a into a register.
	__m128i allmask = _mm_set1_epi32(ALLBITS);
	
	// The divisor in a register.
	__m128i areg = _mm_maskload_epi32(a, allmask);
	__m128i breg = _mm_set1_epi32(b);
	
	// Registers for comparing and masking data.
	__m128i zreg = _mm_set1_epi32(0);

	// Masks for computing elements >= b.
	__m128i ltbmask = _mm_cmplt_epi32(areg, breg);
	__m128i gebmask = _mm_andnot_si128(ltbmask, allmask);
	
	// Masks for computing elements < 0.
	__m128i ltzmask = _mm_cmplt_epi32(areg, zreg);
	__m128i gezmask = _mm_andnot_si128(ltzmask, allmask);

	__m128 outmask = _mm_or_si128(ltzmask, gebmask);
	int some_outside_range = _mm_test_mix_ones_zeros(outmask, allmask);
	
	if (some_outside_range) {
		// Load the array values of a into the register.
		__m128i anreg = _mm_and_si128(areg, ltzmask);
		__m128i bnreg = _mm_and_si128(breg, ltzmask);
		__m128i apreg = _mm_and_si128(areg, gebmask);
		__m128i bpreg = _mm_and_si128(breg, gebmask);
		
		while (some_outside_range) {
			anreg = _mm_add_epi32(anreg, bnreg);
			apreg = _mm_sub_epi32(apreg, bpreg);
			
			ltzmask = _mm_cmplt_epi32(anreg, zreg);
			ltbmask = _mm_cmplt_epi32(apreg, bpreg);
			gebmask = _mm_andnot_si128(ltbmask, allmask);
			
			outmask = _mm_or_si128(ltzmask, gebmask);
			some_outside_range = _mm_test_mix_ones_zeros(outmask, allmask);
		}
		
		ltzmask = _mm_cmplt_epi32(areg, zreg);
		ltbmask = _mm_cmplt_epi32(areg, breg);
		gebmask = _mm_andnot_si128(ltbmask, allmask);
		
		_mm_maskstore_epi32(a, ltzmask, anreg);
		_mm_maskstore_epi32(a, gebmask, apreg);
	}
}

void _mm_mask_modulus(int *a, int b, __m128i loadmask)
{
	__m128i areg = _mm_maskload_epi32(a, loadmask);
	
	// The divisor in a register.
	__m128i breg = _mm_set1_epi32(b);
	
	// Load a into a register.
	__m128i allmask = _mm_set1_epi32(ALLBITS);
	
	// Registers for comparing and masking data.
	__m128i zreg = _mm_set1_epi32(0);

	// Masks for computing elements >= b.
	__m128i ltbmask = _mm_cmplt_epi32(areg, breg);
	__m128i gebmask = _mm_andnot_si128(ltbmask, allmask);
	
	// Masks for computing elements < 0.
	__m128i ltzmask = _mm_cmplt_epi32(areg, zreg);
	__m128i gezmask = _mm_andnot_si128(ltzmask, allmask);

	__m128 outmask = _mm_or_si128(ltzmask, gebmask);
	int some_outside_range = _mm_test_mix_ones_zeros(outmask, allmask);
	
	if (some_outside_range) {
		// Load the array values of a into the register.
		__m128i anreg = _mm_and_si128(areg, ltzmask);
		__m128i apreg = _mm_and_si128(areg, gebmask);
		__m128i bnreg = _mm_and_si128(breg, ltzmask);
		__m128i bpreg = _mm_and_si128(breg, gebmask);
		
		while (some_outside_range) {
			anreg = _mm_add_epi32(anreg, bnreg);
			apreg = _mm_sub_epi32(apreg, bpreg);
			
			ltzmask = _mm_cmplt_epi32(anreg, zreg);
			ltbmask = _mm_cmplt_epi32(apreg, bpreg);
			gebmask = _mm_andnot_si128(ltbmask, allmask);
			
			outmask = _mm_or_si128(ltzmask, gebmask);
			some_outside_range = _mm_test_mix_ones_zeros(outmask, allmask);
		}
		
		ltzmask = _mm_cmplt_epi32(areg, zreg);
		ltbmask = _mm_cmplt_epi32(areg, breg);
		gebmask = _mm_andnot_si128(ltbmask, allmask);
		
		_mm_maskstore_epi32(a, ltzmask, anreg);
		_mm_maskstore_epi32(a, gebmask, apreg);
	}
}


void compute_array_index_modulus(int *indices, int I, int J, int cspan, int M, int N)
{
	int i, j, k = 0;
	int imod, jmod, idxmod;
	
	int ileft, irght;
	int jleft, jrght;
	int jind;
	
	ileft = I - cspan;
	irght = I + cspan;
	jleft = J - cspan;
	jrght = J + cspan;
	
	for (j = jleft; j <= jrght; j++) {
		jind = (j - jleft) * M;
		
		jmod = (j % N) + N;

		for (i = ileft; i <= irght; i++) {
			imod = (i % M) + M;
			
			indices[jind + i - ileft] = J * M + I;
		}
	}
}


void compute_array_index_modulus_intrinsics_m128i(int *indices, int icenter, int cspan, int M)
{
	int i;
	
	int ileft = icenter - cspan;
	int irght = icenter + cspan;
	int iend = irght + CAPACITY(M, BYTES_PER_M128_REGISTER);

	__m128i IA, JA, KA, seqreg, mreg, imask, jmask, allmask;
	int iind[BYTES_PER_M128_REGISTER];
	int jind[BYTES_PER_M128_REGISTER];
	int idx[BYTES_PER_M128_REGISTER];
	
	seqreg = _mm_setr_epi32(0, 1, 2, 3);
	mreg = _mm_set1_epi32(M);
	
	for (j = jleft; j <= jend; j += BYTES_PER_M128_REGISTER) {
			JA = _mm_set1_epi32(j);
			JA = _mm_add_epi32(JA, seqreg);
			JA = _mm_modulus(JA, N);
			JA = _mm_mul_epi32(JA, mreg);

			JA = _mm_set1_epi32(j - jleft - jcutoff);
			JA = _mm_add_epi32(JA, seqreg);
			JA = _mm_mul_epi32(JA, mreg);
			_mm_maskstore_epi32(jind, allmask, JA);
		
		for (i = ileft; i <= iend; i += BYTES_PER_M128_REGISTER) {
			IA = _mm_set1_epi32(i);
			IA = _mm_add_epi32(IA, seqreg);
			IA = _mm_modulus(IA, M);
			KA = _mm_add_epi32(IA, JA);
			_mm_maskstore_epi32(idx, allmask, KA);
			
			IA = _mm_set1_epi32(i - ileft - icutoff);
			IA = _mm_add_epi32(IA, seqreg);
			_mm_maskstore_epi32(iind, allmask, IA);

			indices[jind[0] + iind[0]] = idx[0];
			indices[jind[1] + iind[1]] = idx[1];
			indices[jind[2] + iind[2]] = idx[2];
			indices[jind[3] + iind[3]] = idx[3];
		}
	}
}

void sample_stats(double *runtimes, int n, double *mean, double *std, double *min, double *max)
{
	double sum = 0;
	double sse = 0;
	double currmin, currmax;
	
	for (int k = 0; k < n; k++) {
		sum += runtimes[k];
	}
	
	*mean = sum / n;
	
	for (int k = 0; k < n; k++) {
		sse += (runtimes[k] - *mean) * (runtimes[k] - *mean);
	}
	
	*std = sqrt(sse / (n - 1));
	
	currmin = currmax = runtimes[0];
	
	for (int k = 1; k < n; k++) {
		if (currmin > runtimes[k]) {
			currmin = runtimes[k];
		}
		
		if (currmax < runtimes[k]) {
			currmax = runtimes[k];
		}
	}
	
	*min = currmin;
	*max = currmax;
}

void setrandom(int *array, int M, int N, int lowi, int highi, int lowj, int highj)
{
	int i, j, k;
	int imod = highi - lowi + 1;
	int jmod = highj - lowj + 1;
	int randi, randj;
	int length = M * N;
	
	for (k = 0; k < length; k++) {
		randi = lowi + rand() % imod;
		randj = lowj + rand() % jmod;

		array[k] = randj * M + randi;
	}
}

void setrange(int *array, int M, int N, int I, int J, int cspan)
{
	int i, j;

	for (j = 0; j < N; j++) {
		for (j = 0; j < N; j++) {
			array[j*M + i] = 
		}
	}
}


int array_alloc(int **array, int M, int N, int register_capacity)
{
	int Mcap = CAPACITY(M, register_capacity);
	int Ncap = CAPACITY(N, register_capacity);
	
	*array = calloc(Mcap * Ncap, sizeof(int));
	
	return *array != NULL;
}

void array_free(int **array, int length)
{
	free(*array);
	*array = NULL;
}

int benchmark128(int runs, int Nmin, int Nmax, int Nstep, float dim_ratio)
{
	int m, mcap, ncap;
	int cspan;
	float *runtimes = NULL;
	int *indices = NULL, *array = NULL;
	int I, J;
	
	for (n = Nmin; n <= Nmax; n += Nstep) {
		m = (int)ceilf(n * dim_ratio);
		mcap = CAPACITY(m, BYTES_PER_M128_REGISTER);
		ncap = CAPACITY(n, BYTES_PER_M128_REGISTER);
		
		
		
		cspan = 2 * MIN(m, n) + 1;

		array_alloc(&array, cspan, cspan, BYTES_PER_M128_REGISTER);
		array_alloc(&indices, cspan, cspan, BYTES_PER_M128_REGISTER);
		
		setrandom(array, M, N, -M/2, M/2, -N/2, N/2);
		memcpy(indices, array, Mcap * Ncap * sizeof(int));
		
		compute_array_index_modulus(indices, int I, int J, int cspan, int M, int N)
		
		
		
		
		setrandom(indices, M, N, -M/2, M/2, -N/2, N/2);
		
		
		array_free(&array);
		array_free(&indices);
	}
}

int main(int argc, char *argv)
{
	int M, N;
	int indices;
	srand(0);
	
	return 0;
}







