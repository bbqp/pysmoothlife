#include "intrinsics_utils.h"
#include <immintrin.h>
#include <stdio.h>
	
#define M128_LPERM_MASK (0x39 + ((0x4e) << 8) + ((0x93) << 16) + ((0xe4) << 24))
#define M128_LPERM_TO_IMM8(NPERMS) ((M128_LPERM_MASK >> ((NPERMS) - 1)) & (INT32_ALLBITS))

#define LPERM1 0x39
#define LPERM2 0x4e
#define LPERM3 0x93

#define M128_RPERM_MASK (0x93 + ((0x4e) << 8) + ((0x39) << 16) + ((0xe4) << 24))
#define M128_RPERM_TO_IMM8(NPERMS) ((M128_RPERM_MASK >> ((NPERMS) - 1)) & (INT32_ALLBITS))


//----------------------------------------------------------------------------
// Functions for creating masks.
//----------------------------------------------------------------------------

__m128i _mm_setmask_fromto_epi32(int from, int to)
{
	__m128i mask;
	
	if (from > to || from >= INT32_PER_M128_REG) {
		mask = _mm_set1_epi32(0);
	} else {
		switch (from) {
			case 0:
				switch(to) {
					case 0:
						mask = _mm_setr_epi32(INT32_ALLBITS, 0, 0, 0);
						break;
					case 1:
						mask = _mm_setr_epi32(INT32_ALLBITS, INT32_ALLBITS, 0, 0);
						break;
					case 2:
						mask = _mm_setr_epi32(INT32_ALLBITS, INT32_ALLBITS, INT32_ALLBITS, 0);
						break;
					default:
						mask = _mm_set1_epi32(INT32_ALLBITS);
				}
				break;
			case 1:
				switch(to) {
					case 1:
						mask = _mm_setr_epi32(0, INT32_ALLBITS, 0, 0);
						break;
					case 2:
						mask = _mm_setr_epi32(0, INT32_ALLBITS, INT32_ALLBITS, 0);
						break;
					default:
						mask = _mm_setr_epi32(0, INT32_ALLBITS, INT32_ALLBITS, INT32_ALLBITS);
				}
				break;
			case 2:
				switch(to) {
					case 2:
						mask = _mm_setr_epi32(0, 0, INT32_ALLBITS, 0);
						break;
					default:
						mask = _mm_setr_epi32(0, 0, INT32_ALLBITS, INT32_ALLBITS);
				}
				break;
			case 3:
				mask = _mm_setr_epi32(0, 0, 0, INT32_ALLBITS);
		}
	}
	
	return mask;
}

__m256i _mm256_setmask_fromto_epi32(int from, int to)
{
	__m256i mask;
	
	if (from > to || from >= INT32_PER_M256_REG) {
		mask = _mm256_set1_epi32(0);
	} else {
		switch (from) {
			case 0:
				switch(to) {
					case 0:
						mask = _mm256_setr_epi32(INT32_ALLBITS, 0, 0, 0, 0, 0, 0, 0);
						break;
					case 1:
						mask = _mm256_setr_epi32(INT32_ALLBITS, INT32_ALLBITS, 0, 0, 0, 0, 0, 0);
						break;
					case 2:
						mask = _mm256_setr_epi32(INT32_ALLBITS, INT32_ALLBITS, INT32_ALLBITS, 0, 0, 0, 0, 0);
						break;
					case 3:
						mask = _mm256_setr_epi32(INT32_ALLBITS, INT32_ALLBITS, INT32_ALLBITS, INT32_ALLBITS, 0, 0, 0, 0);
						break;
					case 4:
						mask = _mm256_setr_epi32(INT32_ALLBITS, INT32_ALLBITS, INT32_ALLBITS, INT32_ALLBITS, INT32_ALLBITS, 0, 0, 0);
						break;
					case 5:
						mask = _mm256_setr_epi32(INT32_ALLBITS, INT32_ALLBITS, INT32_ALLBITS, INT32_ALLBITS, INT32_ALLBITS, INT32_ALLBITS, 0, 0);
						break;
					case 6:
						mask = _mm256_setr_epi32(INT32_ALLBITS, INT32_ALLBITS, INT32_ALLBITS, INT32_ALLBITS, INT32_ALLBITS, INT32_ALLBITS, INT32_ALLBITS, 0);
						break;
					default:
						mask = _mm256_set1_epi32(INT32_ALLBITS);
				}
				break;
			case 1:
				switch(to) {
					case 1:
						mask = _mm256_setr_epi32(0, INT32_ALLBITS, 0, 0, 0, 0, 0, 0);
						break;
					case 2:
						mask = _mm256_setr_epi32(0, INT32_ALLBITS, INT32_ALLBITS, 0, 0, 0, 0, 0);
						break;
					case 3:
						mask = _mm256_setr_epi32(0, INT32_ALLBITS, INT32_ALLBITS, INT32_ALLBITS, 0, 0, 0, 0);
						break;
					case 4:
						mask = _mm256_setr_epi32(0, INT32_ALLBITS, INT32_ALLBITS, INT32_ALLBITS, INT32_ALLBITS, 0, 0, 0);
						break;
					case 5:
						mask = _mm256_setr_epi32(0, INT32_ALLBITS, INT32_ALLBITS, INT32_ALLBITS, INT32_ALLBITS, INT32_ALLBITS, 0, 0);
						break;
					case 6:
						mask = _mm256_setr_epi32(0, INT32_ALLBITS, INT32_ALLBITS, INT32_ALLBITS, INT32_ALLBITS, INT32_ALLBITS, INT32_ALLBITS, 0);
						break;
					default:
						mask = _mm256_setr_epi32(0, INT32_ALLBITS, INT32_ALLBITS, INT32_ALLBITS, INT32_ALLBITS, INT32_ALLBITS, INT32_ALLBITS, INT32_ALLBITS);
				}
				break;
			case 2:
				switch(to) {
					case 2:
						mask = _mm256_setr_epi32(0, 0, INT32_ALLBITS, 0, 0, 0, 0, 0);
						break;
					case 3:
						mask = _mm256_setr_epi32(0, 0, INT32_ALLBITS, INT32_ALLBITS, 0, 0, 0, 0);
						break;
					case 4:
						mask = _mm256_setr_epi32(0, 0, INT32_ALLBITS, INT32_ALLBITS, INT32_ALLBITS, 0, 0, 0);
						break;
					case 5:
						mask = _mm256_setr_epi32(0, 0, INT32_ALLBITS, INT32_ALLBITS, INT32_ALLBITS, INT32_ALLBITS, 0, 0);
						break;
					case 6:
						mask = _mm256_setr_epi32(0, 0, INT32_ALLBITS, INT32_ALLBITS, INT32_ALLBITS, INT32_ALLBITS, INT32_ALLBITS, 0);
						break;
					default:
						mask = _mm256_setr_epi32(0, 0, INT32_ALLBITS, INT32_ALLBITS, INT32_ALLBITS, INT32_ALLBITS, INT32_ALLBITS, INT32_ALLBITS);
				}
				break;
			case 3:
				switch(to) {
					case 3:
						mask = _mm256_setr_epi32(0, 0, 0, INT32_ALLBITS, 0, 0, 0, 0);
						break;
					case 4:
						mask = _mm256_setr_epi32(0, 0, 0, INT32_ALLBITS, INT32_ALLBITS, 0, 0, 0);
						break;
					case 5:
						mask = _mm256_setr_epi32(0, 0, 0, INT32_ALLBITS, INT32_ALLBITS, INT32_ALLBITS, 0, 0);
						break;
					case 6:
						mask = _mm256_setr_epi32(0, 0, 0, INT32_ALLBITS, INT32_ALLBITS, INT32_ALLBITS, INT32_ALLBITS, 0);
						break;
					default:
						mask = _mm256_setr_epi32(0, 0, 0, INT32_ALLBITS, INT32_ALLBITS, INT32_ALLBITS, INT32_ALLBITS, INT32_ALLBITS);
				}
				break;
			case 4:
				switch(to) {
					case 4:
						mask = _mm256_setr_epi32(0, 0, 0, 0, INT32_ALLBITS, 0, 0, 0);
						break;
					case 5:
						mask = _mm256_setr_epi32(0, 0, 0, 0, INT32_ALLBITS, INT32_ALLBITS, 0, 0);
						break;
					case 6:
						mask = _mm256_setr_epi32(0, 0, 0, 0, INT32_ALLBITS, INT32_ALLBITS, INT32_ALLBITS, 0);
						break;
					default:
						mask = _mm256_setr_epi32(0, 0, 0, 0, INT32_ALLBITS, INT32_ALLBITS, INT32_ALLBITS, INT32_ALLBITS);
				}
				break;
			case 5:
				switch(to) {
					case 5:
						mask = _mm256_setr_epi32(0, 0, 0, 0, 0, INT32_ALLBITS, 0, 0);
						break;
					case 6:
						mask = _mm256_setr_epi32(0, 0, 0, 0, 0, INT32_ALLBITS, INT32_ALLBITS, 0);
						break;
					default:
						mask = _mm256_setr_epi32(0, 0, 0, 0, 0, INT32_ALLBITS, INT32_ALLBITS, INT32_ALLBITS);
				}
				break;
			case 6:
				switch(to) {
					case 6:
						mask = _mm256_setr_epi32(0, 0, 0, 0, 0, 0, INT32_ALLBITS, 0);
						break;
					default:
						mask = _mm256_setr_epi32(0, 0, 0, 0, 0, 0, INT32_ALLBITS, INT32_ALLBITS);
				}
				break;
			default:
				mask = _mm256_setr_epi32(0, 0, 0, 0, 0, 0, 0, INT32_ALLBITS);
		}
	}
	
	return mask;
}

__m128i _mm_setmask_fromto_epi64(int from, int to)
{
	__m128i mask;
	
	if (from > to || from >= INT64_PER_M128_REG) {
		mask = _mm_set1_epi64x(INT64_ZERO);
	} else {
		switch (from) {
			case 0:
				switch(to) {
					case 0:
						mask = _mm_setr_epi64(_mm_cvtsi64_m64(INT64_ALLBITS), _mm_cvtsi64_m64(INT64_ALLBITS));
						break;
					default:
						mask = _mm_set1_epi64x(INT64_ALLBITS);
				}

				break;
			default:
				mask = _mm_setr_epi64(_mm_cvtsi64_m64(0), _mm_cvtsi64_m64(INT64_ALLBITS));
		}
	}
	
	return mask;
}

__m256i _mm256_setmask_fromto_epi64(int from, int to)
{
	__m256i mask;
	
	if (from > to || from >= INT64_PER_M256_REG) {
		mask = _mm256_set1_epi64x(0);
	} else {
		switch (from) {
			case 0:
				switch(to) {
					case 0:
						mask = _mm256_setr_epi64x(INT64_ALLBITS, 0, 0, 0);
						break;
					case 1:
						mask = _mm256_setr_epi64x(INT64_ALLBITS, INT64_ALLBITS, 0, 0);
						break;
					case 2:
						mask = _mm256_setr_epi64x(INT64_ALLBITS, INT64_ALLBITS, INT64_ALLBITS, 0);
						break;
					default:
						mask = _mm256_set1_epi64x(INT64_ALLBITS);
				}
				break;
			case 1:
				switch(to) {
					case 1:
						mask = _mm256_setr_epi64x(0, INT64_ALLBITS, 0, 0);
						break;
					case 2:
						mask = _mm256_setr_epi64x(0, INT64_ALLBITS, INT64_ALLBITS, 0);
						break;
					default:
						mask = _mm256_setr_epi64x(0, INT64_ALLBITS, INT64_ALLBITS, INT64_ALLBITS);
				}
				break;
			case 2:
				switch(to) {
					case 2:
						mask = _mm256_setr_epi64x(0, 0, INT64_ALLBITS, 0);
						break;
					default:
						mask = _mm256_setr_epi64x(0, 0, INT64_ALLBITS, INT64_ALLBITS);
				}
				break;
			case 3:
				mask = _mm256_setr_epi64x(0, 0, 0, INT64_ALLBITS);
		}
	}
	
	return mask;
}

__m128i _mm_set_mask_epi32(int cutoff)
{
	return _mm_setmask_fromto_epi32(0, cutoff - 1);
}

__m256i _mm256_set_mask_epi32(int cutoff)
{	
	return _mm256_setmask_fromto_epi32(0, cutoff - 1);
}

__m128i _mm_set_mask_epi64(int cutoff)
{
	return _mm_setmask_fromto_epi64(0, cutoff - 1);
}

__m256i _mm256_set_mask_epi64(int cutoff)
{
	return _mm256_setmask_fromto_epi64(0, cutoff - 1);
}

__m128 _mm_set_mask_ps(int cutoff)
{
	return _mm_castsi128_ps(_mm_set_mask_epi32(cutoff));
}

__m256 _mm256_set_mask_ps(int cutoff)
{
	return _mm256_castsi256_ps(_mm256_set_mask_epi32(cutoff));
}

__m128d _mm_set_mask_pd(int cutoff)
{
	return _mm_castsi128_pd(_mm_set_mask_epi64(cutoff));
}

__m256d _mm256_set_mask_pd(int cutoff)
{
	return _mm256_castsi256_pd(_mm256_set_mask_epi64(cutoff));
}

//----------------------------------------------------------------------------
// Functions for setting values of arrays.
//----------------------------------------------------------------------------

void fset_avx2(float *x, int n, float value)
{
	int k;
	int cutoff = n % FLOAT_PER_M256_REG;
	__m256 vreg = _mm256_set1_ps(value);
	__m256i mask;
	
	if (cutoff > 0) {
		mask = _mm256_set_mask_epi32(cutoff);
		_mm256_maskstore_ps(x, mask, vreg);
	}

	for (k = cutoff; k < n; k += FLOAT_PER_M256_REG) {
		_mm256_store_ps(x + k, vreg);
	}
}

void dset_avx2(double *x, int n, double value)
{
	int k;
	int cutoff = n % DOUBLE_PER_M256_REG;
	__m256d vreg = _mm256_set1_pd(value);
	__m256i mask;
	
	if (cutoff > 0) {
		mask = _mm256_set_mask_epi64(cutoff);
		_mm256_maskstore_pd(x, mask, vreg);
	}

	for (k = cutoff; k < n; k += DOUBLE_PER_M256_REG) {
		_mm256_store_pd(x + k, vreg);
	}
}

//----------------------------------------------------------------------------
// Functions for computing sums of elements in registers.
//----------------------------------------------------------------------------

float fdot_avx2(const float *x, const float *y, int n)
{
	__m256 xreg;
	__m256 yreg;
	__m256 preg;
	__m256 sreg = _mm256_set1_ps(0);
	__m256i mask;
	int i;
	int cutoff = n % FLOAT_PER_M256_REG;
	
	if (cutoff > 0) {
		mask = _mm256_set_mask_epi32(cutoff);
		
		xreg = _mm256_maskload_ps(x, mask);
		yreg = _mm256_maskload_ps(y, mask);
		preg = _mm256_mul_ps(xreg, yreg);
		sreg = _mm256_add_ps(sreg, preg);
	}
	
	for (i = cutoff; i < n; i += FLOAT_PER_M256_REG) {
		xreg = _mm256_load_ps(x + i);
		yreg = _mm256_load_ps(y + i);
		preg = _mm256_mul_ps(xreg, yreg);
		sreg = _mm256_add_ps(sreg, preg);
	}
	
	return _mm256_register_sum_ps(sreg);
}

float fdot_indexed_avx2(const float *x, const int *xindices, const float *y, int n)
{
	__m256 xreg;
	__m256 yreg;
	__m256 preg;
	__m256 sreg = _mm256_set1_ps(0);
	__m256i vindex;
	__m256i mask;
	int i;
	int cutoff = n % FLOAT_PER_M256_REG;
	
	if (cutoff > 0) {
		// printf("DEBUG: cutoff -> %d\n", cutoff);
		mask = _mm256_set_mask_epi32(cutoff);
		// printf("DEBUG: Cutoff mask -> ");
		// _mm256_print_register_epi32(mask);
		
		vindex = _mm256_maskload_epi32(xindices, mask);
		// printf("DEBUG: Cutoff vindex -> ");
		// _mm256_print_register_epi32(vindex);
		yreg = _mm256_maskload_ps(y, mask);
		xreg = _mm256_mask_i32gather_ps(sreg, x, vindex, _mm256_castsi256_ps(mask), 4);
		
		preg = _mm256_mul_ps(xreg, yreg);
		sreg = _mm256_add_ps(sreg, preg);
	}
	
	mask = _mm256_set_mask_epi32(INT32_PER_M256_REG);
	
	for (i = cutoff; i < n; i += FLOAT_PER_M256_REG) {
		vindex = _mm256_maskload_epi32(xindices + i, mask);
		// printf("DEBUG: vindex -> ");
		//_mm256_print_register_epi32(vindex);
		yreg = _mm256_load_ps(y + i);
		xreg = _mm256_i32gather_ps(x, vindex, 4);

		preg = _mm256_mul_ps(xreg, yreg);
		sreg = _mm256_add_ps(sreg, preg);
	}
	// printf("DEBUG: \n");
	
	return _mm256_register_sum_ps(sreg);
}

double ddot_avx2(const double *x, const double *y, int n)
{
	__m256d xreg;
	__m256d yreg;
	__m256d preg;
	__m256d sreg = _mm256_set1_pd(0);
	__m256i mask;
	int i;
	int cutoff = n % DOUBLE_PER_M256_REG;
	
	if (cutoff > 0) {
		mask = _mm256_set_mask_epi64(cutoff);
		
		xreg = _mm256_maskload_pd(x, mask);
		yreg = _mm256_maskload_pd(y, mask);
		preg = _mm256_mul_pd(xreg, yreg);
		sreg = _mm256_add_pd(sreg, preg);
	}
	
	for (i = cutoff; i < n; i += DOUBLE_PER_M256_REG) {
		xreg = _mm256_load_pd(x + i);
		yreg = _mm256_load_pd(y + i);
		preg = _mm256_mul_pd(xreg, yreg);
		sreg = _mm256_add_pd(sreg, preg);
	}
	
	return _mm256_register_sum_pd(sreg);
}

float ddot_indexed_avx2(const double *x, const int *xindices, const double *y, int n)
{
	__m256d xreg;
	__m256d yreg;
	__m256d preg;
	__m256d sreg = _mm256_set1_pd(0);
	__m128i vindex;
	__m256i mask;
	__m128i mask128;
	int i;
	int cutoff = n % DOUBLE_PER_M256_REG;
	
	if (cutoff > 0) {
		mask = _mm256_set_mask_epi64(cutoff);
		mask128 = _mm_set_mask_epi32(cutoff);
		
		vindex = _mm_maskload_epi32(xindices, mask128);
		xreg = _mm256_mask_i32gather_pd(sreg, x, vindex, _mm256_castsi256_pd(mask), 8);
		yreg = _mm256_maskload_pd(y, mask);

		preg = _mm256_mul_pd(xreg, yreg);
		sreg = _mm256_add_pd(sreg, preg);
	}
	
	mask128 = _mm_set_mask_epi32(INT32_PER_M128_REG);
	
	for (i = cutoff; i < n; i += DOUBLE_PER_M256_REG) {
		vindex = _mm_maskload_epi32(xindices + i, mask128);
		yreg = _mm256_load_pd(y + i);
		xreg = _mm256_i32gather_pd(x, vindex, 8);
		
		preg = _mm256_mul_pd(xreg, yreg);
		sreg = _mm256_add_pd(sreg, preg);
	}
	
	return _mm256_register_sum_pd(sreg);
}

float _mm_register_sum_ps(__m128 vreg)
{
	vreg = _mm_hadd_ps(vreg, vreg);
	vreg = _mm_hadd_ps(vreg, vreg);
	
	return _mm_cvtss_f32(vreg);
}

float _mm256_register_sum_ps(__m256 vreg)
{
	__m256i idx = _mm256_setr_epi32(0, 1, 4, 5, 2, 3, 6, 7);

	vreg = _mm256_hadd_ps(vreg, vreg);
	vreg = _mm256_permutevar8x32_ps(vreg, idx);
	vreg = _mm256_hadd_ps(vreg, vreg);
	vreg = _mm256_hadd_ps(vreg, vreg);

	return _mm256_cvtss_f32(vreg);
}

float _mm_register_sum_pd(__m128d vreg)
{
	vreg = _mm_hadd_pd(vreg, vreg);
	
	return _mm_cvtsd_f64(vreg);
}

float _mm256_register_sum_pd(__m256d vreg)
{
	vreg = _mm256_hadd_pd(vreg, vreg);
	vreg = _mm256_hadd_pd(vreg, vreg);
	
	return _mm256_cvtsd_f64(vreg);
}

int _mm_count_nonzero_ps(__m128 a)
{
	int cmask = _mm_movemask_ps(a);
	
	return _popcnt32(cmask);
}

int _mm256_count_nonzero_ps(__m256 a)
{
	int cmask = _mm256_movemask_ps(a);
	
	return _popcnt32(cmask);
}

int _mm_count_nonzero_pd(__m128d a)
{
	int cmask = _mm_movemask_pd(a);
	
	return _popcnt32(cmask);
}

int _mm256_count_nonzero_pd(__m256d a)
{
	int cmask = _mm256_movemask_pd(a);
	
	return _popcnt32(cmask);
}

//----------------------------------------------------------------------------
// Mathematical helper functions.
//----------------------------------------------------------------------------

__m256 _mm256_ipow_ps(__m256 a, int n)
{
	__m256 res;
	__m256 one;

	if (n < 0) {
		one = _mm256_set1_ps(1.0f);
		res = _mm256_ipow_ps(a, -n);
		res = _mm256_div_ps(one, res);
	} else if (n > 0) {
		res = _mm256_set1_ps(1.0f);

		for (int i = 0; i < n; i++) {
			res = _mm256_mul_ps(res, a);
		}
	}
	
	return res;
}

//----------------------------------------------------------------------------
// Functions for computing statistics of registers.
//----------------------------------------------------------------------------

float _mm_register_min_ps(__m128 a)
{
	int minidx = 0;
	float currmin = _mm_cvtss_f32(a);
	float nextmin;
	int i;
	
	nextmin = _mm_cvtss_f32(_mm_permute_ps(a, LPERM1));

	if (nextmin < currmin) {
		currmin = nextmin;
	}
	
	nextmin = _mm_cvtss_f32(_mm_permute_ps(a, LPERM2));
		
	if (nextmin < currmin) {
		currmin = nextmin;
	}
	
	nextmin = _mm_cvtss_f32(_mm_permute_ps(a, LPERM3));
		
	if (nextmin < currmin) {
		currmin = nextmin;
	}
	
	return currmin;
}

float _mm256_register_min_ps(__m256 a)
{	
	__m256i idx = _mm256_setr_epi32(0, 1, 2, 3, 4, 5, 6, 7);
	
	int minidx = 0;
	float currmin = _mm256_cvtss_f32(a);
	float nextmin;
	int i;
	
	idx = _mm256_leftperm_epi32(idx, 1);
	nextmin = _mm256_cvtss_f32(_mm256_permutevar8x32_ps(a, idx));
		
	if (nextmin < currmin) {
		currmin = nextmin;
	}
	
	idx = _mm256_leftperm_epi32(idx, 1);
	nextmin = _mm256_cvtss_f32(_mm256_permutevar8x32_ps(a, idx));
		
	if (nextmin < currmin) {
		currmin = nextmin;
	}
	
	
	idx = _mm256_leftperm_epi32(idx, 1);
	nextmin = _mm256_cvtss_f32(_mm256_permutevar8x32_ps(a, idx));
		
	if (nextmin < currmin) {
		currmin = nextmin;
	}
	
	idx = _mm256_leftperm_epi32(idx, 1);
	nextmin = _mm256_cvtss_f32(_mm256_permutevar8x32_ps(a, idx));
		
	if (nextmin < currmin) {
		currmin = nextmin;
	}
	
	idx = _mm256_leftperm_epi32(idx, 1);
	nextmin = _mm256_cvtss_f32(_mm256_permutevar8x32_ps(a, idx));
		
	if (nextmin < currmin) {
		currmin = nextmin;
	}
	
	
	idx = _mm256_leftperm_epi32(idx, 1);
	nextmin = _mm256_cvtss_f32(_mm256_permutevar8x32_ps(a, idx));
		
	if (nextmin < currmin) {
		currmin = nextmin;
	}
	
	idx = _mm256_leftperm_epi32(idx, 1);
	nextmin = _mm256_cvtss_f32(_mm256_permutevar8x32_ps(a, idx));
		
	if (nextmin < currmin) {
		currmin = nextmin;
	}
	
	return currmin;
}

double _mm_register_min_pd(__m128d a)
{
	__m128d a0 = _mm_unpacklo_pd(a, a);
	__m128d a1 = _mm_unpackhi_pd(a, a);
	__m128d min = _mm_min_pd(a0, a1);
	
	return _mm_cvtsd_f64(min);
}

double _mm256_register_min_pd(__m256d a)
{
	const int imm8[3] = {
		0x39, // [0, 1, 2, 3] -> [1, 2, 3, 0] = 0b 00 11 10 01 = 0b(0011)(1001) = 0x39
		0x4e, // [0, 1, 2, 3] -> [2, 3, 0, 1] = 0b 01 00 11 10 = 0b(0100)(1110) = 0x4e
		0x93 // [0, 1, 2, 3] -> [3, 0, 1, 2] = 0b 10 01 00 11 = 0b(1001)(0011) = 0x93
	};
	
	int minidx = 0;
	double currmin = _mm256_cvtsd_f64(a);
	double nextmin;
	int i;
	
	nextmin = _mm256_cvtsd_f64(_mm256_permute4x64_pd(a, LPERM1));
	if (nextmin < currmin) {
		currmin = nextmin;
	}
	
	nextmin = _mm256_cvtsd_f64(_mm256_permute4x64_pd(a, LPERM2));
	if (nextmin < currmin) {
		currmin = nextmin;
	}
	
	nextmin = _mm256_cvtsd_f64(_mm256_permute4x64_pd(a, LPERM3));
	if (nextmin < currmin) {
		currmin = nextmin;
	}
	
	return currmin;
}

int _mm256_masksltnnz_epi64(__m256i a, __m256i mask)
{
	int displacement = 0;
	__m256i azero;

	if (!_mm256_testz_si256(a, a)) {
		// Get the zero positions in the remainder of the register a.
		azero = _mm256_andnot_si256(a, mask);
		
		while (_mm256_extract_epi64(azero, 0)) {
			azero = _mm256_leftshift_epi64(azero, 1);
			mask = _mm256_leftshift_epi64(mask, 1);
			
			displacement++;
		}
	}
	
	return displacement;
}

__m256d _mm256_leftperm_pd(__m256d a, int nperms)
{	
	switch (nperms) {
		case 1:
			_mm256_permute4x64_pd(a, M128_LPERM_TO_IMM8(1));
			break;
		case 2:
			_mm256_permute4x64_pd(a, M128_LPERM_TO_IMM8(2));
			break;
		case 3:
			_mm256_permute4x64_pd(a, M128_LPERM_TO_IMM8(3));
	}
	
	return a;
}

__m256d _mm256_rightperm_pd(__m256d a, int nperms)
{
	switch (nperms) {
		case 1:
			_mm256_permute4x64_pd(a, M128_RPERM_TO_IMM8(1));
			break;
		case 2:
			_mm256_permute4x64_pd(a, M128_RPERM_TO_IMM8(2));
			break;
		case 3:
			_mm256_permute4x64_pd(a, M128_RPERM_TO_IMM8(3));
	}
	
	return a;
}

__m256i _mm256_leftperm_epi64(__m256i a, int nperms)
{
	switch (nperms) {
		case 1:
			_mm256_permute4x64_epi64(a, M128_RPERM_TO_IMM8(1));
			break;
		case 2:
			_mm256_permute4x64_epi64(a, M128_RPERM_TO_IMM8(2));
			break;
		case 3:
			_mm256_permute4x64_epi64(a, M128_RPERM_TO_IMM8(3));
	}

	return a;
}

__m256i _mm256_maskleftperm_epi64(__m256i a, __m256i mask)
{
	return _mm256_and_si256(_mm256_leftperm_epi64(a, 1), mask);
}

__m256i _mm256_leftshift_epi64(__m256i a, int nshifts)
{
	__m256i mask;
	
	if (nshifts > 0) {
		mask = _mm256_set_mask_epi64(INT64_PER_M256_REG - nshifts);
		a = _mm256_maskleftperm_epi64(a, mask);
	}
	
	return a;
}

__m256i _mm256_leftperm_epi32(__m256i a, int nperms)
{
	__m256i lpidx = _mm256_setr_epi32(1, 2, 3, 4, 5, 6, 7, 0);
	int p;
	
	for (p = 0; p < nperms; p++) {
		a = _mm256_permutevar8x32_epi32(a, lpidx);
	}

	return a;
}

//----------------------------------------------------------------------------
// Functions for computing the modulus of elements in registers.
//----------------------------------------------------------------------------

// Assumes b > 0.
__m128i _mm_div_epi32(__m128i a, __m128i b)
{
	__m128i quotient;
	
	return quotient;
}

__m256i _mm256_div_epi32(__m256i a, __m256i b)
{
	__m256i quotient;

	return quotient;
}

__m128i _mm_div_epi64(__m128i a, __m128i b)
{
	__m128i quotient;
	
	return quotient;
}

__m256i _mm256_div_epi64(__m256i a, __m256i b)
{
	__m256i quotient;
	
	return quotient;
}

__m128i _mm_rem_epi32(__m128i a, __m128i b)
{
	__m128i rem;

	return rem;
}

__m256i _mm256_rem_epi32(__m256i a, __m256i b)
{
	__m256i rem;

	return rem;
}

__m128i _mm_rem_epi64(__m128i a, __m128i b)
{
	__m128i rem;

	return rem;
}

__m256i _mm256_rem_epi64(__m256i a, __m256i b)
{
	__m256i rem;

	return rem;
}

__m128i _mm_mod_epi32(__m128i a, __m128i b)
{
	__m128i mod;
	
	return mod;
}

__m256i _mm256_mod_epi32(__m256i a, __m256i b)
{
	__m256i mod;
	
	return mod;
}

__m128i _mm_mod_epi64(__m128i a, __m128i b)
{
	__m128i mod;
	
	return mod;
}

__m256i _mm256_mod_epi64(__m256i a, __m256i b)
{
	__m256i mod;
	
	return mod;
}


//----------------------------------------------------------------------------
// Helper routines for permuting
//----------------------------------------------------------------------------

__m128 _mm_shift_one_left_ps(__m128 a, int index)
{
	__m128i shift_indices;
	
	switch (index) {
		case 0:
			shift_indices = _mm_setr_epi32(3, 1, 2, 0);
			break;
		case 1:
			shift_indices = _mm_setr_epi32(1, 0, 2, 3);
			break;
		case 2:
			shift_indices = _mm_setr_epi32(0, 2, 1, 3);
			break;
		case 3:
			shift_indices = _mm_setr_epi32(0, 1, 3, 2);
	}
	
	return _mm_permutevar_ps(a, shift_indices);
}

__m128 _mm_shift_one_right_ps(__m128 a, int index)
{
	__m128i shift_indices;
	
	switch (index) {
		case 0:
			shift_indices = _mm_setr_epi32(1, 0, 2, 3);
			break;
		case 1:
			shift_indices = _mm_setr_epi32(0, 2, 1, 3);
			break;
		case 2:
			shift_indices = _mm_setr_epi32(0, 1, 3, 2);
			break;
		case 3:
			shift_indices = _mm_setr_epi32(3, 1, 2, 0);
	}
	
	return _mm_permutevar_ps(a, shift_indices);
}

//----------------------------------------------------------------------------
// Helper routines for copying data.
//----------------------------------------------------------------------------

void _mm256_copy1d_epi32(int *dst, const int *src, int n)
{
	int i;
	int cutoff = n % INT32_PER_M256_REG;
	__m256i mask;
	__m256i sreg, dreg;
	
	if (cutoff > 0) {
		mask = _mm256_set_mask_epi32(cutoff);
		
		sreg = _mm256_maskload_epi32(src, mask);
		_mm256_maskstore_epi32(dst, mask, sreg);
	}
	
	mask = _mm256_set_mask_epi32(INT32_ALLBITS);
	
	for (i = cutoff; i < n; i += INT32_PER_M256_REG) {
		sreg = _mm256_maskload_epi32(src, mask);
		_mm256_maskstore_epi32(dst, mask, sreg);
	}
}

void _mm256_copy2d_epi32(int *kind, int nrows, const int *iind, const int *jind, int numi, int numj)
{
	int i, j, jdidx, jsidx;
	int icutoff = numi % INT32_PER_M256_REG;
	int jcutoff = numj % INT32_PER_M256_REG;
	__m256i jreg, kreg;
	__m256i mask;
	
#ifdef CONTIGUOUS_LOOP
	for (j = 0; j < numj; j++) {
		jdidx = j * numi;
		jsidx = jind[j] * nrows;
		jreg = _mm256_set1_epi32(jsidx);
		
		if (icutoff > 0) {
			mask = _mm256_set_mask_epi32(icutoff);
			
			kreg = _mm256_maskload_epi32(iind, mask);
			kreg = _mm256_add_epi32(kreg, jreg);
			_mm256_maskstore_epi32(kind + jdidx, mask, kreg);
		}
		
		mask = _mm256_set_mask_epi32(INT32_ALLBITS);

		for (i = icutoff; i < numi; i += INT32_PER_M256_REG) {
			kreg = _mm256_maskload_epi32(iind + i, mask);
			kreg = _mm256_add_epi32(kreg, jreg);
			_mm256_maskstore_epi32(kind + jdidx + i, mask, kreg);
		}
	}
#else
	if (icutoff > 0) {
		mask = _mm256_set_mask_epi32(icutoff);

		for (j = 0; j < numj; j++) {
			jdidx = j * numi;
			jsidx = jind[j] * nrows;
			jreg = _mm256_set1_epi32(jsidx);
				
			kreg = _mm256_maskload_epi32(iind, mask);
			kreg = _mm256_add_epi32(kreg, jreg);
			_mm256_maskstore_epi32(kind + jdidx, mask, kreg);
		}
	}
	
	mask = _mm256_set_mask_epi32(INT32_ALLBITS);

	for (j = 0; j < numj; j++) {
		jdidx = j * numi;
		jsidx = jind[j] * nrows;
		jreg = _mm256_set1_epi32(jsidx);

		for (i = icutoff; i < numi; i += INT32_PER_M256_REG) {
			kreg = _mm256_maskload_epi32(iind + i, mask);
			kreg = _mm256_add_epi32(kreg, jreg);
			_mm256_maskstore_epi32(kind + jdidx + i, mask, kreg);
		}
	}
#endif
}


void _mm256_copy1d_ps(float *dst, const float *src, int n)
{
	int i;
	int cutoff = n % FLOAT_PER_M256_REG;
	__m256i mask;
	__m256 sreg, dreg;
	
	if (cutoff > 0) {
		mask = _mm256_set_mask_epi32(cutoff);
		
		sreg = _mm256_maskload_ps(src, mask);
		_mm256_maskstore_ps(dst, mask, sreg);
	}
	
	mask = _mm256_set_mask_epi32(INT32_ALLBITS);
	
	for (i = cutoff; i < n; i += INT32_PER_M256_REG) {
		sreg = _mm256_load_ps(src);
		_mm256_store_ps(dst, sreg);
	}
}

void _mm256_copy2d_indexed_ps(float *dst, const float *src, int nrows, const int *iind, const int *jind, int numi, int numj)
{
	int i, j, jdidx, jsidx;
	int icutoff = numi % FLOAT_PER_M256_REG;
	__m256i ireg, jreg, kreg;
	__m256i mask;
	__m256 sreg;
	__m256 zero = _mm256_set1_ps(0);
	
#ifdef CONTIGUOUS_LOOP
	for (j = 0; j < numj; j++) {
		jdidx = j * numi;
		jsidx = jind[j] * nrows;
		jreg = _mm256_set1_epi32(jsidx);
		
		if (icutoff > 0) {
			mask = _mm256_set_mask_epi32(icutoff);
			
			ireg = _mm256_maskload_epi32(iind, mask);
			kreg = _mm256_add_epi32(jreg, ireg);
			
			// Gather the data from the source.
			sreg = _mm256_mask_i32gather_ps(zero, src, kreg, _mm256_castsi256_ps(mask), 4);
			
			// Store the data.
			_mm256_maskstore_ps(dst + jdidx, mask, sreg);
		}
		
		mask = _mm256_set_mask_epi32(INT32_ALLBITS);

		for (i = icutoff; i < numi; i += INT32_PER_M256_REG) {
			ireg = _mm256_maskload_epi32(iind + i, mask);
			kreg = _mm256_add_epi32(jreg, ireg);
			
			// Gather the data from the source.
			sreg = _mm256_i32gather_ps(src, kreg, 4);
			
			// Store the data.
			_mm256_store_ps(dst + jdidx + i, sreg);
		}
	}
#else
	if (icutoff > 0) {
		mask = _mm256_set_mask_epi32(icutoff);

		for (j = 0; j < numj; j++) {
			jdidx = j * numi;
			jsidx = jind[j] * nrows;
			jreg = _mm256_set1_epi32(jsidx);
				
			ireg = _mm256_maskload_epi32(iind, mask);
			kreg = _mm256_add_epi32(jreg, ireg);
			sreg = _mm256_mask_i32gather_ps(zero, src, kreg, _mm256_castsi256_ps(mask), 4);
			_mm256_maskstore_ps(dst + jdidx, mask, sreg);
		}
	}
	
	mask = _mm256_set_mask_epi32(INT32_ALLBITS);

	for (j = 0; j < numj; j++) {
		jdidx = j * numi;
		jsidx = jind[j] * nrows;
		jreg = _mm256_set1_epi32(jsidx);

		for (i = icutoff; i < numi; i += INT32_PER_M256_REG) {
			ireg = _mm256_maskload_epi32(iind + i, mask);
			kreg = _mm256_add_epi32(jreg, ireg);
			sreg = _mm256_i32gather_ps(src, kreg, 4);
			_mm256_store_ps(dst + jdidx + i, sreg);
		}
	}
#endif
}

//----------------------------------------------------------------------------
// Helper routines for printing.
//----------------------------------------------------------------------------

void _mm256_print_register_epi32(__m256i a)
{
	int b[INT32_PER_M256_REG];
	__m256i mask = _mm256_set_mask_epi32(INT32_ALLBITS);
	
	_mm256_maskstore_epi32(b, mask, a);
	printf("%d %d %d %d %d %d %d %d\n", b[0], b[1], b[2], b[3], b[4], b[5], b[6], b[7]);
}

void _mm256_print_register_ps(__m256 a)
{
	float b[FLOAT_PER_M256_REG];
	
	_mm256_store_ps(b, a);
	printf("%f %f %f %f %f %f %f %f\n", b[0], b[1], b[2], b[3], b[4], b[5], b[6], b[7]);
}