#ifndef INTRINSICS_UTILS
#define INTRINSICS_UTILS

#include <immintrin.h>
#include <stdint.h>

#define DOUBLE_PER_M128_REG 2
#define INT64_PER_M128_REG DOUBLE_PER_M128_REG
#define DOUBLE_PER_M256_REG 4
#define INT64_PER_M256_REG DOUBLE_PER_M256_REG

#define FLOAT_PER_M128_REG 4
#define INT32_PER_M128_REG FLOAT_PER_M128_REG
#define FLOAT_PER_M256_REG 8
#define INT32_PER_M256_REG FLOAT_PER_M256_REG

#define INT32_ZERO ((int32_t)0)
#define INT32_ALLBITS ((int32_t)0xFFFFFFFF)
#define INT32_LOWBIT  ((int32_t)0x00000001)
#define INT32_HIGHBIT ((int32_t)0x80000000)

#define INT64_ZERO    ((int64_t)0)
#define INT64_LOWBIT  ((int64_t)0x0000000000000001)
#define INT64_HIGHBIT ((int64_t)0x8000000000000000)
#define INT64_ALLBITS ((int64_t)0xFFFFFFFFFFFFFFFF)


//----------------------------------------------------------------------------
// Functions for creating masks.
//----------------------------------------------------------------------------

__m128i _mm_setmask_fromto_epi32(int, int);
__m256i _mm256_setmask_fromto_epi32(int, int);

__m128i _mm_setmask_fromto_epi64(int, int);
__m256i _mm256_setmask_fromto_epi64(int, int);

__m128i _mm_set_mask_epi32(int);
__m256i _mm256_set_mask_epi32(int);

__m128i _mm_set_mask_epi64(int);
__m256i _mm256_set_mask_epi64(int);

__m128 _mm_set_mask_ps(int);
__m256 _mm256_set_mask_ps(int);

__m128d _mm_set_mask_pd(int);
__m256d _mm256_set_mask_pd(int);

//----------------------------------------------------------------------------
// Functions for setting values of arrays.
//----------------------------------------------------------------------------

void fset_avx2(float *, int, float);
void dset_avx2(double *, int, double);

//----------------------------------------------------------------------------
// Functions for computing sums of elements in registers.
//----------------------------------------------------------------------------

float fdot_avx2(const float *, const float *, int);
float fdot_indexed_avx2(const float *, const int *, const float *, int);
double ddot_avx2(const double *, const double *, int);
float ddot_indexed_avx2(const double *, const int *, const double *, int);

float _mm_register_sum_ps(__m128);
float _mm256_register_sum_ps(__m256);
float _mm_register_sum_pd(__m128d);
float _mm256_register_sum_pd(__m256d);

int _mm_count_nonzero_ps(__m128);
int _mm256_count_nonzero_ps(__m256);
int _mm_count_nonzero_pd(__m128d);
int _mm256_count_nonzero_pd(__m256d);

//----------------------------------------------------------------------------
// Mathematical helper functions.
//----------------------------------------------------------------------------

__m256 _mm256_ipow_ps(__m256, int);

//----------------------------------------------------------------------------
// Functions for computing statistics of registers.
//----------------------------------------------------------------------------

float _mm_register_min_ps(__m128);
float _mm256_register_min_ps(__m256);
double _mm_register_min_pd(__m128d);
double _mm256_register_min_pd(__m256d);
int _mm256_masksltnnz_epi64(__m256i, __m256i);
__m256d _mm256_leftperm_pd(__m256d, int);
__m256d _mm256_rightperm_pd(__m256d, int);
__m256i _mm256_leftperm_epi64(__m256i, int);
__m256i _mm256_maskleftperm_epi64(__m256i, __m256i);
__m256i _mm256_leftshift_epi64(__m256i, int);
__m256i _mm256_leftperm_epi32(__m256i, int);

//----------------------------------------------------------------------------
// Functions for computing division, remainders, and moduli of signed
// integers.
//----------------------------------------------------------------------------

__m128i _mm_div_epi32(__m128i, __m128i);
__m256i _mm256_div_epi32(__m256i, __m256i);
__m128i _mm_div_epi64(__m128i, __m128i);
__m256i _mm256_div_epi64(__m256i, __m256i);

__m128i _mm_rem_epi32(__m128i, __m128i);
__m256i _mm256_rem_epi32(__m256i, __m256i);
__m128i _mm_rem_epi64(__m128i, __m128i);
__m256i _mm256_rem_epi64(__m256i, __m256i);

__m128i _mm_mod_epi32(__m128i, __m128i);
__m256i _mm256_mod_epi32(__m256i, __m256i);
__m128i _mm_mod_epi64(__m128i, __m128i);
__m256i _mm256_mod_epi64(__m256i, __m256i);

//----------------------------------------------------------------------------
// Helper routines for permuting
//----------------------------------------------------------------------------

__m128 _mm_shift_one_left_ps(__m128, int);
__m128 _mm_shift_one_right_ps(__m128, int);

//----------------------------------------------------------------------------
// Helper routines for copying data.
//----------------------------------------------------------------------------

void _mm256_copy1d_epi32(int *, const int *, int);
void _mm256_copy2d_epi32(int *, int, const int *, const int *, int, int);

void _mm256_copy1d_ps(float *, const float *, int);
void _mm256_copy2d_indexed_ps(float *, const float *, int, const int *, const int *, int, int);

//----------------------------------------------------------------------------
// Helper routines for printing.
//----------------------------------------------------------------------------

void _mm256_print_register_epi32(__m256i);
void _mm256_print_register_ps(__m256);

#endif