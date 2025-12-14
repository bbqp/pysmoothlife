#include <stdio.h>
#include <stdlib.h>
#include <immintrin.h>

#define FLOAT_PER_M256_REG 8
#define INT32_PER_M256_REG FLOAT_PER_M256_REG
#define INT32_ALLBITS ((int)0xFFFFFFFF)
#define INT32_HIGHBIT ((int)0x80000000)

int init_array(float *x, int n, float value)
{
	int i;
	
	for (i = 0; i < n; i++) {
		x[i] = value;
	}
}

int init_sequence_array(float *x, int n, float offset)
{
	int i;
	
	for (i = 0; i < n; i++) {
		x[i] = offset + (float)i;
	}
}

int init_random_array(float *x, int n)
{
	int i;
	
	for (i = 0; i < n; i++) {
		x[i] = (float)rand() / RAND_MAX;
	}
}

int init_indices(int *indices, int n)
{
	int i;
	
	for (i = 0; i < n; i++) {
		indices[i] = i;
	}
}

int randomize_indices(int *indices, int n)
{
	int i, j, temp;
	
	for (i = 0; i < n; i++) {
		j = rand() % n;
		
		if (i != j) {
			temp = indices[i];
			indices[i] = indices[j];
			indices[j] = temp;
		}
	}
}

float fdot(const float *x, const float *y, int n)
{
	float sum = 0;
	int i;
	
	for (i = 0; i < n; i++) {
		sum += x[i] * y[i];
	}
	
	return sum;
}

__m256i _mm256_set_mask_epi32(int cutoff)
{
	__m256i mask;
	
	switch (cutoff) {
		case 0:
			mask = _mm256_set1_epi32(0);
			break;
		case 1:
			mask = _mm256_setr_epi32(INT32_ALLBITS, 0, 0, 0, 0, 0, 0, 0);
			break;
		case 2:
			mask = _mm256_setr_epi32(INT32_ALLBITS, INT32_ALLBITS, 0, 0, 0, 0, 0, 0);
			break;
		case 3:
			mask = _mm256_setr_epi32(INT32_ALLBITS, INT32_ALLBITS, INT32_ALLBITS, 0, 0, 0, 0, 0);
			break;
		case 4:
			mask = _mm256_setr_epi32(INT32_ALLBITS, INT32_ALLBITS, INT32_ALLBITS, INT32_ALLBITS, 0, 0, 0, 0);
			break;
		case 5:
			mask = _mm256_setr_epi32(INT32_ALLBITS, INT32_ALLBITS, INT32_ALLBITS, INT32_ALLBITS, INT32_ALLBITS, 0, 0, 0);
			break;
		case 6:
			mask = _mm256_setr_epi32(INT32_ALLBITS, INT32_ALLBITS, INT32_ALLBITS, INT32_ALLBITS, INT32_ALLBITS, INT32_ALLBITS, 0, 0);
			break;
		case 7:
			mask = _mm256_setr_epi32(INT32_ALLBITS, INT32_ALLBITS, INT32_ALLBITS, INT32_ALLBITS, INT32_ALLBITS, INT32_ALLBITS, INT32_ALLBITS, 0);
			break;
		default:
			mask = _mm256_set1_epi32(INT32_ALLBITS);
	}
	
	return mask;
}

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

float _mm256_register_sum_ps(__m256 vreg)
{
	__m256i idx = _mm256_setr_epi32(0, 1, 4, 5, 2, 3, 6, 7);

	vreg = _mm256_hadd_ps(vreg, vreg);
	vreg = _mm256_permutevar8x32_ps(vreg, idx);
	vreg = _mm256_hadd_ps(vreg, vreg);
	vreg = _mm256_hadd_ps(vreg, vreg);

	return _mm256_cvtss_f32(vreg);
}

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
		
		printf("DEBUG -> fdot_avx2(), x: ");
		_mm256_print_register_ps(xreg);
		printf("DEBUG -> fdot_avx2(), y: ");
		_mm256_print_register_ps(yreg);
		printf("DEBUG -> fdot_avx2():\n");
		
		preg = _mm256_mul_ps(xreg, yreg);
		sreg = _mm256_add_ps(sreg, preg);
	}
	
	for (i = cutoff; i < n; i += FLOAT_PER_M256_REG) {
		xreg = _mm256_load_ps(x + i);
		yreg = _mm256_load_ps(y + i);
		
		printf("DEBUG -> fdot_avx2(), x: ");
		_mm256_print_register_ps(xreg);
		printf("DEBUG -> fdot_avx2(), y: ");
		_mm256_print_register_ps(yreg);
		
		preg = _mm256_mul_ps(xreg, yreg);
		sreg = _mm256_add_ps(sreg, preg);
	}
	printf("DEBUG -> fdot_avx2():\n");
	
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
		mask = _mm256_set_mask_epi32(cutoff);
		vindex = _mm256_maskload_epi32(xindices, mask);
		
		printf("DEBUG -> fdot_indexed_avx2(), vindex: ");
		_mm256_print_register_epi32(vindex);
		
		yreg = _mm256_maskload_ps(y, mask);
		xreg = _mm256_mask_i32gather_ps(sreg, x, vindex, _mm256_castsi256_ps(mask), 4);
		
		printf("DEBUG -> fdot_indexed_avx2(), x: ");
		_mm256_print_register_ps(xreg);
		printf("DEBUG -> fdot_indexed_avx2(), y: ");
		_mm256_print_register_ps(yreg);
		
		preg = _mm256_mul_ps(xreg, yreg);
		sreg = _mm256_add_ps(sreg, preg);
		printf("DEBUG -> fdot_indexed_avx2():\n");
	}
	
	mask = _mm256_set_mask_epi32(INT32_PER_M256_REG);
	
	for (i = cutoff; i < n; i += FLOAT_PER_M256_REG) {
		vindex = _mm256_maskload_epi32(xindices + i, mask);
		
		printf("DEBUG -> fdot_indexed_avx2(), vindex: ");
		_mm256_print_register_epi32(vindex);
		
		yreg = _mm256_load_ps(y + i);
		xreg = _mm256_i32gather_ps(x, vindex, 4);
		
		printf("DEBUG -> fdot_indexed_avx2(), x: ");
		_mm256_print_register_ps(xreg);
		printf("DEBUG -> fdot_indexed_avx2(), y: ");
		_mm256_print_register_ps(yreg);

		preg = _mm256_mul_ps(xreg, yreg);
		sreg = _mm256_add_ps(sreg, preg);
	}
	printf("DEBUG -> fdot_indexed_avx2():\n");
	
	return _mm256_register_sum_ps(sreg);
}

int main(int argc, char *argv[])
{
	int n = 10;
	float *x = NULL, *y = NULL;
	int *idx = NULL;
	float dot[4];
	float relerr[3];
	
	x = calloc(2*n, sizeof(float));
	idx = calloc(n, sizeof(int));
	
	if (x != NULL && idx != NULL) {
		srand(0);
		
		y = x + n;
		init_sequence_array(x, n, 1.0);
		init_array(y, n, 1);
		init_indices(idx, n);
		
		dot[0] = fdot(x, y, n);
		dot[1] = fdot_avx2(x, y, n);
		dot[2] = fdot_indexed_avx2(x, idx, y, n);
		
		randomize_indices(idx, n);
		
		printf("Random Indices:\n");
		for (int i = 0; i < n; i++) {
			printf("%d ", idx[i]);
		}
		printf("\n");
		
		dot[3] = fdot_indexed_avx2(x, idx, y, n);
		
		for (int i = 1; i < 4; i++) {
			relerr[i - 1] = (dot[i] - dot[0]) / dot[0] * 100;
		}
		
		printf("Serial Dot Product:                 %g\n", dot[0]);
		printf("AVX2 Dot Product:                   %g\n", dot[1]);
		printf("Indexed AVX2 Dot Product:           %g\n", dot[2]);
		printf("Random Indexed AVX2 Dot Product:    %g\n\n", dot[3]);
		
		printf("AVX2 Relative Error:                %+f%%\n", relerr[0]);
		printf("Indexed AVX2 Relative Error:        %+f%%\n", relerr[1]);
		printf("Random Indexed AVX2 Relative Error: %+f%%\n", relerr[2]);
	}
	
	free(x);
	free(idx);
	
	return 0;
}















