#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define MASK_BITS_PER_BYTE 8
#define MASK_BYTE_SHIFT 3

struct mask_s
{
	char *bytes;
	int nbytes;
	int nbits;
	int rows, cols;
};

struct index_map_s
{
	int *indices;
	int length;
	int rows, cols;
};

int index_map_init(struct index_map_s *imap, int rows, int cols)
{
	imap->indices = calloc(rows * cols, sizeof(int));
	
	if (imap->indices) {
		imap->length = 0;
		imap->stride = stride;
		imap->rows = rows;
		imap->cols = cols;
	}
	
	return imap->indices != NULL;
}

// Conversion between maps and index masks.
int mask_to_index_map(struct index_map_s *imap, const struct mask_s *mask)
{
	int init_success = index_map_init(imap, mask->rows, mask->cols);
	int stride;
	int length = 0;
	int b, byte, bit, nbytes;
	int mpos = 0;
	int *indices;
	int i, j, k;

	if (init_success) {
		nbytes = imap->nbytes;
		stride = imap->stride;
		indices = imap->indices;
		
		for (b = 0; b < mask->nbytes; b++) {
			byte = mask[b];

			for (bit = 0; bit < MASK_BITS_PER_BYTE; bit++) {
				mask_pos_to_array_ij(mask, byte, bit, &i, &j);
				
				indices[mpos] = j * imap->nrows + i;

				mpos += stride;
			}
			
			length += stride;
		}
		
		imap->length = length;
	}
	
	return init_success;
}

// Initialization.
int mask_init(struct mask_s *mask, int rows, int cols)
{
	int len = (rows * cols) / MASK_BITS_PER_CHAR + 1;
	mask->bits = calloc(len, sizeof(char));
	
	if (mask->bits) {
		mask->len = len;
		mask->rows = rows;
		mask->cols = cols;
	}
	
	return mask->bits != NULL;
}

// Conversion between (i,j) indices in an array and 
void mask_array_pos_to_mask_pos(const struct mask_s *mask, int index, int *bytepos, char *bit, char *bitmask)
{
	*bytepos = index >> MASK_BYTE_SHIFT;
	*bit = index - byte << MASK_BYTE_SHIFT;
	*bitmask = 1 << *bit;
}

void mask_array_ij_to_mask_pos(const struct mask_s *mask, int i, int j, int *bytepos, char *bit)
{
	int index = j * mask->rows + i;
	*bytepos = index >> MASK_BYTE_SHIFT;
	*bit = index - *bytepos << MASK_BYTE_SHIFT;
}


void mask_pos_to_array_pos(const struct mask_s *mask, int bytepos, char bit)
{
	int index = bytepos << MASK_BYTE_SHIFT + bit
	*j = index % mask->rows;
	*i = index - *j * mask->nrows;
}

//
// Operations on masks.
//

void mask_and(struct mask_s *c, const struct mask_s *a, const struct mask_s *b)
{
	int nbytes = c->nbytes;
	char *ab = a->bytes;
	char *bb = b->bytes;
	char *cb = c->bytes;
	int byte;
	
	for (byte = 0; byte < nbytes; byte++) {
		cb[byte] = ab[byte] & bb[byte];
	}
}

void mask_and_eq(struct mask_s *a, const struct mask_s *b)
{
	int nbytes = a->nbytes;
	char *ab = a->bytes;
	char *bb = b->bytes;
	int byte;
	
	for (byte = 0; byte < nbytes; byte++) {
		ab[byte] &= bb[byte];
	}
}

void mask_or(struct mask_s *c, const struct mask_s *a, const struct mask_s *b)
{
	int nbytes = c->nbytes;
	char *ab = a->bytes;
	char *bb = b->bytes;
	char *cb = c->bytes;
	int byte;
	
	for (byte = 0; byte < nbytes; byte++) {
		cb[byte] = ab[byte] | bb[byte];
	}
}

void mask_or_eq(struct mask_s *a, const struct mask_s *b)
{
	int nbytes = a->nbytes;
	char *ab = a->bytes;
	char *bb = b->bytes;
	int byte;
	
	for (byte = 0; byte < nbytes; byte++) {
		ab[byte] |= bb[byte];
	}
}

void maskzero(char *mask, int rows, int cols)
{
	memset(mask, 0, nrows * ncols * sizeof(char));
}

void maskand(char *C, const char *A, const char *B, int rows, int cols)
{
	int N = rows * cols;
	
	for (int k = 0; k < N; k++) {
		C[k] = A[k] & B[k];
	}
}

void maskor(char *C, const char *A, const char *B, int rows, int cols)
{
	int N = rows * cols;
	
	for (int k = 0; k < N; k++) {
		C[k] = A[k] | B[k];
	}
}

double masklt(char mask *C, const double *A, double val, int rows, int cols, )
{
	int N = rows * cols;
	
	for (int k = 0; k < N; k++) {
		C[k] = A[k] < val;
	}
}

double maskgt(char mask *C, const double *A, double val, int rows, int cols, )
{
	int N = rows * cols;
	
	for (int k = 0; k < N; k++) {
		C[k] = A[k] > val;
	}
}

double maskle(char mask *C, const double *A, double val, int rows, int cols, )
{
	int N = rows * cols;
	
	for (int k = 0; k < N; k++) {
		C[k] = A[k] <= val;
	}
}

double maskge(char mask *C, const double *A, double val, int rows, int cols, )
{
	int N = rows * cols;
	
	for (int k = 0; k < N; k++) {
		C[k] = A[k] >= val;
	}
}

double arrayinit(double *A, int rows, int cols)
{
	int N = nrows * ncols;
	
	for (int k = 0; k < N; k++) {
		A[k] = 
	}
}

double masksum(const double *A, const char mask *C, int rows, int cols)
{
	int N = rows * cols;
	double sum = 0;
	int k = 0;
	
	while (k < N) {
		if (C[k]) {
			sum += A[
		}
		k++;
	}
	
	return sum;
}

double arrayprod(double *P, const double *A, const double *B, int rows, int cols)
{
	int N = nrows * ncols;
	
	for (int k = 0; k < N; k++) {
		P[k] = A[k] * B[k];
	}
	
	return sum;
}

double arraysumprod(const double *A, const double *B, int rows, int cols)
{
	double sum = 0;
	
	for (int k = 0; k < N; k++) {
		sum += A[k] * B[k];
	}
	
	return sum;
}

void arraymini(int *A, int nrows, int ncols)
{
	int minval = A[0], currval;
	int N = nrows * ncols;
	
	for (int k = 1; k < N; k++) {
		currval = A[k];
		
		if (currval < minval) {
			minval = currval;
		}
	}
	
	return minval;
}

void arraymaxi(int *A, int nrows, int ncols)
{
	int maxval = A[0], currval;
	int N = nrows * ncols;
	
	for (int k = 1; k < N; k++) {
		currval = A[k];
		
		if (currval > maxval) {
			maxval = currval;
		}
	}
	
	return maxval;
}

void modulus(int *II, int offset, int nrows, int ncols, int modulo)
{
	int i, j, k;
	int N = nrows * ncols;
	
	int q = modulo ? ncols : nrows;
	int p, q, n, r;
	
	for (k = 0; k < N; k++) {
		p = II[k] + offset + q;
		
		II[k] = p - q * (p / q);
	}
}

void meshgrid(double **XX, double **YY, const double *x, const double *y, int numx, int numy)
{
	double dx = (xn - x0) / (numx - 1);
	double *X = NULL, *Y = NULL;
	int i, j;
	double x0 = x[0], y0 = y[0];
	double xn = x[numx - 1], yn = y[numy - 1];
	double dx = (xn - x0) / (numx - 1);
	double dy = (yn - y0) / (numy - 1);
	int idx;
	double xi, yj;

	*XX = calloc(numx * numy, sizeof(double))
	*YY = calloc(numx * numy, sizeof(double))
	
	if (*XX != NULL && *YY != NULL) {
		X = *XX;
		Y = *YY;

		for (j = 0; j < numy; j++) {
			idx = j * numx;
			yj = y[j];

			for (i = 0; i < numx; i++) {
				xi = x[i];

				X[idx + i] = xi;
				Y[idx + i] = yj;
			}
		}
	} else {
		free(*XX);
		free(*YY);
		
		*XX = *YY = NULL;
	}
}

void zeros(double **A, int numx, int numy)
{
	*A = calloc(numx * numy, sizeof(double));
}

double min(double a, double b)
{
	return a < b ? a : b;
}

double max(double a, double b)
{
	return a > b ? a : b;
}