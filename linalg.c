#include <stdio.h>
#include "matrix.c"

void swap(float *a, float *b) {
	float temp = *a;
	*a = *b;
	*b = temp;
}



/*
Gauss-Jordan with full pivoting, transcribed from "Numerical Recipies in C"

Let's assume we use the matrix_f32 struct, so our data is set up in memory differently
than in the book (we just have one block of memory, not an array of pointers)
*/
int GaussJordan(matrix_f32 *a, matrix_f32 *b) {
	assert(a->cols == a->rows);
	size_t n = a->rows;
	int icol = 0, irow = 0;
	float big = 0, pivinv = 0;
	//todo(AION) : for now we malloc everything, we can come back
	// and write another function to heap allocate for small n
	int *indxc = malloc(n * sizeof(int));
	int *indxr = malloc(n * sizeof(int));
	int *ipiv = calloc(n, sizeof(int));
	for (size_t i = 0; i < n; i++) {
		big = 0;
		for (size_t j = 0; j < n; j++) {
			if (ipiv[j] != 1) {
				for (size_t k = 0; k < n; k++) {
					if (ipiv[k] == 0) {
						if (fabs(MatrixGetItem(a, j, k)) >= big) {
							big = fabs(MatrixGetItem(a, j, k));
							irow = j;
							icol = k;
						}
					} else {
						if (ipiv[k] > 1) {
							perror("GaussJordan: Singular Matrix");
							return 1;
						}
					}
				}
			}
		}
		++(ipiv[icol]);
		if (irow != icol) {
			for (size_t l = 0; l < n; l++) {
				swap(MatrixGetAddr(a, irow, l), MatrixGetAddr(a, icol, l));
			}
			if (b) {
				for (size_t l = 0; l < b->cols; l++) {
					swap(MatrixGetAddr(b, irow, l), MatrixGetAddr(b, icol, l));
				}
			}
		}
		indxr[i] = irow;
		indxc[i] = icol;
		if (MATRIXITEM(a, icol, icol) == 0) {
			perror("GaussJordan: Singular Matrix");
			return 1;
		}
		pivinv = 1 / MATRIXITEM(a, icol, icol);
		MATRIXITEM(a, icol, icol) = 1;
		for (size_t l = 0; l < n; l++) {
			MATRIXITEM(a, icol, l) *= pivinv;
		}
		if (b) {
			for (size_t l = 0; l < b->cols; l++) {
				MATRIXITEM(b, icol, l) *= pivinv;
			}
		}
		for (size_t ll = 0; ll < n; ll++) {
			if (ll != (size_t)icol) {
				float dum = MATRIXITEM(a, ll, icol);
				MATRIXITEM(a, ll, icol) = 0;
				for (size_t l = 0; l < n; l++) {
					MATRIXITEM(a, ll, l) -= MATRIXITEM(a, icol, l) * dum;
				}
				if (b) {
					for (size_t l = 0; l < b->cols; l++) {
						MATRIXITEM(b, ll, l) -= MATRIXITEM(b, icol, l) * dum;
					}
				}
			}
		}
	}
	for (int l = n - 1; l >= 0; l--) {
		if (indxr[l] != indxc[l]) {
			for (size_t k = 0; k < n; k++) {
				swap(MatrixGetAddr(a, k, indxr[l]), MatrixGetAddr(a, k, indxc[l]));
			}
		}
	}
	free(ipiv);
	free(indxr);
	free(indxc);
	return 0;
}
