#include <stdio.h>
#include "matrix.c"


void Swap(float *a, float *b) {
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
	// and write another function to heap allocate for small n (don't know if this is easy/possible in C)
	int *indxc = malloc(n * sizeof(int));
	int *indxr = malloc(n * sizeof(int));
	int *ipiv = calloc(n, sizeof(int));
	for (size_t i = 0; i < n; i++) {
		big = 0;
		for (size_t j = 0; j < n; j++) {
			if (ipiv[j] != 1) {
				for (size_t k = 0; k < n; k++) {
					if (ipiv[k] == 0) {
						if (ABS(MatrixGetItem(a, j, k)) >= big) {
							big = ABS(MatrixGetItem(a, j, k));
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
				Swap(MatrixGetAddr(a, irow, l), MatrixGetAddr(a, icol, l));
			}
			if (b) {
				for (size_t l = 0; l < b->cols; l++) {
					Swap(MatrixGetAddr(b, irow, l), MatrixGetAddr(b, icol, l));
				}
			}
		}
		indxr[i] = irow;
		indxc[i] = icol;
		if (MATRIX_ITEM(a, icol, icol) == 0) {
			perror("GaussJordan: Singular Matrix");
			return 1;
		}
		pivinv = 1 / MATRIX_ITEM(a, icol, icol);
		MATRIX_ITEM(a, icol, icol) = 1;
		// todo(AION): check that this is stuff is properly vectorized by the compiler
		for (size_t l = 0; l < n; l++) {
			MATRIX_ITEM(a, icol, l) *= pivinv;
		}
		if (b) {
			for (size_t l = 0; l < b->cols; l++) {
				MATRIX_ITEM(b, icol, l) *= pivinv;
			}
		}
		for (size_t ll = 0; ll < n; ll++) {
			if (ll != (size_t)icol) {
				float dum = MATRIX_ITEM(a, ll, icol);
				MATRIX_ITEM(a, ll, icol) = 0;
				for (size_t l = 0; l < n; l++) {
					MATRIX_ITEM(a, ll, l) -= MATRIX_ITEM(a, icol, l) * dum;
				}
				if (b) {
					for (size_t l = 0; l < b->cols; l++) {
						MATRIX_ITEM(b, ll, l) -= MATRIX_ITEM(b, icol, l) * dum;
					}
				}
			}
		}
	}
	for (int l = n - 1; l >= 0; l--) {
		if (indxr[l] != indxc[l]) {
			for (size_t k = 0; k < n; k++) {
				Swap(MatrixGetAddr(a, k, indxr[l]), MatrixGetAddr(a, k, indxc[l]));
			}
		}
	}
	free(ipiv);
	free(indxr);
	free(indxc);
	return 0;
}


/*
LU decomposition of input matrix A. Decomposition will be done inplace. Transcribed from "Numerical recipies in C".

As above, we use our own matrix struct, not their pointers to pointers stuff because I think that's stupid.
*/
int LUFactorize(matrix_f32 *a, size_t *index) {
	assert(a->cols == a->rows);
	size_t n = a->cols;
	float *scaling_values = malloc(n * sizeof(float));  // what we scaled each row by (1 / pivot)
	size_t imax = 0;  // row where we found the latest pivot element
	float temp;  // placeholder for possible max values
	for (size_t i = 0; i < n; i++) {
		float largest_element = 0;
		for (size_t j = 0; j < n; j++) {
			if ((temp = ABS(MATRIX_ITEM(a, i, j))) > largest_element) {
				largest_element = temp;
			}
		}
		if (largest_element == 0) {
			perror("Singular matrix in LUFactorize");
			return 1;
		}
		scaling_values[i] = 1 / largest_element;
	}

	for (size_t j = 0; j < n; j++) {
		for (size_t i = 0; i < j; i++) {
			float sum = MATRIX_ITEM(a, i, j);
			for (size_t k = 0; k < i; k++) {
				sum -= MATRIX_ITEM(a, i, k) * MATRIX_ITEM(a, k, j);
			}
			MATRIX_ITEM(a, i, j) = sum;
		}
		float largest_element = 0;
		for (size_t i = j; i < n; i++) {
			float sum = MATRIX_ITEM(a, i, j);
			for (size_t k = 0; k < j; k++) {
				sum -= MATRIX_ITEM(a, i, k) * MATRIX_ITEM(a, k, j);
			}
			MATRIX_ITEM(a, i, j) = sum;

			float temp;
			if ((temp = scaling_values[i] * ABS(sum)) >= largest_element) {
				largest_element = temp;
				imax = i;
			}
		}
		if (j != imax) {
			// todo(AION): maybe this can be done parallelized? Not sure
			// how SIMD would work for reassignment
			for (size_t k = 0; k < n; k++) {
				Swap(MatrixGetAddr(a, imax, k), MatrixGetAddr(a, j, k));
			}
			scaling_values[imax] = scaling_values[j];
		}
		index[j] = imax;
		// In the original version, they have an extra step to floor a[i][i] to some minimal value here. 
		// if (ABS(MATRIX_ITEM(a, j, j)) < MIN_FLOAT) {
		// 	MATRIX_ITEM(a, j, j) = SIGN(MATRIX_ITEM(a, j, j)) * MIN_FLOAT;
		// }
		if (j < n - 1) {
			temp = 1 / MATRIX_ITEM(a, j, j);
			for (size_t i = j + 1; i < n; i++) {
				MATRIX_ITEM(a, i, j) *= temp;
			}
		}
	}
	free(scaling_values);
	return 0;
}



/*
Once we have LU decomposition, we can use forward/backward substitution to solve systems.

Adapted from "Numerical Recipies in C", again with our own matrix class and without passing around
all the extra data about what rows where reversed.
We assume matrix A is passed in already as an LU decomposition (in compact form), and b is viewwed as a matrix
of column vectors where we want to solve (inplace) the equation Ax = b

*/

int LUBackSub(matrix_f32 *a, matrix_f32 *b) {
	assert(a->rows == a->cols);
	size_t n = a->rows;
	for (size_t j = 0; j < b->cols; j++) {
		size_t nonzero_idx = 0; // keep track of the first row with nonzero element of column b[:][j]
		for (size_t i = 0; i < n; i++) {
			float sum = MATRIX_ITEM(b, i, j);
			if (nonzero_idx) {
				for (size_t k = nonzero_idx; k < n; k++) {
					sum -= MATRIX_ITEM(a, i, k) * MATRIX_ITEM(b, k, j);
				}
			} else if (sum) {
				nonzero_idx = i;
			}
		}
		for (size_t i = n - 1; i > 0; i--) {
			float sum = MATRIX_ITEM(b, i, j);
			for (size_t k = i + 1; k < n; k++) {
				sum -= MATRIX_ITEM(a, i, k) * MATRIX_ITEM(b, k, j);
			}
			MATRIX_ITEM(b, i, j) = sum / MATRIX_ITEM(a, i, i);
		}
	}
	return 0;
}
