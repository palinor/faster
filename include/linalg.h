#pragma once
#include "matrix.h"
void Swap(float *a, float *b);
/*
Gauss-Jordan with full pivoting, transcribed from "Numerical Recipies in C"

Let's assume we use the matrix_f32 struct, so our data is set up in memory differently
than in the book (we just have one block of memory, not an array of pointers)
*/
int GaussJordan(matrix_f32 *a, matrix_f32 *b);

/*
LU decomposition of input matrix A. Decomposition will be done inplace. Transcribed from "Numerical recipies in C".

As above, we use our own matrix struct, not their pointers to pointers stuff because I think that's stupid.
*/
int LUFactorize(matrix_f32 *a, size_t *index);

/*
Once we have LU decomposition, we can use forward/backward substitution to solve systems.

Adapted from "Numerical Recipies in C", again with our own matrix class and without passing around
all the extra data about what rows where reversed.
We assume matrix A is passed in already as an LU decomposition (in compact form), and b is viewwed as a matrix
of column vectors where we want to solve (inplace) the equation Ax = b

*/
int LUBackSub(matrix_f32 *a, matrix_f32 *b);
