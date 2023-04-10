#pragma once
#include <stdlib.h>
#include "matrix.h"
int LogErrorExact(float expectedOutput, float output, float *inputs, size_t nInputs);
int LogErrorLeq(float expectedOutput, float output, float errorTol, float *inputs, size_t nInputs);
void RunTestFunction(int testF(void), char *fName);
int TestFSqrt();
int TestFPow();
void TestSwap();
int TestGaussJordan();
int TestLU();
void MatrixLUPopulate(matrix_f32 *lu, matrix_f32 *l, matrix_f32 *u);
int TestLUFactorize();

int TestsMain();
