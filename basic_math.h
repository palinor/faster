#pragma once
#include "linalg.h"

#define GET_BIT(x, i) (((x) >> (i)) & 1)
#define PRINT_RES(func, ...) ((printf("Input : %f, Output: %f\n", (__VA_ARGS__), (func)(__VA_ARGS__))))

float FSqrt(float x);
float FPow(float base, int exponent);
float Exp(float exponent);
