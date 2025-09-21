
/**
 * (C) Copyright 2024 Aion Feehan. All Rights Reserved.
 *
 * This software is provided 'as-is', without any express or implied warranty. In no event will the authors
 * be held liable for any damages arising from the use of this software.
 *
 */
#include <cassert>
#include <cstdlib>
#include <cstring>
#include <cmath>

#ifndef uint
typedef unsigned int uint;
#endif

#include "arena_allocator.h"

 /**
  * Structure to represent polynomial functions: f(x) -> a_0 + a_1 * X + ... + a_n * X ^ n
  */
struct PolynomialFloat {
    /** pointer to coefficient memory */
    float *coefficients;
    /** size of coefficient memory */
    uint number_of_coefficients;
    /**
     * Evaluate the polynomial function at some point x
     */
    float operator()(float x) {
        float *this_coefficient = coefficients + number_of_coefficients - 1;
        float result = *this_coefficient;
        for (int deg_idx = number_of_coefficients - 1; deg_idx > 1; --deg_idx) {
            result = fmaf(result, x, *this_coefficient--);
        }
        result += *this_coefficient;
        return result;
    };

};

/**
 * Add 2 polynomials together.
 * We assume the result has already been allocated, to allow for freedom of use in what memory
 * regions this function is called on. Sets result = A + B in polynomial terms.
 *
 * @param result : pointer to pre-allocated memory where we store the result.
 * @param A : first polynomial to sum
 * @param B : second polynomial
 */
void PolynomialFloatAdd(PolynomialFloat *result, PolynomialFloat *A, PolynomialFloat *B);
/**
 * Multiply 2 polynomials together (algebraically).
 * Sets the result equal to the algenraic representation of A x B (we assume result has been pre-allocated).
 *
 * @param result : pointer to pre-allocated result memory block
 * @param A : left polynomial to multiply
 * @param B : right polynomial to multiply
 */
void PolynomialFloatMultiply(PolynomialFloat *result, PolynomialFloat *A, PolynomialFloat *B);
/**
 * Evaluate a polynomial at a given point x
 * This is equivalent to calling the () method on PolynomialFloat itself
 *
 * @param polynome : struct containing the polynomial coefficients to evaluate
 * @param x : point at which we evaluate the polynomial function
 */
float PolynomialFloatEvaluate(PolynomialFloat *polynome, float x);
/**
 * Returns the struct containing the formal derivative (done algebraically) of the polynomial input.
 * We assume the result has been pre-allocated.
 *
 * @param result : pointer to where we should return the value
 * @param input : polynomial to be differentiated
 */
void PolynomialFloatDerivative(PolynomialFloat *result, PolynomialFloat *input);
/**
 * Compute the formal primitive of a polynomial (done algebraically).
 * If A(X) = a_0 + a_1 * X + ... + a_n X^n, we return B(X) = constant + a_0 * X + a_1 / 2 * X^2 + ... + a_n / (n + 1) * X^n.
 * We assume the result has been pre-allocated.
 *
 * @param result : pointer to where we should return the value
 * @param input : polynomial we want to compute the primitive for
 * @param constant : constant to uniquely define the primitive. Will be the first coefficient in the result (its evaluation at 0)
 */
void PolynomialFloatPrimitive(PolynomialFloat *result, PolynomialFloat *input, float constant);

#define POLYNOMIAL_EXP_FUNCTION_SIZE 4
#define POLYNOMIAL_EXP_FUNCTION_NB_COEFS 4

/**
 * Struct to represent a function that look like \f[ f(x) = \sum_i^n P_i(x) \times e^{\alpha_i * x} \f]
 * where the \f$ P_i \f$ are polynomials and \f$ \alpha_i \f$ are constants.
 *
 * This class of functions is stable by addition, multiplication, derivitive and primitive, and we have algebraic formulas
 * for each operation. This makes it a convenient building block for representing our term structures.
 *
 * The elements in the item are always kept sorted in ascending order of the exp_coef, up to a (hardcoded) precision of 1e-9. Two coefficients
 * are assumed to be equal at this precision level, and are treated as equal when e.g. summing two functions together.
 *
 * To make memory management easier, we allocate room for 4 polynomials of 4 coefficients each.
 */
struct PolynomialExpFunctionSumFloat {
    float exp_coefs[POLYNOMIAL_EXP_FUNCTION_SIZE];
    PolynomialFloat polynomials[POLYNOMIAL_EXP_FUNCTION_SIZE];
    uint number_of_polynomials;

    float operator()(float x) {
        float result = 0;
        float *this_exp_coef = exp_coefs;
        PolynomialFloat *this_polynomial = polynomials;
        for (uint i = 0; i < number_of_polynomials; ++i) {
            result += (*this_polynomial++)(x) * exp((*this_exp_coef++) * x);
        }
        return result;
    }
};


bool ListIsSortedFloat(float *list_to_sort, uint list_length);
void InsertSortInplaceFloat(float *list_to_sort, uint list_length);
void InsertSortInplaceAndTrackFloat(float *list_to_sort, uint list_length, uint *tracking_list);
uint CountDistinctFloatSortedArrays(float *list_1, uint list_1_length, float *list_2, uint list_2_length);
/**
 * Add two functions together.
 * We assume that both A and B have their elements sorted by increasing value of exp_coef, and that result has already been allocated.
 */
void PolynomialExpFunctionSumFloatAdd(PolynomialExpFunctionSumFloat *result, PolynomialExpFunctionSumFloat *A, PolynomialExpFunctionSumFloat *B);
/**
 * Multiply 2 functions together and set the value in result.
 *
 */
void PolynomialExpFunctionSumFloatMultiply(PolynomialExpFunctionSumFloat *result, PolynomialExpFunctionSumFloat *A, PolynomialExpFunctionSumFloat *B);
void PolynomialExpFunctionSumFloatDerivative(PolynomialExpFunctionSumFloat *result, PolynomialExpFunctionSumFloat *input);
void PolynomialExpFloatSinglePrimitive(PolynomialFloat *result, PolynomialFloat *input, float exp_coef);
void PolynomialExpFunctionSumFloatPrimitive(PolynomialExpFunctionSumFloat *result, PolynomialExpFunctionSumFloat *input);
void PolynomialExpFunctionSumFloatToExponential(PolynomialExpFunctionSumFloat *result, PolynomialExpFunctionSumFloat *input);

struct PiecewiseFunctionFloat {
    uint term_structure_size;
    float *time_term_structure;
    PolynomialExpFunctionSumFloat *piecewise_functions;

    float operator()(float x) {
        float min_time = time_term_structure[0];
        float max_time = time_term_structure[term_structure_size - 1];
        assert(min_time < x);
        assert(max_time > x);
        assert(term_structure_size > 0);
        PolynomialExpFunctionSumFloat *function = piecewise_functions;
        float *this_time = time_term_structure + 1;
        while (*this_time < x) {
            ++this_time;
            ++function;
        }
        return (*function)(x);
    }
};

void PiecewiseFunctionFloatAlloc(PiecewiseFunctionFloat *result, uint term_structure_size, Arena *arena);
void PiecewiseFunctionFloatAdd(PiecewiseFunctionFloat *result, PiecewiseFunctionFloat *A, PiecewiseFunctionFloat *B);
void PiecewiseFunctionFloatMultiply(PiecewiseFunctionFloat *result, PiecewiseFunctionFloat *A, PiecewiseFunctionFloat *B);
void PiecewiseFunctionFloatAddConstant(PiecewiseFunctionFloat *result, float constant);
void PiecewiseFunctionFloatMultiplyByConstant(PiecewiseFunctionFloat *result, float constant);
void PiecewiseFunctionFloatDerivative(PiecewiseFunctionFloat *result, PiecewiseFunctionFloat *input);
void PiecewiseFunctionFloatPrimitive(PiecewiseFunctionFloat *result, PiecewiseFunctionFloat *input);
float PiecewiseFunctionFloatIntegral(PiecewiseFunctionFloat *function, float integral_start, float integral_end);
/**
 * Create a new function that looks like g: t -> exp(Integral_start^t f(s) ds)
 */
void PiecewiseFunctionFloatGetExponentialForwardIntegral(PiecewiseFunctionFloat *result, PiecewiseFunctionFloat *input, float start);
void PiecewiseFunctionFloatGetExponentialBackwardIntegral(PiecewiseFunctionFloat *result, PiecewiseFunctionFloat *input, float end);
void TestPolynomials();
void TestPolynomialExpFunctionSumFloat();
void TestConstantPiecewiseFunction();
void TestPiecewiseFunctionFloatExponential();
