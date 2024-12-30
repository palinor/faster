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

#include "arena_allocator.h"
#include "piecewise_functions.h"



/**
 * Add 2 polynomials together.
 * We assume the result has already been allocated, to allow for freedom of use in what memory
 * regions this function is called on. Sets result = A + B in polynomial terms.
 *
 * @param result : pointer to pre-allocated memory where we store the result.
 * @param A : first polynomial to sum
 * @param B : second polynomial
 */
void PolynomialFloatAdd(PolynomialFloat *result, PolynomialFloat *A, PolynomialFloat *B) {
    PolynomialFloat *shorter_polynomial = A->number_of_coefficients < B->number_of_coefficients ? A : B;
    PolynomialFloat *longer_polynomial = A->number_of_coefficients < B->number_of_coefficients ? B : A;
    result->number_of_coefficients = longer_polynomial->number_of_coefficients;
    float *this_coefficient = result->coefficients;
    float *shorter_coef = shorter_polynomial->coefficients;
    float *longer_coef = longer_polynomial->coefficients;
    if (!this_coefficient) return;
    for (uint coef_idx = 0; coef_idx < shorter_polynomial->number_of_coefficients; ++coef_idx) {
        *this_coefficient++ = *shorter_coef++ + *longer_coef++;
    }
    float *coefficients_to_copy = longer_polynomial->coefficients + shorter_polynomial->number_of_coefficients;
    uint copy_length = (longer_polynomial->number_of_coefficients - shorter_polynomial->number_of_coefficients) * sizeof(float);
    memcpy(this_coefficient, coefficients_to_copy, copy_length);
}

/**
 * Multiply 2 polynomials together (algebraically).
 * Sets the result equal to the algenraic representation of A x B (we assume result has been pre-allocated).
 *
 * @param result : pointer to pre-allocated result memory block
 * @param A : left polynomial to multiply
 * @param B : right polynomial to multiply
 */
void PolynomialFloatMultiply(PolynomialFloat *result, PolynomialFloat *A, PolynomialFloat *B) {
    // deg(result) = deg(A) + deg(B)
    // deg(A) = A->number_of_coefficients - 1
    if (!result->coefficients) return;
    PolynomialFloat temp_result;
    temp_result.number_of_coefficients = A->number_of_coefficients + B->number_of_coefficients - 1;
    assert (temp_result.number_of_coefficients > 0);
    if ((A->number_of_coefficients < 2) != (B->number_of_coefficients < 2)) {
        temp_result.number_of_coefficients--;
    }
    temp_result.coefficients = (float *)alloca(temp_result.number_of_coefficients * sizeof(float));
    memset(temp_result.coefficients, 0, temp_result.number_of_coefficients * sizeof(float));
    for (uint a_idx = 0; a_idx < A->number_of_coefficients; ++a_idx) {
        for (uint b_idx = 0; b_idx < B->number_of_coefficients; ++b_idx) {
            temp_result.coefficients[a_idx + b_idx] += A->coefficients[a_idx] * B->coefficients[b_idx];
        }
    }
    result->number_of_coefficients = temp_result.number_of_coefficients;
    memcpy(result->coefficients, temp_result.coefficients, result->number_of_coefficients * sizeof(float));
}

/**
 * Evaluate a polynomial at a given point x
 * This is equivalent to calling the () method on PolynomialFloat itself
 *
 * @param polynome: struct containing the polynomial coefficients to evaluate
 * @param x: point at which we evaluate the polynomial function
 */
float PolynomialFloatEvaluate(PolynomialFloat *polynome, float x) {
    float result = 0;
    float x_n = 1;
    float *this_coefficient = polynome->coefficients;
    for (uint deg_idx = 0; deg_idx < polynome->number_of_coefficients; ++deg_idx) {
        result += *this_coefficient++ * x_n;
        x_n *= x;
    }
    return result;
}

/**
 * Returns the struct containing the formal derivative (done algebraically) of the polynomial input.
 * We assume the result has been pre-allocated.
 *
 * @param result: pointer to where we should return the value
 * @param input: polynomial to be differentiated
 */
void PolynomialFloatDerivative(PolynomialFloat *result, PolynomialFloat *input) {
    if (!result->coefficients) return;
    if (input->number_of_coefficients == 1) {
        result->number_of_coefficients = 1;
        result->coefficients[0] = 0;
        return;
    }
    result->number_of_coefficients = input->number_of_coefficients - 1;
    float *this_coef = input->coefficients + 1;
    float *result_coef = result->coefficients;
    for (uint coef_idx = 1; coef_idx < input->number_of_coefficients; ++coef_idx) {
        *result_coef++ = *this_coef++ * (float)coef_idx;
    }
}

/**
 * Compute the formal primitive of a polynomial (done algebraically).
 * If A(X) = a_0 + a_1 * X + ... + a_n X^n, we return B(X) = constant + a_0 * X + a_1 / 2 * X^2 + ... + a_n / (n + 1) * X^n.
 * We assume the result has been pre-allocated.
 *
 * @param result : pointer to where we should return the value
 * @param input : polynomial we want to compute the primitive for
 * @param constant : constant to uniquely define the primitive. Will be the first coefficient in the result (its evaluation at 0)
 */
void PolynomialFloatPrimitive(PolynomialFloat *result, PolynomialFloat *input, float constant) {
    if (!result->coefficients) return;
    result->number_of_coefficients = input->number_of_coefficients + 1;
    result->coefficients[0] = constant;
    float *this_coef = input->coefficients;
    float *result_coef = result->coefficients + 1;
    for (uint coef_idx = 0; coef_idx < input->number_of_coefficients; ++coef_idx) {
        *result_coef++ = *this_coef++ / (float)(coef_idx + 1);
    }
}


bool ListIsSortedFloat(float *list_to_sort, uint list_length) {
    float *this_element = list_to_sort + 1;
    float *last_element = list_to_sort;
    bool is_sorted = *last_element < *this_element;
    uint list_idx = 0;
    while (is_sorted && (list_idx < list_length)) {
        last_element = this_element++;
        is_sorted &= (*last_element < *this_element);
        ++list_idx;
    }
    return is_sorted;
}

void InsertSortInplaceFloat(float *list_to_sort, uint list_length) {
    for (uint i = 1; i < list_length; ++i) {
        float this_value = list_to_sort[i];
        uint j = i;
        while ((list_to_sort[j] > this_value) && (j > 0)) {
            list_to_sort[j] = list_to_sort[j - 1];
            --j;
        }
        list_to_sort[j] = this_value;
    }
    return;
}


void InsertSortInplaceAndTrackFloat(float *list_to_sort, uint list_length, uint *tracking_list) {
    for (uint i = 0; i < list_length; ++i) {
        tracking_list[i] = i;
    }
    for (uint i = 1; i < list_length; ++i) {
        float this_value = list_to_sort[i];
        uint j = i;
        while ((list_to_sort[j] > this_value) && (j > 0)) {
            list_to_sort[j] = list_to_sort[j - 1];
            tracking_list[j] = tracking_list[j - 1];
            --j;
        }
    }
}


uint CountDistinctFloatSortedArrays(float *list_1, uint list_1_length, float *list_2, uint list_2_length) {
    uint distinct_count = list_1_length + list_2_length;
    float *element_1 = list_1;
    for (uint i = 0; i < list_1_length; ++i) {
        float *element_2 = list_2;
        while (*element_2 < *element_1) {
            ++element_2;
        }
        if (fabsf(*element_2 - *element_1) < 1e-9f) {
            --distinct_count;
        }
        ++element_1;
    }
    return distinct_count;
}

/**
 * Add two functions together.
 * We assume that both A and B have their elements sorted by increasing value of exp_coef, and that result has already been allocated.
 */
void PolynomialExpFunctionSumFloatAdd(PolynomialExpFunctionSumFloat *result, PolynomialExpFunctionSumFloat *A, PolynomialExpFunctionSumFloat *B) {
    uint n_distinct_exp_coefs = CountDistinctFloatSortedArrays(A->exp_coefs, A->number_of_polynomials, B->exp_coefs, B->number_of_polynomials);
    assert (n_distinct_exp_coefs <= POLYNOMIAL_EXP_FUNCTION_SIZE);
    PolynomialFloat *a_polynomial = A->polynomials;
    float *a_exp_coef = A->exp_coefs;
    PolynomialFloat *b_polynomial = B->polynomials;
    float *b_exp_coef = B->exp_coefs;

    PolynomialExpFunctionSumFloat temp_result;
    temp_result.number_of_polynomials = n_distinct_exp_coefs;
    for (uint i = 0; i < n_distinct_exp_coefs; ++i) {
        temp_result.polynomials[i].coefficients = (float *)alloca(sizeof(float) * POLYNOMIAL_EXP_FUNCTION_NB_COEFS);
        memset(temp_result.polynomials[i].coefficients, 0, sizeof(float) * POLYNOMIAL_EXP_FUNCTION_NB_COEFS);
    }
    PolynomialFloat *result_polynomial = temp_result.polynomials;
    float *result_exp_coef = temp_result.exp_coefs;
    uint a_exp_count = 0;
    uint b_exp_count = 0;
    for (uint i = 0; i < n_distinct_exp_coefs; ++i) {
        if (!(a_exp_count < A->number_of_polynomials)) {
            *result_exp_coef = *b_exp_coef++;
            result_polynomial->number_of_coefficients = b_polynomial->number_of_coefficients;
            memcpy(result_polynomial->coefficients, (b_polynomial++)->coefficients, sizeof(float) * result_polynomial->number_of_coefficients);
            ++b_exp_count;
        } else if (!(b_exp_count < B->number_of_polynomials)) {
            *result_exp_coef = *a_exp_coef++;
            result_polynomial->number_of_coefficients = a_polynomial->number_of_coefficients;
            memcpy(result_polynomial->coefficients, (a_polynomial++)->coefficients, sizeof(float) * result_polynomial->number_of_coefficients);
            ++a_exp_count;
        } else if (fabs(*a_exp_coef - *b_exp_coef) < 1e-9f) {
            PolynomialFloatAdd(result_polynomial, a_polynomial, b_polynomial);
            *result_exp_coef = *a_exp_coef;
            ++a_polynomial;
            ++a_exp_coef;
            ++b_polynomial;
            ++b_exp_coef;
            ++a_exp_count;
            ++b_exp_count;
        } else if (*a_exp_coef < *b_exp_coef) {
            *result_exp_coef = *a_exp_coef++;
            result_polynomial->number_of_coefficients = a_polynomial->number_of_coefficients;
            memcpy(result_polynomial->coefficients, (a_polynomial++)->coefficients, sizeof(float) * result_polynomial->number_of_coefficients);
            ++a_exp_count;
        } else {
            *result_exp_coef = *b_exp_coef++;
            result_polynomial->number_of_coefficients = b_polynomial->number_of_coefficients;
            memcpy(result_polynomial->coefficients, (b_polynomial++)->coefficients, sizeof(float) * result_polynomial->number_of_coefficients);
            ++b_exp_count;
        }
        ++result_exp_coef;
        ++result_polynomial;
    }
    result->number_of_polynomials = n_distinct_exp_coefs;
    memcpy(result->exp_coefs, temp_result.exp_coefs, n_distinct_exp_coefs * sizeof(float));
    PolynomialFloat *source_polynomial = temp_result.polynomials;
    PolynomialFloat *dest_polynomial = result->polynomials;
    for (uint i = 0; i < n_distinct_exp_coefs; ++i) {
        uint n_coefs = source_polynomial->number_of_coefficients;
        dest_polynomial->number_of_coefficients = n_coefs;
        memcpy((dest_polynomial++)->coefficients, (source_polynomial++)->coefficients, n_coefs * sizeof(float));
    }
}


/**
 * Multiply 2 functions together and set the value in result.
 *
 */
void PolynomialExpFunctionSumFloatMultiply(PolynomialExpFunctionSumFloat *result, PolynomialExpFunctionSumFloat *A, PolynomialExpFunctionSumFloat *B) {
    uint coef_list_size = 1;
    float scratch_exp_coef_list[2 * POLYNOMIAL_EXP_FUNCTION_SIZE];
    memset(scratch_exp_coef_list, 0, 2 * POLYNOMIAL_EXP_FUNCTION_SIZE * sizeof(float));
    assert(2 * POLYNOMIAL_EXP_FUNCTION_SIZE > A->number_of_polynomials + B->number_of_polynomials);
    float *a_coef = A->exp_coefs;
    scratch_exp_coef_list[0] = *a_coef + B->exp_coefs[0];
    for (uint i = 0; i < A->number_of_polynomials; ++i) {
        float *b_coef = B->exp_coefs;
        for (uint j = 0; j < B->number_of_polynomials; ++j) {
            bool is_duplicate = false;
            float this_coef = *a_coef + *b_coef;
            float *scratch_exp_coef = scratch_exp_coef_list;
            for (uint k = 0; k < coef_list_size; ++k) {
                if (fabsf(this_coef - *scratch_exp_coef++) < 1e-9f) {
                    is_duplicate = true;
                    break;
                }
            }
            if (!is_duplicate) {
                scratch_exp_coef_list[coef_list_size++ - 1] = this_coef;
            }
            ++b_coef;
        }
        ++a_coef;
    }

    for (uint i = 1; i < coef_list_size; ++i) {
        float this_value = scratch_exp_coef_list[i];
        uint j = i;
        while ((j > 0) && (scratch_exp_coef_list[j - 1] > this_value)) {
            scratch_exp_coef_list[j] = scratch_exp_coef_list[j - 1];
            scratch_exp_coef_list[j - 1] = this_value;
            --j;
        }
    }

    PolynomialFloat temp_polynomial_results[POLYNOMIAL_EXP_FUNCTION_SIZE];
    // init all the polynomials to 1.
    assert (coef_list_size <= POLYNOMIAL_EXP_FUNCTION_SIZE);
    for (uint i = 0; i < coef_list_size; ++i) {
        temp_polynomial_results[i].coefficients = (float *)alloca(POLYNOMIAL_EXP_FUNCTION_NB_COEFS * sizeof(float));
        memset(temp_polynomial_results[i].coefficients, 0, POLYNOMIAL_EXP_FUNCTION_NB_COEFS * sizeof(float));
        temp_polynomial_results[i].coefficients[0] = 1;
        temp_polynomial_results[i].number_of_coefficients = 1;
    }

    for (uint i = 0; i < A->number_of_polynomials; ++i) {
        for (uint j = 0; j < B->number_of_polynomials; ++j) {
            float exp_coef = A->exp_coefs[i] + B->exp_coefs[j];
            uint exp_coef_index = 0;
            float *this_exp_coef = scratch_exp_coef_list;
            while (fabsf(exp_coef - *this_exp_coef++) > 1e-9f) {
                ++exp_coef_index;
            }
            PolynomialFloatMultiply(temp_polynomial_results + exp_coef_index, temp_polynomial_results + exp_coef_index, A->polynomials + i);
            PolynomialFloatMultiply(temp_polynomial_results + exp_coef_index, temp_polynomial_results + exp_coef_index, B->polynomials + j);
        }
    }
    for (uint i = 0; i < coef_list_size; ++i) {
        result->polynomials[i].number_of_coefficients = temp_polynomial_results[i].number_of_coefficients;
        memcpy(result->polynomials[i].coefficients, temp_polynomial_results[i].coefficients, POLYNOMIAL_EXP_FUNCTION_NB_COEFS * sizeof(float));
    }
    memcpy(result->exp_coefs, scratch_exp_coef_list, coef_list_size * sizeof(float));
    result->number_of_polynomials = coef_list_size;
}


void PolynomialExpFunctionSumFloatDerivative(PolynomialExpFunctionSumFloat *result, PolynomialExpFunctionSumFloat *input) {
    PolynomialFloat derivative;
    PolynomialFloat exp_product;
    float coefficient_list_1[POLYNOMIAL_EXP_FUNCTION_NB_COEFS];
    float coefficient_list_2[POLYNOMIAL_EXP_FUNCTION_NB_COEFS];
    derivative.coefficients = coefficient_list_1;
    exp_product.coefficients = coefficient_list_2;

    PolynomialFloat *this_polynomial = input->polynomials;
    float *this_exp_coef = input->exp_coefs;
    for (uint i = 0; i < input->number_of_polynomials; ++i) {
        PolynomialFloatDerivative(&derivative, this_polynomial);
        float *this_coef = this_polynomial->coefficients;
        exp_product.number_of_coefficients = this_polynomial->number_of_coefficients;
        for (uint j = 0; j < this_polynomial->number_of_coefficients; ++j) {
            exp_product.coefficients[j] = (*this_exp_coef) * (*this_coef++);
        }
        PolynomialFloatAdd(result->polynomials + i, &derivative, &exp_product);
    }
    memcpy(result->exp_coefs, input->exp_coefs, input->number_of_polynomials * sizeof(float));
    result->number_of_polynomials = input->number_of_polynomials;
}

void PolynomialExpFloatSinglePrimitive(PolynomialFloat *result, PolynomialFloat *input, float exp_coef) {
    if (fabsf(exp_coef) < 1e-9) {
        PolynomialFloatPrimitive(result, input, 0);
        return;
    }
    if (input->number_of_coefficients == 1) {
        result->number_of_coefficients = 1;
        result->coefficients[0] = input->coefficients[0] / exp_coef;
        return;
    }
    PolynomialFloat temp_1;
    PolynomialFloat temp_2;
    float temp_coefs_1[POLYNOMIAL_EXP_FUNCTION_NB_COEFS];
    float temp_coefs_2[POLYNOMIAL_EXP_FUNCTION_NB_COEFS];
    memset(temp_coefs_1, 0, POLYNOMIAL_EXP_FUNCTION_NB_COEFS * sizeof(float));
    memset(temp_coefs_2, 0, POLYNOMIAL_EXP_FUNCTION_NB_COEFS * sizeof(float));
    temp_1.coefficients = temp_coefs_1;
    temp_2.coefficients = temp_coefs_2;
    temp_1.number_of_coefficients = input->number_of_coefficients;
    PolynomialFloatDerivative(&temp_1, input);
    for (uint i = 0; i < input->number_of_coefficients; ++i) {
        temp_1.coefficients[i] /= exp_coef;
    }
    PolynomialExpFloatSinglePrimitive(&temp_2, &temp_1, exp_coef);

    for (uint i = 0; i < input->number_of_coefficients; ++i) {
        result->coefficients[i] += input->coefficients[i] / exp_coef - temp_2.coefficients[i];
    }
    result->number_of_coefficients = input->number_of_coefficients;

}

void PolynomialExpFunctionSumFloatPrimitive(PolynomialExpFunctionSumFloat *result, PolynomialExpFunctionSumFloat *input) {
    result->number_of_polynomials = input->number_of_polynomials;
    memcpy(result->exp_coefs, input->exp_coefs, input->number_of_polynomials * sizeof(float));
    PolynomialFloat *result_polynomial = result->polynomials;
    PolynomialFloat *input_polynomial = input->polynomials;
    float *exp_coef = input->exp_coefs;
    for (uint i = 0; i < input->number_of_polynomials; ++i) {
        PolynomialExpFloatSinglePrimitive(result_polynomial, input_polynomial, *exp_coef);
        ++result_polynomial;
        ++input_polynomial;
        ++exp_coef;
    }
}

void PolynomialExpFunctionSumFloatToExponential(PolynomialExpFunctionSumFloat *result, PolynomialExpFunctionSumFloat *input) {
    assert(input->number_of_polynomials == 1);
    assert(fabsf(input->exp_coefs[0]) < 1e-9);
    PolynomialFloat input_polynomial = input->polynomials[0];
    assert(input_polynomial.number_of_coefficients < 3);
    result->number_of_polynomials = 1;
    result->polynomials[0].coefficients[0] = expf(input_polynomial.coefficients[0]);
    result->polynomials[0].number_of_coefficients = 1;
    if (input_polynomial.number_of_coefficients == 2) {
        result->exp_coefs[0] = input_polynomial.coefficients[1];
    }
    result->polynomials[0].coefficients[1] = 0;
    result->polynomials[0].coefficients[2] = 0;
    result->polynomials[0].coefficients[3] = 0;
}


void PiecewiseFunctionFloatAlloc(PiecewiseFunctionFloat *result, uint term_structure_size, Arena *arena) {
    result->term_structure_size = term_structure_size;
    result->time_term_structure = (float *)ArenaGetMemory(term_structure_size * sizeof(float), arena);
    result->piecewise_functions = (PolynomialExpFunctionSumFloat *)ArenaGetMemory(term_structure_size * sizeof(PolynomialExpFunctionSumFloat), arena);
    PolynomialExpFunctionSumFloat *this_function = result->piecewise_functions;
    for (uint i = 0; i < term_structure_size; ++i) {
        for (uint j = 0; j < POLYNOMIAL_EXP_FUNCTION_SIZE; ++j) {
            this_function->polynomials[j].coefficients = (float *)ArenaGetMemory(POLYNOMIAL_EXP_FUNCTION_NB_COEFS * sizeof(float), arena);
        }
        ++this_function;
    }
}

void PiecewiseFunctionFloatAdd(PiecewiseFunctionFloat *result, PiecewiseFunctionFloat *A, PiecewiseFunctionFloat *B) {
    assert(A->term_structure_size == B->term_structure_size);
    uint term_structure_size = A->term_structure_size;
    float *a_time = A->time_term_structure;
    float *b_time = B->time_term_structure;
    PolynomialExpFunctionSumFloat *result_function = result->piecewise_functions;
    PolynomialExpFunctionSumFloat *a_function = A->piecewise_functions;
    PolynomialExpFunctionSumFloat *b_function = B->piecewise_functions;

    for (uint i = 0; i < term_structure_size; ++i) {
        assert(fabsf(*a_time++ - *b_time++) < 1e-9f);
        PolynomialExpFunctionSumFloatAdd(result_function++, a_function++, b_function++);
    }
    result->term_structure_size = term_structure_size;
}

void PiecewiseFunctionFloatMultiply(PiecewiseFunctionFloat *result, PiecewiseFunctionFloat *A, PiecewiseFunctionFloat *B) {
    assert(A->term_structure_size == B->term_structure_size);
    result->term_structure_size = A->term_structure_size;
    float *a_time = A->time_term_structure;
    float *b_time = B->time_term_structure;
    PolynomialExpFunctionSumFloat *result_function = result->piecewise_functions;
    PolynomialExpFunctionSumFloat *a_function = A->piecewise_functions;
    PolynomialExpFunctionSumFloat *b_function = B->piecewise_functions;

    for (uint i = 0; i < result->term_structure_size; ++i) {
        assert(fabsf(*a_time++ - *b_time++) < 1e-9f);
        PolynomialExpFunctionSumFloatMultiply(result_function++, a_function++, b_function++);
    }
}

void PiecewiseFunctionFloatAddConstant(PiecewiseFunctionFloat *result, float constant) {
    PolynomialExpFunctionSumFloat *this_function = result->piecewise_functions;
    for (uint i = 0; i < result->term_structure_size; ++i) {
        PolynomialFloat *this_polynomial = this_function->polynomials;
        for (uint j = 0; j < this_function->number_of_polynomials; ++j) {
            this_polynomial->coefficients[0] += constant;
            ++this_polynomial;
        }
        ++this_function;
    }
}

void PiecewiseFunctionFloatMultiplyByConstant(PiecewiseFunctionFloat *result, float constant) {
    PolynomialExpFunctionSumFloat *this_function = result->piecewise_functions;
    for (uint i = 0; i < result->term_structure_size; ++i) {
        PolynomialFloat *this_polynomial = this_function->polynomials;
        for (uint j = 0; j < this_function->number_of_polynomials; ++j) {
            float *this_coef = this_polynomial->coefficients;
            for (uint k = 0; k < this_polynomial->number_of_coefficients; ++k) {
                *this_coef++ *= constant;
            }
            ++this_polynomial;
        }
        ++this_function;
    }
}

void PiecewiseFunctionFloatDerivative(PiecewiseFunctionFloat *result, PiecewiseFunctionFloat *input) {
    result->term_structure_size = input->term_structure_size;
    memcpy(result->time_term_structure, input->time_term_structure, input->term_structure_size * sizeof(float));
    PolynomialExpFunctionSumFloat *result_function = result->piecewise_functions;
    PolynomialExpFunctionSumFloat *input_function = input->piecewise_functions;
    for (uint i = 0; i < result->term_structure_size; ++i) {
        PolynomialExpFunctionSumFloatDerivative(result_function++, input_function++);
    }
}

void PiecewiseFunctionFloatPrimitive(PiecewiseFunctionFloat *result, PiecewiseFunctionFloat *input) {
    result->term_structure_size = input->term_structure_size;
    memcpy(result->time_term_structure, input->time_term_structure, input->term_structure_size * sizeof(float));
    PolynomialExpFunctionSumFloat *result_function = result->piecewise_functions;
    PolynomialExpFunctionSumFloat *input_function = input->piecewise_functions;
    for (uint i = 0; i < result->term_structure_size; ++i) {
        PolynomialExpFunctionSumFloatPrimitive(result_function++, input_function++);
    }
}

float PiecewiseFunctionFloatIntegral(PiecewiseFunctionFloat *function, float integral_start, float integral_end) {
    if (fabsf(integral_end - integral_start) < 1e-9f) return 0;
    if (integral_start > integral_end) {
        return -PiecewiseFunctionFloatIntegral(function, integral_end, integral_start);
    }
    float function_time_start = function->time_term_structure[0];
    float function_time_end = function->time_term_structure[function->term_structure_size - 1];
    assert(integral_start <= function_time_end);
    assert(integral_start >= function_time_start);
    assert(integral_end <= function_time_end);
    assert(integral_end >= function_time_start);

    float result = 0;
    PolynomialExpFunctionSumFloat *piecewise_function = function->piecewise_functions;
    float *this_time = function->time_term_structure;
    uint term_structure_idx = 0;
    for (uint i = 0; i < function->term_structure_size; ++i) {
        ++this_time;
        if (*this_time > integral_start) {
            break;
        }
        ++piecewise_function;
        ++term_structure_idx;
    }
    PolynomialExpFunctionSumFloat temp_primitive;
    for (uint i = 0; i < POLYNOMIAL_EXP_FUNCTION_SIZE; ++i) {
        temp_primitive.polynomials[i].coefficients = (float *)alloca(POLYNOMIAL_EXP_FUNCTION_NB_COEFS * sizeof(float));
        temp_primitive.polynomials[i].number_of_coefficients = 0;
    }
    PolynomialExpFunctionSumFloatPrimitive(&temp_primitive, piecewise_function);
    result += temp_primitive(*this_time) - temp_primitive(integral_start);
    float last_time;
    for (uint i = term_structure_idx; i < function->term_structure_size; ++i) {
        last_time = *this_time++;
        ++piecewise_function;
        if (*this_time > integral_end - 1e-5f) {
            break;
        }
        PolynomialExpFunctionSumFloatPrimitive(&temp_primitive, piecewise_function);
        result += temp_primitive(*this_time) - temp_primitive(last_time);
    }

    PolynomialExpFunctionSumFloatPrimitive(&temp_primitive, piecewise_function);
    result += temp_primitive(integral_end) - temp_primitive(last_time);
    return result;
}

/**
 * Create a new function that looks like g: t -> exp(Integral_start^t f(s) ds)
 */
void PiecewiseFunctionFloatGetExponentialForwardIntegral(PiecewiseFunctionFloat *result, PiecewiseFunctionFloat *input, float start) {
    result->term_structure_size = input->term_structure_size;
    if (result != input) {
        memcpy(result->time_term_structure, input->time_term_structure, input->term_structure_size * sizeof(float));
    }
    PolynomialExpFunctionSumFloat temp_primitive;

    for (uint i = 0; i < POLYNOMIAL_EXP_FUNCTION_SIZE; ++i) {
        temp_primitive.polynomials[i].coefficients = (float *)alloca(POLYNOMIAL_EXP_FUNCTION_NB_COEFS * sizeof(float));
        memset(temp_primitive.polynomials[i].coefficients, 0, POLYNOMIAL_EXP_FUNCTION_NB_COEFS * sizeof(float));
    }
    float *this_time = input->time_term_structure;
    PolynomialExpFunctionSumFloat *input_function = input->piecewise_functions;
    PolynomialExpFunctionSumFloat *result_function = result->piecewise_functions;
    for (uint i = 0; i < input->term_structure_size - 1; ++i) {
        float integral_to_this_time = expf(PiecewiseFunctionFloatIntegral(input, start, *this_time));
        PolynomialExpFunctionSumFloatPrimitive(&temp_primitive, input_function++);
        float primitive_at_this_time = temp_primitive(*this_time++);
        result_function->number_of_polynomials = 1;
        result_function->exp_coefs[0] = temp_primitive.polynomials[0].coefficients[1];
        result_function->polynomials[0].coefficients[0] = expf(-primitive_at_this_time) * integral_to_this_time;
        result_function->polynomials[0].number_of_coefficients = 1;
        ++result_function;
    }
}

void PiecewiseFunctionFloatGetExponentialBackwardIntegral(PiecewiseFunctionFloat *result, PiecewiseFunctionFloat *input, float end) {
    PiecewiseFunctionFloatGetExponentialForwardIntegral(result, input, end);
    PolynomialExpFunctionSumFloat *result_function = result->piecewise_functions;
    for (uint term_structure_idx = 0; term_structure_idx < input->term_structure_size; ++term_structure_idx) {
        PolynomialFloat *this_polynomial = result_function->polynomials;
        for (uint polynome_idx = 0; polynome_idx < result_function->number_of_polynomials; ++polynome_idx) {
            float *this_coef = this_polynomial->coefficients;
            for (uint coef_idx = 0; coef_idx < this_polynomial->number_of_coefficients; ++coef_idx) {
                *this_coef = -(*this_coef++);
            }
            ++this_polynomial;
        }
        ++result_function;
    }
}

void TestPolynomials() {
    PolynomialFloat polynomial;
    polynomial.number_of_coefficients = 3;
    float coefficients[3] = { 1, 1, 1 };
    polynomial.coefficients = coefficients;
    assert(polynomial(0) == 1);
    assert(polynomial(1) == 3);

    PolynomialFloat derivative;
    derivative.coefficients = (float *)malloc(4 * sizeof(float));
    PolynomialFloatDerivative(&derivative, &polynomial);
    assert(derivative.number_of_coefficients == 2);
    assert(derivative.coefficients[0] == 1);
    assert(derivative.coefficients[1] == 2);
    free(derivative.coefficients);

    PolynomialFloat primitive;
    primitive.coefficients = (float *)malloc(4 * sizeof(float));
    PolynomialFloatPrimitive(&primitive, &polynomial, 1);
    assert(primitive.number_of_coefficients == 4);
    assert(primitive.coefficients[0] == 1.0f);
    assert(primitive.coefficients[1] == 1.0f);
    assert(primitive.coefficients[2] == 0.5f);
    assert(primitive.coefficients[3] == 1.0f / 3.0f);
    free(primitive.coefficients);

}

void TestPolynomialExpFunctionSumFloat() {
    PolynomialExpFunctionSumFloat function;
    PolynomialFloat polynomial_1;
    PolynomialFloat polynomial_2;
    polynomial_1.number_of_coefficients = 3;
    polynomial_2.number_of_coefficients = 3;
    float coefs_1[3] = { 1, 1, 1 };
    float coefs_2[3] = { 2, 2, 1 };
    polynomial_1.coefficients = coefs_1;
    polynomial_2.coefficients = coefs_2;

    function.exp_coefs[0] = -1;
    function.number_of_polynomials = 1;
    function.polynomials[0] = polynomial_1;

    assert(function(0) == 1);
    assert(fabsf(function(1) - 3.0f * expf(-1)) < 1e-6);

    PolynomialExpFunctionSumFloat derivative;
    PolynomialExpFunctionSumFloat primitive;

    float derivative_coefs[POLYNOMIAL_EXP_FUNCTION_NB_COEFS];
    float primitive_coefs[POLYNOMIAL_EXP_FUNCTION_NB_COEFS];
    derivative.polynomials[0].coefficients = derivative_coefs;
    primitive.polynomials[0].coefficients = primitive_coefs;
    PolynomialExpFunctionSumFloatPrimitive(&primitive, &function);
    PolynomialExpFunctionSumFloatDerivative(&derivative, &primitive);

    assert(derivative.exp_coefs[0] == -1);
    assert(derivative.number_of_polynomials == 1);
    assert(derivative.polynomials[0].coefficients[0] == 1);
    assert(derivative.polynomials[0].coefficients[1] == 1);
    assert(derivative.polynomials[0].coefficients[2] == 1);

}

void TestConstantPiecewiseFunction() {
    PiecewiseFunctionFloat function;
    float term_structure[8] = { -1, 0, 1, 2, 3, 4, 5, 6 };
    PolynomialExpFunctionSumFloat piecewise_functions[8];
    function.term_structure_size = 8;
    function.time_term_structure = term_structure;
    function.piecewise_functions = piecewise_functions;
    float constant_array[4] = { 1, 0, 0, 0 };
    for (uint i = 0; i < 8; ++i) {
        piecewise_functions[i].exp_coefs[0] = 0;
        piecewise_functions[i].exp_coefs[1] = 0;
        piecewise_functions[i].exp_coefs[2] = 0;
        piecewise_functions[i].exp_coefs[3] = 0;

        piecewise_functions[i].number_of_polynomials = 1;
        piecewise_functions[i].polynomials[0].coefficients = constant_array;
        piecewise_functions[i].polynomials[0].number_of_coefficients = 1;
    }

    assert(PiecewiseFunctionFloatIntegral(&function, 0, 1) == 1.0f);
    assert(PiecewiseFunctionFloatIntegral(&function, 0, 2.5) == 2.5f);
    assert(PiecewiseFunctionFloatIntegral(&function, -0.5, 1) == 1.5f);
    assert(PiecewiseFunctionFloatIntegral(&function, 0, 4) == 4.0f);

    for (uint i = 0; i < 8; ++i) {
        piecewise_functions[i].polynomials[0].coefficients = (float *)alloca(4 * sizeof(float));
        piecewise_functions[i].polynomials[0].coefficients[0] = (float)term_structure[i];
    }

    assert(PiecewiseFunctionFloatIntegral(&function, 0, 1) == 0.0f);
    assert(PiecewiseFunctionFloatIntegral(&function, 0, 2.5) == 2.0f);

}

void TestPiecewiseFunctionFloatExponential() {

    PiecewiseFunctionFloat function;
    float term_structure[8] = { -1, 0, 1, 2, 3, 4, 5, 6 };
    PolynomialExpFunctionSumFloat piecewise_functions[8];
    function.term_structure_size = 8;
    function.time_term_structure = term_structure;
    function.piecewise_functions = piecewise_functions;
    float constant_array[4] = { 1, 0, 0, 0 };
    for (uint i = 0; i < 8; ++i) {
        piecewise_functions[i].exp_coefs[0] = 0;
        piecewise_functions[i].exp_coefs[1] = 0;
        piecewise_functions[i].exp_coefs[2] = 0;
        piecewise_functions[i].exp_coefs[3] = 0;

        piecewise_functions[i].number_of_polynomials = 1;
        piecewise_functions[i].polynomials[0].coefficients = constant_array;
        piecewise_functions[i].polynomials[0].number_of_coefficients = 1;
    }

    PiecewiseFunctionFloat exp_function;
    exp_function.time_term_structure = (float *)alloca(8 * sizeof(float));
    exp_function.piecewise_functions = (PolynomialExpFunctionSumFloat *)alloca(8 * sizeof(PolynomialExpFunctionSumFloat));
    for (uint i = 0; i < 8; ++i) {
        for (uint j = 0; j < POLYNOMIAL_EXP_FUNCTION_SIZE; ++j) {
            exp_function.piecewise_functions[i].polynomials[j].coefficients = (float *)alloca(POLYNOMIAL_EXP_FUNCTION_NB_COEFS * sizeof(float));
        }
    }
    PiecewiseFunctionFloatGetExponentialForwardIntegral(&exp_function, &function, 0);

    float integral_value = PiecewiseFunctionFloatIntegral(&exp_function, 0, 1);
    float integral_target = expf(1) - 1;
    assert(integral_value == integral_target);
    integral_value = PiecewiseFunctionFloatIntegral(&exp_function, 0, 4.5);
    integral_target = expf(4.5) - 1;
    assert(integral_value == integral_target);

}

/*
int main() {
    TestPolynomials();
    TestPolynomialExpFunctionSumFloat();
    TestConstantPiecewiseFunction();
    TestPiecewiseFunctionFloatExponential();
}
*/
