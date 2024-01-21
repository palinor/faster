//
//  yield_curve.h
//  example_apple_metal
//
//  Created by Aion Feehan on 1/18/24.
//  Copyright © 2024 Warren Moore. All rights reserved.
//

#ifndef yield_curve_h
struct YieldCurveFloat{
    float *interpolation_times = nullptr;
    float *short_rates = nullptr;
    size_t number_of_points;
};

float yieldCurveFloatGetInterestRate(YieldCurveFloat *yield_curve, float time_to_maturity);
float discountFactorSpotFromYieldCurveFloat(YieldCurveFloat *yield_curve, float time_to_maturity);
float discountFactorForwardFromYieldCurveFloat(YieldCurveFloat *yield_curve, float time_to_start, float time_to_end);
int yieldCurveFloatStripFras(YieldCurveFloat *yield_curve, float *forward_rates, float *expiry_times_in_years, size_t number_of_instruments);
struct YieldCurveSwapF {
    size_t number_of_fixed_payments;
    size_t number_of_floating_payments;
    float swap_rate;
    float *fixed_payment_times = nullptr;
    float *floating_payment_times = nullptr;
};

int yieldCurveSwapFInit(YieldCurveSwapF *swap, size_t number_of_fixed_payments, size_t number_of_floating_payments, float swap_rate);

inline void yieldCurveSwapFFree(YieldCurveSwapF *swap);

int yieldCurveFloatStripSwaps(YieldCurveFloat *yield_curve, YieldCurveSwapF *swap_list, float *swap_rates, size_t number_of_swaps, size_t number_of_already_stripped_fras);

#define yield_curve_h


#endif /* yield_curve_h */
