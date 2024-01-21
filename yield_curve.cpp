#include <math.h>
#include <stdlib.h>

struct YieldCurveFloat{
	float *interpolation_times = nullptr;
	float *short_rates = nullptr;
	size_t number_of_points;
};

float yieldCurveFloatGetInterestRate(YieldCurveFloat *yield_curve, float time_to_maturity) {
	float *this_maturity_time = yield_curve->interpolation_times;
	float *this_short_rate = yield_curve->short_rates;
	for (size_t i = 0; i < yield_curve->number_of_points; ++i) {
		if (*this_maturity_time++ < time_to_maturity) {
			++this_short_rate;
		}
		else {
			break;
		}
	}
	return *this_short_rate;
}

float discountFactorSpotFromYieldCurveFloat(YieldCurveFloat *yield_curve, float time_to_maturity) {
	float *this_maturity_time = yield_curve->interpolation_times;
	float *this_short_rate = yield_curve->short_rates;
	float log_discount_factor = 0;
	for (size_t i = 0; i < yield_curve->number_of_points; ++i) {
		if (*(this_maturity_time + 1) < time_to_maturity) {
			log_discount_factor += *this_short_rate++ * (*(this_maturity_time + 1) - *this_maturity_time++);
		}
		else {
			log_discount_factor += *this_short_rate * (time_to_maturity - *this_maturity_time);
			break;
		}
	}
	return expf(-log_discount_factor);
}

float discountFactorForwardFromYieldCurveFloat(YieldCurveFloat *yield_curve, float time_to_start, float time_to_end) {
	return discountFactorSpotFromYieldCurveFloat(yield_curve, time_to_end) / discountFactorSpotFromYieldCurveFloat(yield_curve, time_to_end);
}

int yieldCurveFloatStripFras(YieldCurveFloat *yield_curve, float *forward_rates, float *expiry_times_in_years, size_t number_of_instruments) {
	float first_expiry_time = 0;
	float second_expiry_time = 0;
	if ((yield_curve->short_rates) && (yield_curve->number_of_points != number_of_instruments)) {
		yield_curve->short_rates = (float *)realloc(yield_curve->short_rates, number_of_instruments * sizeof(float));
	}
	if ((yield_curve->interpolation_times) && (yield_curve->number_of_points != number_of_instruments)) {
		yield_curve->interpolation_times = (float *)realloc(yield_curve->interpolation_times, number_of_instruments * sizeof(float));
	}
	yield_curve->number_of_points = number_of_instruments;
	if (!yield_curve->short_rates) yield_curve->short_rates = (float *)malloc(number_of_instruments * sizeof(float));
	if (!yield_curve->interpolation_times) yield_curve->interpolation_times = (float *)malloc(number_of_instruments * sizeof(float));
	
	float *this_forward_rate = forward_rates;
	float *this_expiry_time = expiry_times_in_years;
	float *yield_curve_short_rate = yield_curve->short_rates;
	float *yield_curve_expiry_time = yield_curve->interpolation_times;
	for (size_t i = 0; i < number_of_instruments; ++i) {
		second_expiry_time = *this_expiry_time;
		*yield_curve_expiry_time++ = *this_expiry_time++;
		float coverage = second_expiry_time - first_expiry_time;
		*yield_curve_short_rate++ = logf(1 + coverage * (*this_forward_rate++)) / coverage;
		first_expiry_time = second_expiry_time;
	}
	return 0;
}

struct YieldCurveSwapF {
	size_t number_of_fixed_payments;
	size_t number_of_floating_payments;
	float swap_rate;
	float *fixed_payment_times = nullptr;
	float *floating_payment_times = nullptr;
};

int yieldCurveSwapFInit(YieldCurveSwapF *swap, size_t number_of_fixed_payments, size_t number_of_floating_payments, float swap_rate) {
	swap->swap_rate = swap_rate;
	if ((swap->fixed_payment_times) && (swap->number_of_fixed_payments != number_of_fixed_payments)) {
		swap->fixed_payment_times = (float *)realloc(swap->fixed_payment_times, number_of_fixed_payments * sizeof(float));
	}
	if (!swap->fixed_payment_times) swap->fixed_payment_times = (float *)malloc(number_of_fixed_payments * sizeof(float));

	if ((swap->floating_payment_times) && (swap->number_of_floating_payments != number_of_floating_payments)) {
		swap->floating_payment_times = (float *)realloc(swap->floating_payment_times, number_of_floating_payments * sizeof(float));
	}
	if (!swap->floating_payment_times) swap->floating_payment_times = (float *)malloc(number_of_floating_payments * sizeof(float));
	swap->number_of_fixed_payments = number_of_fixed_payments;
	swap->number_of_floating_payments = number_of_floating_payments;
	float *fixed_payment_time = swap->fixed_payment_times;
	float *floating_payment_time = swap->floating_payment_times;
	for (size_t i = 0; i < number_of_fixed_payments; ++i) {
		*fixed_payment_time++ = (float)(i + 1);
	}
	for (size_t i = 0; i < number_of_floating_payments; ++i) {
		*floating_payment_time++ = (float)(i + 1);
	}
	return 0;
}


inline void yieldCurveSwapFFree(YieldCurveSwapF *swap) {
	free(swap->fixed_payment_times);
	free(swap->floating_payment_times);
}

int yieldCurveFloatStripSwaps(YieldCurveFloat *yield_curve, YieldCurveSwapF *swap_list, float *swap_rates, size_t number_of_swaps, size_t number_of_already_stripped_fras) {
	float first_expiry_time = 0;
	float second_expiry_time = 0;
	size_t number_of_total_instruments = number_of_swaps + number_of_already_stripped_fras;
	if ((yield_curve->short_rates) && (yield_curve->number_of_points != number_of_total_instruments)) {
		yield_curve->short_rates = (float *)realloc(yield_curve->short_rates, number_of_total_instruments * sizeof(float));
	}
	if ((yield_curve->interpolation_times) && (yield_curve->number_of_points != number_of_total_instruments)) {
		yield_curve->interpolation_times = (float *)realloc(yield_curve->interpolation_times, number_of_total_instruments * sizeof(float));
	}
	yield_curve->number_of_points = number_of_total_instruments;
	if (!yield_curve->short_rates) {
		if (number_of_already_stripped_fras > 0) {
			return -1;
		}
		else {
			yield_curve->short_rates = (float *)malloc(number_of_swaps * sizeof(float));
		}
	}
	if (!yield_curve->interpolation_times) {
		if (number_of_already_stripped_fras > 0) {
			return -1;
		}
		else {
			yield_curve->interpolation_times = (float *)malloc(number_of_swaps * sizeof(float));
		}
	}
	float *this_swap_rate = swap_rates;
	float *previous_swap_rate = swap_rates;
	float *yield_curve_short_rate = yield_curve->short_rates + number_of_already_stripped_fras;
	float *yield_curve_expiry_time = yield_curve->interpolation_times + number_of_already_stripped_fras;
	float fixed_leg_sum = 0;
	size_t number_of_known_swap_discount_factors = 0;
	YieldCurveSwapF *first_swap = swap_list;
	YieldCurveSwapF *last_swap = first_swap + number_of_swaps - 1;
	size_t last_fra_payment_time = yield_curve->interpolation_times[number_of_already_stripped_fras - 1];
	float *this_swap_payment_time = first_swap->fixed_payment_times;
	for (size_t i = 0; i < first_swap->number_of_fixed_payments; ++i) {
		if (last_fra_payment_time > *this_swap_payment_time++) ++number_of_known_swap_discount_factors;
	}
	this_swap_payment_time = first_swap->fixed_payment_times;
	float last_swap_payment_time = 0;
	float this_discount_factor = 1;
	for (size_t i = 0; i < number_of_known_swap_discount_factors; ++i) {
		float coverage = (*this_swap_payment_time - last_swap_payment_time);
		this_discount_factor = discountFactorSpotFromYieldCurveFloat(yield_curve, *this_swap_payment_time);
		fixed_leg_sum += coverage * this_discount_factor;
		last_swap_payment_time = *this_swap_payment_time;
		++this_swap_payment_time;
	}
	float last_discount_factor = this_discount_factor;
	float last_swap_rate = 1;
	YieldCurveSwapF *this_swap = first_swap;
	for (size_t i = 0; i < number_of_swaps; ++i) {
		this_swap_payment_time = this_swap->fixed_payment_times + this_swap->number_of_fixed_payments - 1;
		fixed_leg_sum *= *this_swap_rate / last_swap_rate;
		float coverage = (*this_swap_payment_time - last_swap_payment_time);
		float this_discount_factor = (1 - fixed_leg_sum) / (1 + coverage * (*this_swap_rate));
		*yield_curve_expiry_time++ = *this_swap_payment_time;
		*yield_curve_short_rate++ = logf(last_discount_factor / this_discount_factor) / coverage;
		fixed_leg_sum += coverage * this_discount_factor * (*this_swap_rate);
		last_swap_rate = *this_swap_rate++;
		last_swap_payment_time = *this_swap_payment_time;
		++this_swap;
	}
	return 0;
}

