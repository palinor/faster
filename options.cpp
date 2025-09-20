/**
 * (C) Copyright 2024 Aion Feehan. All Rights Reserved.
 *
 * This software is provided 'as-is', without any express or implied warranty. In no event will the authors
 * be held liable for any damages arising from the use of this software.
 *
 */

#include <math.h>
#ifdef __APPLE__
#include "piecewise_functions.h"
#else
#include "piecewise_functions.cpp"
#endif

enum class OptionType {
	CALL,
	PUT,
};

enum class OptionModelType {
	BLACK_SCHOLES,
	SHIFTED_SABR,
	COUNT
};

enum class OptionVolType {
	NORMAL,
	LOGNORMAL
};

enum class ModelPrecision {
	FLOAT,
	DOUBLE
};


struct OptionModelHeader {
	OptionVolType vol_type;
	OptionModelType model_type;

};

struct OptionModelParametersShiftedSabrFloat {
	OptionModelHeader model_header;
	float sigma_0;
	float alpha;
	float beta;
	float rho;
	float zeta;
};

struct ShiftedSabrParamGrid {
	float *times_to_expiry = nullptr;
	float *underlying_maturity_length = nullptr;
	OptionModelParametersShiftedSabrFloat *parameters = nullptr;
	uint number_of_expiries = 0;
	uint number_of_underlyings = 0;
};




void InterpolateShiftedSabrParametersFloat(void *target, float expiry, float underlying_length, ShiftedSabrParamGrid *param_grid) {
	OptionModelHeader *result_header = (OptionModelHeader *)target;
	OptionModelParametersShiftedSabrFloat *result = (OptionModelParametersShiftedSabrFloat *)target;
	uint expiry_idx_a = 0;
	uint expiry_idx_b = 0;
	uint underlying_idx_a = 0;
	uint underlying_idx_b = 0;
	float *this_expiry = param_grid->times_to_expiry;
	float *this_underlying_length = param_grid->underlying_maturity_length;
	uint expiry_idx = 0;
	uint underlying_idx = 0;
	for (expiry_idx = 0; expiry_idx < param_grid->number_of_expiries; ++expiry_idx) {
		expiry_idx_b = expiry_idx;
		if (*this_expiry > expiry) break;
		++this_expiry;
		expiry_idx_a = expiry_idx_b;
	}
	for (underlying_idx = 0; underlying_idx < param_grid->number_of_underlyings; ++underlying_idx) {
		underlying_idx_b = underlying_idx;
		if (*this_underlying_length > underlying_length) break;
		++this_underlying_length;
		underlying_idx_a = underlying_idx_b;
	}
	OptionModelParametersShiftedSabrFloat *params_expiry_a_underlying_a = param_grid->parameters + expiry_idx_a * param_grid->number_of_underlyings + underlying_idx_a;
	OptionModelParametersShiftedSabrFloat *params_expiry_a_underlying_b = param_grid->parameters + expiry_idx_a * param_grid->number_of_underlyings + underlying_idx_b;
	OptionModelParametersShiftedSabrFloat *params_expiry_b_underlying_a = param_grid->parameters + expiry_idx_b * param_grid->number_of_underlyings + underlying_idx_a;
	OptionModelParametersShiftedSabrFloat *params_expiry_b_underlying_b = param_grid->parameters + expiry_idx_b * param_grid->number_of_underlyings + underlying_idx_b;

	result_header->vol_type = params_expiry_a_underlying_a->model_header.vol_type;
	result_header->model_type = params_expiry_a_underlying_a->model_header.model_type;
	float expiry_time_a = param_grid->times_to_expiry[expiry_idx_a];
	float expiry_time_b = param_grid->times_to_expiry[expiry_idx_b];
	float underlying_length_a = param_grid->underlying_maturity_length[underlying_idx_a];
	float underlying_length_b = param_grid->underlying_maturity_length[underlying_idx_b];
	float underlying_interpolation_coefficient = underlying_idx_a == underlying_idx_b ? 0 : (underlying_length - underlying_length_a) / (underlying_length_b - underlying_length_a);
	float expiry_interpolation_coefficient = expiry_idx_a == expiry_idx_b ? 0 : (expiry - expiry_time_a) / (expiry_time_b - expiry_time_a);

	float sigma_0_expiry_a = params_expiry_a_underlying_a->sigma_0 + underlying_interpolation_coefficient * params_expiry_a_underlying_b->sigma_0;
	float sigma_0_expiry_b = params_expiry_b_underlying_a->sigma_0 + underlying_interpolation_coefficient * params_expiry_b_underlying_b->sigma_0;
	result->sigma_0 = sigma_0_expiry_a + expiry_interpolation_coefficient * sigma_0_expiry_b;

	float alpha_expiry_a = params_expiry_a_underlying_a->alpha + underlying_interpolation_coefficient * params_expiry_a_underlying_b->alpha;
	float alpha_expiry_b = params_expiry_b_underlying_a->alpha + underlying_interpolation_coefficient * params_expiry_b_underlying_b->alpha;
	result->alpha = alpha_expiry_a + expiry_interpolation_coefficient * alpha_expiry_b;

	float beta_expiry_a = params_expiry_a_underlying_a->beta + underlying_interpolation_coefficient * params_expiry_a_underlying_b->beta;
	float beta_expiry_b = params_expiry_b_underlying_a->beta + underlying_interpolation_coefficient * params_expiry_b_underlying_b->beta;
	result->beta = beta_expiry_a + expiry_interpolation_coefficient * beta_expiry_b;

	float rho_expiry_a = params_expiry_a_underlying_a->rho + underlying_interpolation_coefficient * params_expiry_a_underlying_b->rho;
	float rho_expiry_b = params_expiry_b_underlying_a->rho + underlying_interpolation_coefficient * params_expiry_b_underlying_b->rho;
	result->rho = rho_expiry_a + expiry_interpolation_coefficient * rho_expiry_b;

	float zeta_expiry_a = params_expiry_a_underlying_a->zeta + underlying_interpolation_coefficient * params_expiry_a_underlying_b->zeta;
	float zeta_expiry_b = params_expiry_b_underlying_a->zeta + underlying_interpolation_coefficient * params_expiry_b_underlying_b->zeta;
	result->zeta = zeta_expiry_a + expiry_interpolation_coefficient * zeta_expiry_b;
}

float OptionForwardPremiumNormalVol(float normal_vol, float forward, float strike, float time_to_maturity_in_years, OptionType option_type) {
	const float one_over_sqrt_2pi = 0.3989422804;
	float moneyness = 0;
	switch (option_type) {
	case OptionType::CALL: {
		moneyness = forward - strike;
	} break;
	case OptionType::PUT: {
		moneyness = strike - forward;
	} break;
	}
	float d_1 = moneyness / (normal_vol * sqrtf(time_to_maturity_in_years));
	float n_d1 = (1 + erff(d_1 / sqrtf(2))) / 2;
	float result = moneyness * n_d1 + (normal_vol * sqrtf(time_to_maturity_in_years)) * one_over_sqrt_2pi * expf(-d_1 * d_1 / 2);
	return result;
}

float OptionForwardPremiumLognormalVol(float lognormal_vol, float forward, float strike, float time_to_maturity_in_years, OptionType option_type) {
	float moneyness = logf(forward / strike);
	float standard_dev = lognormal_vol * sqrtf(time_to_maturity_in_years);
	float d_1 = moneyness / standard_dev + standard_dev / 2;
	float d_2 = d_1 - standard_dev;
	switch (option_type) {
	case OptionType::CALL: {
		float n_d1 = (1 + erff(d_1 / sqrtf(2))) / 2;
		float n_d2 = (1 + erff(d_2 / sqrtf(2))) / 2;
		return forward * n_d1 - strike * n_d2;
	} break;
	case OptionType::PUT: {
		float n_d1 = (1 + erff(-d_1 / sqrtf(2))) / 2;
		float n_d2 = (1 + erff(-d_2 / sqrtf(2))) / 2;
		return n_d2 * strike - n_d1 * forward;
	} break;
	}
	return 0;
}

float OptionForwardImpliedLognormalVolBrent(float option_forward_premium, float forward, float strike, float time_to_maturity_in_years, OptionType option_type, float tolerance) {
	float left_bracket_vol = 0.000001f;
	float right_bracket_vol = 3.f;
	float option_forward_premium_left_bracket = OptionForwardPremiumLognormalVol(left_bracket_vol, forward, strike, time_to_maturity_in_years, option_type);
	float option_forward_premium_right_bracket = OptionForwardPremiumLognormalVol(right_bracket_vol, forward, strike, time_to_maturity_in_years, option_type);
	if ((option_forward_premium_left_bracket - option_forward_premium) * (option_forward_premium_right_bracket - option_forward_premium) > 0) {
		return -1;
	}
	float closest_vol = 0.0f;
	float closest_premium = 0.0f;
	float left_error = fabsf(option_forward_premium_left_bracket - option_forward_premium);
	float right_error = fabsf(option_forward_premium_right_bracket - option_forward_premium);
	if (left_error < right_error) {
		closest_vol = left_bracket_vol;
		closest_premium = option_forward_premium_left_bracket;
	} else {
		closest_vol = right_bracket_vol;
		closest_premium = option_forward_premium_right_bracket;
	}
	float solution = 0;
	bool we_use_bisection_method = true;
	float last_closest_vol = closest_vol;
	// brent
	while ((left_error > tolerance) && (right_error > tolerance) && (fabs(right_bracket_vol - left_bracket_vol) > tolerance)) {
		if ((option_forward_premium_left_bracket != closest_premium) && (option_forward_premium_right_bracket != closest_premium)) {
			solution = left_bracket_vol * option_forward_premium_right_bracket * closest_premium / ((option_forward_premium_left_bracket - option_forward_premium_right_bracket) * (option_forward_premium_left_bracket - closest_premium))
				+ right_bracket_vol * option_forward_premium_left_bracket * closest_premium / ((option_forward_premium_right_bracket - option_forward_premium_left_bracket) * (option_forward_premium_right_bracket - closest_premium))
				+ closest_vol * option_forward_premium_left_bracket * option_forward_premium_right_bracket / ((closest_premium - option_forward_premium_left_bracket) * (closest_premium - option_forward_premium_right_bracket));
		} else {
			solution = right_bracket_vol - option_forward_premium_right_bracket * (right_bracket_vol - left_bracket_vol) / (option_forward_premium_right_bracket - option_forward_premium_left_bracket);
		}
		bool condition_1 = (solution < (3 * left_bracket_vol + right_bracket_vol) / 4) || (right_bracket_vol < solution);
		bool condition_2 = we_use_bisection_method && (fabsf(solution - right_bracket_vol) > fabsf(right_bracket_vol - closest_vol) / 2);
		bool condition_3 = !we_use_bisection_method && (fabsf(solution - right_bracket_vol) > fabsf(closest_vol - last_closest_vol) / 2);
		bool condition_4 = we_use_bisection_method && (fabsf(right_bracket_vol - closest_vol) < tolerance);
		bool condition_5 = !we_use_bisection_method && (fabsf(closest_vol - last_closest_vol) < tolerance);
		if (condition_1 || condition_2 || condition_3 || condition_4 || condition_5) {
			solution = (left_bracket_vol + right_bracket_vol) / 2;
			we_use_bisection_method = true;
		} else {
			we_use_bisection_method = false;
		}
		float solution_premium = OptionForwardPremiumLognormalVol(solution, forward, strike, time_to_maturity_in_years, option_type);
		last_closest_vol = closest_vol;
		if ((option_forward_premium_left_bracket - option_forward_premium) * (solution_premium - option_forward_premium) < 0) {
			right_bracket_vol = solution;
			option_forward_premium_right_bracket = solution_premium;
		} else {
			left_bracket_vol = solution;
			option_forward_premium_left_bracket = solution_premium;
		}

	}
	return solution;
}


float OptionForwardImpliedNormalVolBrent(float option_forward_premium, float forward, float strike, float time_to_maturity_in_years, OptionType option_type, float tolerance) {
	float left_bracket_vol = 0.0f;
	float right_bracket_vol = 0.3f;
	float option_forward_premium_left_bracket = OptionForwardPremiumNormalVol(left_bracket_vol, forward, strike, time_to_maturity_in_years, option_type);
	float option_forward_premium_right_bracket = OptionForwardPremiumNormalVol(right_bracket_vol, forward, strike, time_to_maturity_in_years, option_type);
	if ((option_forward_premium_left_bracket - option_forward_premium) * (option_forward_premium_right_bracket - option_forward_premium) > 0) {
		return -1;
	}
	float closest_vol = 0.0f;
	float closest_premium = 0.0f;
	float left_error = fabsf(option_forward_premium_left_bracket - option_forward_premium);
	float right_error = fabsf(option_forward_premium_right_bracket - option_forward_premium);
	if (left_error < right_error) {
		closest_vol = left_bracket_vol;
		closest_premium = option_forward_premium_left_bracket;
	} else {
		closest_vol = right_bracket_vol;
		closest_premium = option_forward_premium_right_bracket;
	}
	float solution = 0;
	bool we_use_bisection_method = true;
	float last_closest_vol = closest_vol;
	// brent
	while ((left_error > tolerance) && (right_error > tolerance) && (fabs(right_bracket_vol - left_bracket_vol) > tolerance)) {
		if ((option_forward_premium_left_bracket != closest_premium) && (option_forward_premium_right_bracket != closest_premium)) {
			solution = left_bracket_vol * option_forward_premium_right_bracket * closest_premium / ((option_forward_premium_left_bracket - option_forward_premium_right_bracket) * (option_forward_premium_left_bracket - closest_premium))
				+ right_bracket_vol * option_forward_premium_left_bracket * closest_premium / ((option_forward_premium_right_bracket - option_forward_premium_left_bracket) * (option_forward_premium_right_bracket - closest_premium))
				+ closest_vol * option_forward_premium_left_bracket * option_forward_premium_right_bracket / ((closest_premium - option_forward_premium_left_bracket) * (closest_premium - option_forward_premium_right_bracket));
		} else {
			solution = right_bracket_vol - option_forward_premium_right_bracket * (right_bracket_vol - left_bracket_vol) / (option_forward_premium_right_bracket - option_forward_premium_left_bracket);
		}
		bool condition_1 = (solution < (3 * left_bracket_vol + right_bracket_vol) / 4) || (right_bracket_vol < solution);
		bool condition_2 = we_use_bisection_method && (fabsf(solution - right_bracket_vol) > fabsf(right_bracket_vol - closest_vol) / 2);
		bool condition_3 = !we_use_bisection_method && (fabsf(solution - right_bracket_vol) > fabsf(closest_vol - last_closest_vol) / 2);
		bool condition_4 = we_use_bisection_method && (fabsf(right_bracket_vol - closest_vol) < tolerance);
		bool condition_5 = !we_use_bisection_method && (fabsf(closest_vol - last_closest_vol) < tolerance);
		if (condition_1 || condition_2 || condition_3 || condition_4 || condition_5) {
			solution = (left_bracket_vol + right_bracket_vol) / 2;
			we_use_bisection_method = true;
		} else {
			we_use_bisection_method = false;
		}
		float solution_premium = OptionForwardPremiumNormalVol(solution, forward, strike, time_to_maturity_in_years, option_type);
		last_closest_vol = closest_vol;
		if ((option_forward_premium_left_bracket - option_forward_premium) * (solution_premium - option_forward_premium) < 0) {
			right_bracket_vol = solution;
			option_forward_premium_right_bracket = solution_premium;
		} else {
			left_bracket_vol = solution;
			option_forward_premium_left_bracket = solution_premium;
		}

	}
	return solution;
}

/*
* Implied Normal Volatility, Peter Jackel, 2017
*/
float OptionForwardImpliedNormalVolJackel(float premium, float forward, float strike, float time_to_maturity_in_years, OptionType option_type) {
	float theta = option_type == OptionType::CALL ? -1.0f : 1.0f;
	float moneyness = forward - strike > 0 ? forward - strike : 0;
	float psi_star = -fabsf(premium - (theta * moneyness)) / fabsf(strike - forward);
	const float psi_barrier = -0.001882039271f;
	const float one_over_sqrt_2pi = 0.3989422804f;
	float x_bar = 0.0f;
	if (psi_star < psi_barrier) {
		float g = 1 / (psi_star - 0.5f);
		float g_2 = g * g;
		const float numerator_coef_1 = 0.032114372355f;
		const float numerator_coef_2 = 0.016969777977f;
		const float numerator_coef_3 = 2.6207332461e-3f;
		const float numerator_coef_4 = 9.6066952861e-5f;
		const float denominator_coef_1 = 0.6635646938f;
		const float denominator_coef_2 = 0.14528712196f;
		const float denominator_coef_3 = 0.010472855461f;
		float xi_numerator = numerator_coef_1 - g_2 * (numerator_coef_2 - g_2 * (numerator_coef_3 - g_2 * numerator_coef_4));
		float xi_denominator = 1 - g_2 * (denominator_coef_1 - g_2 * (denominator_coef_2 - denominator_coef_3 * g_2));
		float xi = xi_numerator / xi_denominator;
		x_bar = g * (one_over_sqrt_2pi + xi * g_2);
	} else {
		float h = sqrtf(-logf(-psi_star));
		const float numerator_coef_1 = 9.4883409779f;
		const float numerator_coef_2 = 9.6320903635f;
		const float numerator_coef_3 = 0.58556997323f;
		const float numerator_coef_4 = 2.1464093351f;
		const float denominator_coef_1 = 0.65174820867f;
		const float denominator_coef_2 = 1.5120247828f;
		const float denominator_coef_3 = 6.6437847132e-5f;
		float x_bar_numerator = numerator_coef_1 - h * (numerator_coef_2 - h * (numerator_coef_3 + numerator_coef_4 * h));
		float x_bar_denominator = 1 - h * (denominator_coef_1 + h * (denominator_coef_2 + denominator_coef_3 * h));
		x_bar = x_bar_numerator / x_bar_denominator;
	}
	float p_x_bar = expf(-x_bar * x_bar / 2);
	float psi_x_bar = (1 + erff(x_bar / sqrtf(2))) / 2 + p_x_bar / x_bar;
	float q = (psi_x_bar - psi_star) / p_x_bar;
	float result_numerator = 3 * q * x_bar * x_bar * (2 - q * x_bar * (2 + x_bar * x_bar));
	float result_denominator = 6 + q * x_bar * (-12 + x_bar * (6 * q + x_bar * (-6 + q * x_bar * (3 + x_bar * x_bar))));
	float x_star = x_bar + result_numerator / result_denominator;
	return fabsf(forward - strike) / fabsf(x_star * sqrtf(time_to_maturity_in_years));

}

float ShiftedSabrImpliedNormalVol(float forward, float strike, float time_to_maturity_in_years, float sigma_0, float alpha, float beta, float rho, float zeta) {
	float shifted_forward = forward + zeta;
	float shifted_strike = strike + zeta;
	float moneyness = forward - strike;
	float i_0 = 0;
	float i_1 = 0;
	float z = 0;
	if (beta > 0.99) {
		z = alpha * moneyness / sigma_0;
	} else {
		z = alpha / sigma_0 * (powf(shifted_forward, 1 - beta) - powf(shifted_strike, 1 - beta)) / (1 - beta);
	}
	if (fabsf(moneyness) < 1e-10) {
		i_0 = sigma_0 * powf(shifted_strike, beta);
	} else {
		i_0 = alpha * moneyness / logf(
			(sqrtf(1 - 2 * rho * z + z * z) + z - rho) / (1 - rho)
		);
	}
	i_1 = (((beta - 1) * (beta - 1) - 1) * sigma_0 * sigma_0 / powf((shifted_forward + shifted_strike) / 2, 2 * (1 - beta))
		+ 6 * rho * alpha * sigma_0 * beta / powf((shifted_forward + shifted_strike) / 2, 1 - beta)
		+ 2 * alpha * alpha - 3 * rho * rho * alpha * alpha
		) / 24;
	return i_0 * (1 + i_1 * time_to_maturity_in_years);
}


void ShiftedSabrNormalSmile(float *output_vols, float forward, float *strikes, uint n_strikes, float time_to_maturity_in_years, OptionModelParametersShiftedSabrFloat *shifted_sabr_params) {
	float *this_output_vol = output_vols;
	float *this_strike = strikes;
	for (uint i = 0; i < n_strikes; ++i) {
		*this_output_vol++ = ShiftedSabrImpliedNormalVol(
			forward,
			*this_strike++,
			time_to_maturity_in_years,
			shifted_sabr_params->sigma_0,
			shifted_sabr_params->alpha,
			shifted_sabr_params->beta,
			shifted_sabr_params->rho,
			shifted_sabr_params->zeta
		);
	}
}

float CashLevel(float forward, float fixed_swap_leg_coverage, float cms_maturity_in_years) {
	float y = (1 + fixed_swap_leg_coverage * forward);
	return (1 - powf(y, -cms_maturity_in_years)) / forward;
}

float DCashLevel(float forward, float fixed_swap_leg_coverage, float cms_maturity_in_years) {
	float y = (1 + fixed_swap_leg_coverage * forward);
	return cms_maturity_in_years * fixed_swap_leg_coverage * powf(y, -cms_maturity_in_years - 1) / forward - (1 - powf(y, -cms_maturity_in_years)) / (forward * forward);
}

float D2CashLevel(float forward, float fixed_swap_leg_coverage, float cms_maturity_in_years) {
	float y = (1 + fixed_swap_leg_coverage * forward);
	return -cms_maturity_in_years * (cms_maturity_in_years + 1) * fixed_swap_leg_coverage * fixed_swap_leg_coverage * powf(y, -cms_maturity_in_years - 2) / forward
		- 2 * cms_maturity_in_years * fixed_swap_leg_coverage * powf(y, -cms_maturity_in_years - 1) / (forward * forward)
		+ 2 * (1 - powf(y, -cms_maturity_in_years)) / (forward * forward * forward);
}


float CmsReplicationIntegralFunction(float strike, float fixed_swap_leg_coverage, float cms_maturity_in_years) {
	if (fabsf(strike) < 1e-7) {
		return cms_maturity_in_years * (cms_maturity_in_years + 1) * (cms_maturity_in_years + 2)
			* fixed_swap_leg_coverage * fixed_swap_leg_coverage * fixed_swap_leg_coverage / 3;
	}
	float cash_level = CashLevel(strike, fixed_swap_leg_coverage, cms_maturity_in_years);
	return (
		2 * DCashLevel(strike, fixed_swap_leg_coverage, cms_maturity_in_years) * (strike / cash_level - 1)
		- strike * D2CashLevel(strike, fixed_swap_leg_coverage, cms_maturity_in_years)
		) / (cash_level * cash_level);
}


float CmsReplicationForwardPremium(float forward, float *normal_vols, float *strikes, uint n_strikes, float time_to_payment_in_years, float fixed_swap_leg_coverage, float cms_maturity_in_years) {
	float cash_level = CashLevel(forward, fixed_swap_leg_coverage, cms_maturity_in_years);
	float result = forward;
	uint number_of_puts = 0;
	float *this_strike = strikes;
	if ((*this_strike > forward) || ((*(this_strike + n_strikes - 1)) < forward)) {
		float starting_strike = forward - 3 * (*normal_vols);
		float ending_strike = forward + 3 * (*normal_vols);
		float strike_step = (ending_strike - starting_strike) / ((float)n_strikes);
		for (uint strike_idx = 0; strike_idx < n_strikes; ++strike_idx) {
			*this_strike++ = starting_strike + (float)strike_step * (float)strike_idx;
		}
		this_strike = strikes;
	}
	while (*this_strike++ < forward) {
		++number_of_puts;
	}
	float put_integral = 0;
	float call_integral = 0;
	this_strike = strikes;
	float last_strike = *this_strike;
	float *this_vol = normal_vols;

	float cm_replication_value = CmsReplicationIntegralFunction(
		*this_strike,
		fixed_swap_leg_coverage,
		cms_maturity_in_years
	);
	float last_call_price = OptionForwardPremiumNormalVol(
		*(this_vol + number_of_puts),
		forward,
		*(this_strike + number_of_puts),
		time_to_payment_in_years,
		OptionType::CALL
	);
	float last_call_integration_point = cm_replication_value * last_call_price;
	float last_put_price = OptionForwardPremiumNormalVol(
		*this_vol++,
		forward,
		*this_strike++,
		time_to_payment_in_years,
		OptionType::PUT
	);
	float last_put_integration_point = cm_replication_value * last_put_price;

	for (uint i = 1; i < n_strikes; ++i) {
		float cm_function_value = CmsReplicationIntegralFunction(
			*this_strike,
			fixed_swap_leg_coverage,
			cms_maturity_in_years
		);
		if (i < number_of_puts) {
			float strike_distance = (*this_strike - last_strike);
			last_strike = *this_strike;
			float put_price = OptionForwardPremiumNormalVol(
				*this_vol++,
				forward,
				*this_strike++,
				time_to_payment_in_years,
				OptionType::PUT
			);
			float integration_value = cm_function_value * put_price;
			put_integral += (last_put_integration_point + integration_value) * strike_distance / 2;
			last_put_integration_point = integration_value;
		} else if (i == number_of_puts) {
			++this_strike;
			++this_vol;
		} else {
			float strike_distance = (*this_strike - last_strike);
			last_strike = *this_strike;
			float call_price = OptionForwardPremiumNormalVol(
				*this_vol++,
				forward,
				*this_strike++,
				time_to_payment_in_years,
				OptionType::PUT
			);
			float integration_value = cm_function_value * call_price;
			call_integral += (last_call_integration_point + integration_value) * strike_distance / 2;
			last_call_integration_point = integration_value;
		}
	}

	result += (put_integral + call_integral) / cash_level;
	return result;
}


struct CharacteristicFunctionProjectionConstants {

	float time_to_start;
	float time_to_end;
	float strike;
	float discount_factor;
	float forward;

};

enum class FunctionContextType {
	NONE,
	VANILLA_CALL_PRICING,
	GET_VOL
};

struct FunctionContextHeader {
	FunctionContextType context_type = FunctionContextType::NONE;

};

struct GetVolContext {
	FunctionContextHeader context_header = { FunctionContextType::GET_VOL };
};

struct VanillaCallPricingFunctionContext {
	FunctionContextHeader context_header = { FunctionContextType::VANILLA_CALL_PRICING };
	CharacteristicFunctionProjectionConstants constants;
	float (*get_vol)(float, float, float, void *);
	GetVolContext get_vol_context;
};


// this is the theta(z) function we want to project on to the Chebyshev series
float VanillaCallPricingFunction(
	float x,
	void *context
) {
	FunctionContextHeader *context_header = (FunctionContextHeader *)context;
	if (context_header->context_type != FunctionContextType::VANILLA_CALL_PRICING) {
		return -1;
	}
	VanillaCallPricingFunctionContext *pricing_context = (VanillaCallPricingFunctionContext *)(context);
	CharacteristicFunctionProjectionConstants constants = pricing_context->constants;
	float equivalent_strike = expf(-x) * constants.strike - 1;
	equivalent_strike /= (constants.time_to_end - constants.time_to_start);
	float time_to_expiry = constants.time_to_end;
	float sigma = pricing_context->get_vol(equivalent_strike, constants.time_to_end, constants.forward, &(pricing_context->get_vol_context));
	float result = OptionForwardPremiumNormalVol(sigma, constants.forward, equivalent_strike, constants.time_to_end, OptionType::CALL);
	result /= constants.discount_factor;
	return result;
}


void ProjectVanillaCallPricingFunction(
	float *projection_coefficients,
	CharacteristicFunctionProjectionConstants *constants,
	float (*get_vol)(float, float, float, void *),
	void *get_vol_context,
	size_t chebyshev_degree,
	size_t chebyshev_precision
) {
	PolynomialFloat *chebyshev_polynomials;
	PolynomialFloat projected_result;
	projected_result.coefficients = projection_coefficients;
	projected_result.number_of_coefficients = chebyshev_degree + 1;
	PolynomialFloat32ChebyshevMallocAllocate(chebyshev_polynomials, chebyshev_degree);
	PolynomialFloat32CreateChebyshev(chebyshev_polynomials, chebyshev_degree);

	float *chebyshev_coefficients = (float *)malloc(sizeof(float) * (chebyshev_degree + 1));
	VanillaCallPricingFunctionContext pricing_context;
	pricing_context.constants = *constants;
	pricing_context.get_vol = get_vol;
	pricing_context.get_vol_context = *(GetVolContext*)get_vol_context;
	PolynomialFloat32ProjectChebyshev(chebyshev_coefficients, chebyshev_degree, VanillaCallPricingFunction, &pricing_context, chebyshev_precision);
	PolynomialFloat32ChebyshevProjectionPolynomial(&projected_result, chebyshev_coefficients, chebyshev_polynomials, chebyshev_degree);

}




struct LGM1FTermStructureFloat {

	uint number_of_points;
	float r_0;
	float *tenors;
	PiecewiseFunctionFloat *lambda_term_structure;
	PiecewiseFunctionFloat *sigma_term_structure;
	PiecewiseFunctionFloat *mu_term_structure;
	float **lambdas; // these are views to the coefficients in the term structure elements
	float **mus;
	float **sigmas;

};

void LGM1FSetConstantLambda(LGM1FTermStructureFloat *result, float lambda) {
	PiecewiseFunctionFloat *lambda_term_structure = result->lambda_term_structure;
	PiecewiseFunctionFloat *sigma_term_structure = result->sigma_term_structure;
	PiecewiseFunctionFloat *mu_term_structure = result->mu_term_structure;
	PolynomialExpFunctionSumFloat *this_lambda_function = lambda_term_structure->piecewise_functions;
	PolynomialExpFunctionSumFloat *this_mu_function = mu_term_structure->piecewise_functions;
	PolynomialExpFunctionSumFloat *this_sigma_function = sigma_term_structure->piecewise_functions;
	for (uint i = 0; i < lambda_term_structure->term_structure_size; ++i) {
		this_lambda_function->number_of_polynomials = 1;
		this_mu_function->number_of_polynomials = 1;
		this_sigma_function->number_of_polynomials = 1;
		this_lambda_function->exp_coefs[0] = 0;
		this_sigma_function->exp_coefs[0] = 0;
		this_mu_function->exp_coefs[0] = 0;
		this_lambda_function->polynomials[0].number_of_coefficients = 1;
		this_mu_function->polynomials[0].number_of_coefficients = 1;
		this_sigma_function->polynomials[0].number_of_coefficients = 1;
		this_lambda_function->polynomials[0].coefficients[0] = lambda;
		this_sigma_function->polynomials[0].coefficients[0] = 0;
		this_mu_function->polynomials[0].coefficients[0] = 0;
		++this_lambda_function;	
		++this_mu_function;	
		++this_sigma_function;	
	}
}

void LGM1FSetLambdaTermStructure(LGM1FTermStructureFloat *result, float *lambdas) {
	PiecewiseFunctionFloat *lambda_term_structure = result->lambda_term_structure;
	PolynomialExpFunctionSumFloat *this_function = lambda_term_structure->piecewise_functions;
	float *this_source_lambda = lambdas;
	for (uint i = 0; i < lambda_term_structure->term_structure_size; ++i) {
		this_function->number_of_polynomials = 1;
		this_function->exp_coefs[0] = 0;
		this_function->polynomials[0].coefficients[0] = *this_source_lambda++;
		++this_function;
	}
}


void LGM1FCapitalLambda(PiecewiseFunctionFloat *result, PiecewiseFunctionFloat *lambda_term_structure, float T_end) {
	PiecewiseFunctionFloat temp_backward;
	temp_backward.term_structure_size = lambda_term_structure->term_structure_size;
	temp_backward.piecewise_functions = (PolynomialExpFunctionSumFloat *)alloca(lambda_term_structure->term_structure_size * sizeof(PolynomialExpFunctionSumFloat));
	temp_backward.time_term_structure = (float *)alloca(lambda_term_structure->term_structure_size * sizeof(float));
	memcpy(temp_backward.time_term_structure, lambda_term_structure->time_term_structure, lambda_term_structure->term_structure_size * sizeof(float));
	PolynomialExpFunctionSumFloat *this_source_function = lambda_term_structure->piecewise_functions;
	PolynomialExpFunctionSumFloat *this_destination_function = temp_backward.piecewise_functions;
	for (uint i = 0; i < lambda_term_structure->term_structure_size; ++i) {
		this_destination_function->number_of_polynomials = this_source_function->number_of_polynomials;
		memcpy(this_destination_function->exp_coefs, this_source_function->exp_coefs, this_source_function->number_of_polynomials * sizeof(float));
		PolynomialFloat *this_source_polynomial = this_source_function->polynomials;
		PolynomialFloat *this_destination_polynomial = this_destination_function->polynomials;
		for (uint j = 0; j < this_source_function->number_of_polynomials; ++j) {
			this_destination_polynomial->coefficients = (float *)alloca(POLYNOMIAL_EXP_FUNCTION_NB_COEFS * sizeof(float));
			this_destination_polynomial->number_of_coefficients = this_source_polynomial->number_of_coefficients;
			float *this_source_coef = this_source_polynomial->coefficients;
			float *this_destination_coef = this_destination_polynomial->coefficients;
			for (uint k = 0; k < this_source_polynomial->number_of_coefficients; ++k) {
				*this_destination_coef++ = -1 * (*this_source_coef++);
			}
			++this_destination_polynomial;
			++this_source_polynomial;
		}
		++this_destination_function;
		++this_source_function;
	}
	// at this point temp_backward is f : t -> -lambda (piecewise)
	// we turn it in to g : t -> exp(-int_0^t lambda_v dv)
	PiecewiseFunctionFloatGetExponentialForwardIntegral(result, &temp_backward, 0);
	// and then integrate to get G its primitive
	PiecewiseFunctionFloatPrimitive(result, result);
	// now we use temp_backward to store exp(int_0^t lambda_v dv)
	PiecewiseFunctionFloatGetExponentialForwardIntegral(&temp_backward, lambda_term_structure, 0);
	float epsilon = 1e-6f; // just to make sure we stay on the left side of the interval
	float this_integral_value = (*result)(T_end - epsilon);
	// we evaluate the integral from t to T_end as G(T_end) - G(t)
	PiecewiseFunctionFloatMultiplyByConstant(result, -1);
	PiecewiseFunctionFloatAddConstant(result, this_integral_value);
	// and multiply by the remaining function
	PiecewiseFunctionFloatMultiply(result, result, &temp_backward);
}


void LGM1FDiscountFactorVolProcess(PiecewiseFunctionFloat *result, PiecewiseFunctionFloat *lambda_term_structure, PiecewiseFunctionFloat *sigma_term_structure, float T_end) {
	LGM1FCapitalLambda(result, lambda_term_structure, T_end);
	PiecewiseFunctionFloatMultiply(result, result, sigma_term_structure);
	PiecewiseFunctionFloatMultiplyByConstant(result, -1);
}


float LGM1FDiscountFactorDrift(LGM1FTermStructureFloat *lgm_term_structure, float start_time, float end_time) {
	PiecewiseFunctionFloat drift_piecewise_function;
	PiecewiseFunctionFloat lambda_x_mu_piecewise_function;
	drift_piecewise_function.term_structure_size = lgm_term_structure->lambda_term_structure->term_structure_size;
	lambda_x_mu_piecewise_function.term_structure_size = lgm_term_structure->lambda_term_structure->term_structure_size;
	drift_piecewise_function.piecewise_functions = (PolynomialExpFunctionSumFloat *)alloca(drift_piecewise_function.term_structure_size * sizeof(PolynomialExpFunctionSumFloat));
	lambda_x_mu_piecewise_function.piecewise_functions = (PolynomialExpFunctionSumFloat *)alloca(drift_piecewise_function.term_structure_size * sizeof(PolynomialExpFunctionSumFloat));
	for (uint i = 0; i < drift_piecewise_function.term_structure_size; ++i) {
		for (uint j = 0; j < POLYNOMIAL_EXP_FUNCTION_SIZE; ++j) {
			drift_piecewise_function.piecewise_functions[i].polynomials[j].coefficients = (float *)alloca(POLYNOMIAL_EXP_FUNCTION_NB_COEFS * sizeof(float));
			lambda_x_mu_piecewise_function.piecewise_functions[i].polynomials[j].coefficients = (float *)alloca(POLYNOMIAL_EXP_FUNCTION_NB_COEFS * sizeof(float));
		}
	}

	LGM1FCapitalLambda(&drift_piecewise_function, lgm_term_structure->lambda_term_structure, end_time);
	float constant_term = lgm_term_structure->r_0 * drift_piecewise_function(start_time);
	PiecewiseFunctionFloatMultiply(&lambda_x_mu_piecewise_function, lgm_term_structure->lambda_term_structure, lgm_term_structure->mu_term_structure);
	PiecewiseFunctionFloatMultiply(&drift_piecewise_function, &drift_piecewise_function, &drift_piecewise_function);
	return PiecewiseFunctionFloatIntegral(&drift_piecewise_function, start_time, end_time) + constant_term;
}

float LGM1FDiscountFactorVariance(LGM1FTermStructureFloat *lgm_term_structure, float start_time, float end_time) {
	PiecewiseFunctionFloat variance_piecewise_function;
	variance_piecewise_function.term_structure_size = lgm_term_structure->lambda_term_structure->term_structure_size;
	variance_piecewise_function.term_structure_size = lgm_term_structure->lambda_term_structure->term_structure_size;
	variance_piecewise_function.piecewise_functions = (PolynomialExpFunctionSumFloat *)alloca(variance_piecewise_function.term_structure_size * sizeof(PolynomialExpFunctionSumFloat));
	for (uint i = 0; i < variance_piecewise_function.term_structure_size; ++i) {
		for (uint j = 0; j < POLYNOMIAL_EXP_FUNCTION_SIZE; ++j) {
			variance_piecewise_function.piecewise_functions[i].polynomials[j].coefficients = (float *)alloca(POLYNOMIAL_EXP_FUNCTION_NB_COEFS * sizeof(float));
		}
	}

	LGM1FCapitalLambda(&variance_piecewise_function, lgm_term_structure->lambda_term_structure, end_time);
	PiecewiseFunctionFloatMultiply(&variance_piecewise_function, &variance_piecewise_function, lgm_term_structure->sigma_term_structure);
	PiecewiseFunctionFloatMultiply(&variance_piecewise_function, &variance_piecewise_function, &variance_piecewise_function);
	return PiecewiseFunctionFloatIntegral(&variance_piecewise_function, start_time, end_time);
}


float LGM1FDiscountFactorTForwardChangeOfMeasure(LGM1FTermStructureFloat *lgm_term_structure, float start_time, float end_time, float t_forward) {
	PiecewiseFunctionFloat t_forward_drift_function;
	PiecewiseFunctionFloat capital_lambda_t_forward;
	t_forward_drift_function.term_structure_size = lgm_term_structure->lambda_term_structure->term_structure_size;
	capital_lambda_t_forward.term_structure_size = lgm_term_structure->lambda_term_structure->term_structure_size;
	t_forward_drift_function.term_structure_size = lgm_term_structure->lambda_term_structure->term_structure_size;
	capital_lambda_t_forward.term_structure_size = lgm_term_structure->lambda_term_structure->term_structure_size;
	t_forward_drift_function.piecewise_functions = (PolynomialExpFunctionSumFloat *)alloca(t_forward_drift_function.term_structure_size * sizeof(PolynomialExpFunctionSumFloat));
	capital_lambda_t_forward.piecewise_functions = (PolynomialExpFunctionSumFloat *)alloca(t_forward_drift_function.term_structure_size * sizeof(PolynomialExpFunctionSumFloat));
	for (uint i = 0; i < t_forward_drift_function.term_structure_size; ++i) {
		for (uint j = 0; j < POLYNOMIAL_EXP_FUNCTION_SIZE; ++j) {
			t_forward_drift_function.piecewise_functions[i].polynomials[j].coefficients = (float *)alloca(POLYNOMIAL_EXP_FUNCTION_NB_COEFS * sizeof(float));
			capital_lambda_t_forward.piecewise_functions[i].polynomials[j].coefficients = (float *)alloca(POLYNOMIAL_EXP_FUNCTION_NB_COEFS * sizeof(float));
		}
	}

	LGM1FCapitalLambda(&t_forward_drift_function, lgm_term_structure->lambda_term_structure, end_time);
	LGM1FCapitalLambda(&capital_lambda_t_forward, lgm_term_structure->lambda_term_structure, t_forward);
	PiecewiseFunctionFloatMultiply(&t_forward_drift_function, &t_forward_drift_function, lgm_term_structure->sigma_term_structure);
	PiecewiseFunctionFloatMultiply(&t_forward_drift_function, &t_forward_drift_function, &t_forward_drift_function);
	PiecewiseFunctionFloatMultiply(&t_forward_drift_function, &t_forward_drift_function, &capital_lambda_t_forward);
	return - PiecewiseFunctionFloatIntegral(&t_forward_drift_function, start_time, end_time);
}

float LGM1FDiscountFactorValue(LGM1FTermStructureFloat *lgm_term_structure, float start_time, float end_time) {
	float log_discount_factor_drift = LGM1FDiscountFactorDrift(lgm_term_structure, start_time, end_time);
	float log_discount_factor_variance = LGM1FDiscountFactorVariance(lgm_term_structure, start_time, end_time);
	return expf(-log_discount_factor_drift + log_discount_factor_variance / 2.0f);
}


float LGM1FGetTermCapletLogVariance(LGM1FTermStructureFloat *lgm_term_structure, float time_to_expiry, float underlying_length) {
	PiecewiseFunctionFloat vol_process_t1, vol_process_t2;
	vol_process_t1.term_structure_size = lgm_term_structure->lambda_term_structure->term_structure_size;
	vol_process_t2.term_structure_size = lgm_term_structure->lambda_term_structure->term_structure_size;
	vol_process_t1.term_structure_size = lgm_term_structure->lambda_term_structure->term_structure_size;
	vol_process_t2.term_structure_size = lgm_term_structure->lambda_term_structure->term_structure_size;
	vol_process_t1.piecewise_functions = (PolynomialExpFunctionSumFloat *)alloca(vol_process_t1.term_structure_size * sizeof(PolynomialExpFunctionSumFloat));
	vol_process_t2.piecewise_functions = (PolynomialExpFunctionSumFloat *)alloca(vol_process_t2.term_structure_size * sizeof(PolynomialExpFunctionSumFloat));
	vol_process_t1.time_term_structure = (float *)alloca(vol_process_t1.term_structure_size * sizeof(float));
	vol_process_t2.time_term_structure = (float *)alloca(vol_process_t2.term_structure_size * sizeof(float));
	for (uint i = 0; i < vol_process_t1.term_structure_size; ++i) {
		vol_process_t1.piecewise_functions[i].number_of_polynomials = 0;
		vol_process_t2.piecewise_functions[i].number_of_polynomials = 0;
		for (uint j = 0; j < POLYNOMIAL_EXP_FUNCTION_SIZE; ++j) {
			vol_process_t1.piecewise_functions[i].polynomials[j].coefficients = (float *)alloca(POLYNOMIAL_EXP_FUNCTION_NB_COEFS * sizeof(float));
			vol_process_t1.piecewise_functions[i].polynomials[j].number_of_coefficients = 0;
			vol_process_t1.piecewise_functions[i].exp_coefs[j] = 0;
			vol_process_t2.piecewise_functions[i].polynomials[j].coefficients = (float *)alloca(POLYNOMIAL_EXP_FUNCTION_NB_COEFS * sizeof(float));
			vol_process_t2.piecewise_functions[i].polynomials[j].number_of_coefficients = 0;
			vol_process_t2.piecewise_functions[i].exp_coefs[j] = 0;
		}
	}

	LGM1FDiscountFactorVolProcess(&vol_process_t1, lgm_term_structure->lambda_term_structure, lgm_term_structure->sigma_term_structure, time_to_expiry);
	LGM1FDiscountFactorVolProcess(&vol_process_t2, lgm_term_structure->lambda_term_structure, lgm_term_structure->sigma_term_structure, time_to_expiry + underlying_length);

	PiecewiseFunctionFloatMultiplyByConstant(&vol_process_t1, -1);
	PiecewiseFunctionFloatAdd(&vol_process_t1, &vol_process_t1, &vol_process_t2); // at this point vol_process_t1 is the vol process of the shifted libor

	PiecewiseFunctionFloatMultiply(&vol_process_t1, &vol_process_t1, &vol_process_t1);
	return PiecewiseFunctionFloatIntegral(&vol_process_t1, 0, time_to_expiry);

}


float LGM1FGetCompoundedCapletLognormalVol(LGM1FTermStructureFloat *lgm_term_structure, float time_to_expiry, float underlying_length, uint number_of_compoundings) {
	float daily_coverage = underlying_length / (float)number_of_compoundings;
	PiecewiseFunctionFloat vol_process_ti, vol_process_tend;
	vol_process_ti.term_structure_size = lgm_term_structure->lambda_term_structure->term_structure_size;
	vol_process_tend.term_structure_size = lgm_term_structure->lambda_term_structure->term_structure_size;
	vol_process_ti.term_structure_size = lgm_term_structure->lambda_term_structure->term_structure_size;
	vol_process_tend.term_structure_size = lgm_term_structure->lambda_term_structure->term_structure_size;
	vol_process_ti.piecewise_functions = (PolynomialExpFunctionSumFloat *)alloca(vol_process_ti.term_structure_size * sizeof(PolynomialExpFunctionSumFloat));
	vol_process_tend.piecewise_functions = (PolynomialExpFunctionSumFloat *)alloca(vol_process_tend.term_structure_size * sizeof(PolynomialExpFunctionSumFloat));
	for (uint i = 0; i < vol_process_ti.term_structure_size; ++i) {
		for (uint j = 0; j < POLYNOMIAL_EXP_FUNCTION_SIZE; ++j) {
			vol_process_ti.piecewise_functions[i].polynomials[j].coefficients = (float *)alloca(POLYNOMIAL_EXP_FUNCTION_NB_COEFS * sizeof(float));
			vol_process_tend.piecewise_functions[i].polynomials[j].coefficients = (float *)alloca(POLYNOMIAL_EXP_FUNCTION_NB_COEFS * sizeof(float));
		}
	}

	LGM1FDiscountFactorVolProcess(
		&vol_process_tend, 
		lgm_term_structure->lambda_term_structure, 
		lgm_term_structure->sigma_term_structure, 
		time_to_expiry + underlying_length
		);
	float integral_start = 0;
	float integral_end = time_to_expiry;
	float result = 0;
	for (uint compounding_idx = 0; compounding_idx < number_of_compoundings; ++compounding_idx) {
		LGM1FDiscountFactorVolProcess(
			&vol_process_ti, 
			lgm_term_structure->lambda_term_structure, 
			lgm_term_structure->sigma_term_structure,
			time_to_expiry + compounding_idx * daily_coverage
			);
		PiecewiseFunctionFloatMultiplyByConstant(&vol_process_tend, -1);
		PiecewiseFunctionFloatAdd(&vol_process_ti, &vol_process_ti, &vol_process_tend);
		PiecewiseFunctionFloatMultiply(&vol_process_ti, &vol_process_ti, &vol_process_ti);
		result += PiecewiseFunctionFloatIntegral(&vol_process_ti, integral_start, integral_end);
		integral_start = integral_end;
		integral_end += daily_coverage;
	}
	return result;
}

void LGM1FSetMuTermStructureFromDiscountFactors(LGM1FTermStructureFloat *lgm_term_structure, float *log_discount_factors, float *d_log_discount_factors, float *d2_log_discount_factors) {
	float *this_time = lgm_term_structure->tenors;
	for (uint term_structure_idx = 0; term_structure_idx < lgm_term_structure->number_of_points; ++term_structure_idx) {
		float *this_mu = lgm_term_structure->mu_term_structure->piecewise_functions[term_structure_idx].polynomials[0].coefficients;

		++this_time;
	}
}

struct CapletVolParams {
	float normal_vol;
	float lognormal_vol;
	float coverage;
	float time_to_expiry;
	float forward;
};

void LGM1FFitSigmaTermStructureFromTermCaplets(
	LGM1FTermStructureFloat *lgm_term_structure, 
	CapletVolParams *vol_params,
	float bisection_precision,
	uint max_number_of_iterations
	) {
	
	CapletVolParams *this_vol_params = vol_params;
	const float sigma_upper_bound = 1;
	const float sigma_lower_bound = 0;
	uint number_of_vols = lgm_term_structure->number_of_points - 2;
	for (uint vol_idx = 0; vol_idx < number_of_vols; ++vol_idx) {
		float *sigma_to_fit = lgm_term_structure->sigma_term_structure->piecewise_functions[vol_idx].polynomials->coefficients;
		float target_vol = this_vol_params->normal_vol;
		float time_to_expiry = this_vol_params->time_to_expiry;
		float target_variance = target_vol * target_vol * time_to_expiry;
		float coverage = this_vol_params->coverage;
		float forward = this_vol_params->forward;
		uint number_of_iterations = 0;
		float error = INFINITY;
		float upper_sigma = sigma_upper_bound;
		float lower_sigma = sigma_lower_bound;
		while ((number_of_iterations < max_number_of_iterations) && (error > bisection_precision)) {
			float middle_sigma = (upper_sigma - lower_sigma) / 2;
			*sigma_to_fit = middle_sigma;
			float middle_log_variance = LGM1FGetTermCapletLogVariance(lgm_term_structure, time_to_expiry, this_vol_params->coverage);
			float middle_variance = expf(middle_log_variance) - 1;
			middle_variance *= (1 + coverage * forward) * (1 + coverage * forward) / (coverage * coverage);
			error = fabsf(middle_variance - target_variance);
			if (middle_variance < target_variance) {
				lower_sigma = middle_sigma;
			}
			else {
				upper_sigma = middle_sigma;
			}
			++number_of_iterations;
		}
		++this_vol_params;
	}
}

float NormalVolToLognormalVolFloat(float normal_vol, float forward, float time_to_expiry) {
	float call_price = OptionForwardPremiumNormalVol(normal_vol, forward, forward, time_to_expiry, OptionType::CALL);
	float implied_lognormal_vol = OptionForwardImpliedLognormalVolBrent(call_price, forward, forward, time_to_expiry, OptionType::CALL, 2e-7); 
	return implied_lognormal_vol;
}

/**
 * Set up to apply the Jamshidian trick to get the swaption vol
 */
float LGM1FGetSwaptionLognormalVol(LGM1FTermStructureFloat *lgm_term_structure, float time_to_expiry, float underlying_length, uint fixed_pay_frequency) {
    return 0;
}
