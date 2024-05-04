#include <math.h>


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

struct OptionModelParametersShiftedSabrFloat{
	OptionModelHeader *model_header;
	float sigma_0;
	float alpha;
	float beta;
	float rho;
	float zeta;
};

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


float OptionForwardImpliedNormalVolBrent(float option_forward_premium, float forward, float strike, float time_to_maturity_in_years, OptionType option_type, float tolerance) {
	float left_bracket_vol = 0;
	float right_bracket_vol = 0.3;
	float option_forward_premium_left_bracket = OptionForwardPremiumNormalVol(left_bracket_vol, forward, strike, time_to_maturity_in_years, option_type);
	float option_forward_premium_right_bracket = OptionForwardPremiumNormalVol(right_bracket_vol, forward, strike, time_to_maturity_in_years, option_type);
	if (option_forward_premium_left_bracket * option_forward_premium_right_bracket > 0) {
		return -1;
	}
	float closest_vol = 0;
	float closest_premium = 0;
	float left_error = fabsf(option_forward_premium_left_bracket - option_forward_premium);
	float right_error = fabsf(option_forward_premium_right_bracket - option_forward_premium);
	if (left_error < right_error) {
		closest_vol = left_bracket_vol;
		closest_premium = option_forward_premium_left_bracket;
	}
	else {
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
		}
		else {
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
		}
		else {
			we_use_bisection_method = false;
		}
		float solution_premium = OptionForwardPremiumNormalVol(solution, forward, strike, time_to_maturity_in_years, option_type);
		last_closest_vol = closest_vol;
		if (option_forward_premium_left_bracket * solution_premium < 0) {
			right_bracket_vol = solution;
			option_forward_premium_right_bracket = solution_premium;
		}
		else {
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
	// We need the option to be out of the money
	// So if we are in the money, we convert to other option type by call put parity
	switch (option_type) {
	case OptionType::CALL: {
		if (forward < strike) {
			premium += strike - forward;
			option_type = OptionType::PUT;
		}
	} break;
	case OptionType::PUT: {
		if (strike < forward) {
			premium += forward - strike;
			option_type = OptionType::CALL;
		}
	} break;
	default: break;
	}
	float theta = option_type == OptionType::CALL ? -1 : 1;
	float psi_star = theta * premium / (strike - forward);
	const float psi_barrier = -0.001882039271;
	const float one_over_sqrt_2pi = 0.3989422804;
	float x_bar = 0;
	if (psi_star < psi_barrier) {
		float g = 1 / (psi_star - 0.5);
		float g_2 = g * g;
		const float numerator_coef_1 = 0.032114372355;
		const float numerator_coef_2 = 0.016969777977;
		const float numerator_coef_3 = 0.002207332461;
		const float numerator_coef_4 = 9.6066952861e-5;
		const float denominator_coef_1 = 0.6635646938;
		const float denominator_coef_2 = 0.14528712196;
		const float denominator_coef_3 = 0.010472855461;
		float xi_numerator = numerator_coef_1 - g_2 * (numerator_coef_2 - g_2 * (numerator_coef_3 - g_2 * numerator_coef_4));
		float xi_denominator = 1 - g_2 * (denominator_coef_1 - g_2 * (denominator_coef_2 - denominator_coef_3 * g_2));
		float xi = xi_numerator / xi_denominator;
		x_bar = g * (one_over_sqrt_2pi + xi * g_2);
	}
	else {
		float h = sqrtf(-logf(-psi_star));
		const float numerator_coef_1 = 9.4883409779;
		const float numerator_coef_2 = 9.6320903635;
		const float numerator_coef_3 = 0.58556997323;
		const float numerator_coef_4 = 2.1464093351;
		const float denominator_coef_1 = 0.65174820867;
		const float denominator_coef_2 = 1.5120247828;
		const float denominator_coef_3 = 6.6437847132e-5;
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
	}
	else {
		z = alpha / sigma_0 * (powf(shifted_forward, 1 - beta) - powf(shifted_strike, 1 - beta)) / (1 - beta);
	}
	if (fabsf(moneyness) < 1e-10) {
		i_0 = sigma_0 * powf(shifted_strike, beta);
	}
	else {
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


void ShiftedSabrNormalSmile(float *output_vols, float forward, float *strikes, size_t n_strikes, float time_to_maturity_in_years, OptionModelParametersShiftedSabrFloat *shifted_sabr_params) {
	float *this_output_vol = output_vols;
	float *this_strike = strikes;
	for (size_t i = 0; i < n_strikes; ++i) {
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
	float cash_level = CashLevel(strike, fixed_swap_leg_coverage, cms_maturity_in_years);
	return (
		2 * DCashLevel(strike, fixed_swap_leg_coverage, cms_maturity_in_years) * (strike / cash_level - 1)
		- strike * D2CashLevel(strike, fixed_swap_leg_coverage, cms_maturity_in_years)
		) / (cash_level * cash_level);
}


float CmsReplicationForwardPremium(float forward, float *normal_vols, float *strikes, size_t n_strikes, float time_to_payment_in_years, float fixed_swap_leg_coverage, float cms_maturity_in_years) {
	float cash_level = CashLevel(forward, fixed_swap_leg_coverage, cms_maturity_in_years);
	float result = forward;
	size_t number_of_puts = 0;
	float *this_strike = strikes;
	while (*this_strike++ < forward) {
		++number_of_puts;
	}
	float put_integral = 0;
	float call_integral = 0;
	this_strike = strikes;
	float last_strike = *this_strike;
	float *this_vol = normal_vols;

	float last_call_price = OptionForwardPremiumNormalVol(
		*(this_vol + number_of_puts),
		forward,
		*(this_strike + number_of_puts),
		time_to_payment_in_years,
		OptionType::CALL
	);
	float last_put_price = OptionForwardPremiumNormalVol(
		*this_vol++,
		forward,
		*this_strike++,
		time_to_payment_in_years,
		OptionType::PUT
		);

	for (size_t i = 1; i < n_strikes; ++i) {
		if (i < number_of_puts) {
			float strike_distance = (*this_strike - last_strike);
			float put_price = OptionForwardPremiumNormalVol(
				*this_vol++,
				forward,
				*this_strike++,
				time_to_payment_in_years,
				OptionType::PUT
			);
			last_strike = *this_strike;
			put_integral += (last_put_price + put_price) * strike_distance / 2;
			last_put_price = put_price;
		}
		else if (i == number_of_puts) {
			last_strike = *(this_strike + 1);
			++this_strike;
			++this_vol;
		}
		else {
			float strike_distance = (*this_strike - last_strike);
			float call_price = OptionForwardPremiumNormalVol(
				*this_vol++,
				forward,
				*this_strike++,
				time_to_payment_in_years,
				OptionType::PUT
			);
			last_strike = *this_strike;
			call_integral += (last_call_price + call_price) * strike_distance / 2;
			last_call_price = call_price;
		}
	}

	result += (put_integral + call_integral) / cash_level;
	return result;
}



