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

float optionForwardPremiumNormalVol(float normal_vol, float forward, float strike, float time_to_maturity_in_years, OptionType option_type) {
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

float optionForwardPremiumLognormalVol(float lognormal_vol, float forward, float strike, float time_to_maturity_in_years, OptionType option_type) {
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


float optionForwardImpliedNormalVol(float option_forward_premium, float forward, float strike, float time_to_maturity_in_years, OptionType option_type, float tolerance) {
	float left_bracket_vol = 0;
	float right_bracket_vol = 0.3;
	float option_forward_premium_left_bracket = optionForwardPremiumNormalVol(left_bracket_vol, forward, strike, time_to_maturity_in_years, option_type);
	float option_forward_premium_right_bracket = optionForwardPremiumNormalVol(right_bracket_vol, forward, strike, time_to_maturity_in_years, option_type);
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
		float solution_premium = optionForwardPremiumNormalVol(solution, forward, strike, time_to_maturity_in_years, option_type);
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

float shiftedSabrImpliedNormalVol(float forward, float strike, float time_to_maturity_in_years, float sigma_0, float alpha, float beta, float rho, float zeta) {
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


void shiftedSabrNormalSmile(float *output_vols, float forward, float *strikes, size_t n_strikes, float time_to_maturity_in_years, OptionModelParametersShiftedSabrFloat *shifted_sabr_params) {
	float *this_output_vol = output_vols;
	float *this_strike = strikes;
	for (size_t i = 0; i < n_strikes; ++i) {
		*this_output_vol++ = shiftedSabrImpliedNormalVol(
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

float cashLevel(float forward, float fixed_swap_leg_coverage, float cms_maturity_in_years) {
	float y = (1 + fixed_swap_leg_coverage * forward);
	return (1 - powf(y, -cms_maturity_in_years)) / forward;
}

float dCashLevel(float forward, float fixed_swap_leg_coverage, float cms_maturity_in_years) {
	float y = (1 + fixed_swap_leg_coverage * forward);
	return cms_maturity_in_years * fixed_swap_leg_coverage * powf(y, -cms_maturity_in_years - 1) / forward - (1 - powf(y, -cms_maturity_in_years)) / (forward * forward);
}

float d2CashLevel(float forward, float fixed_swap_leg_coverage, float cms_maturity_in_years) {
	float y = (1 + fixed_swap_leg_coverage * forward);
	return -cms_maturity_in_years * (cms_maturity_in_years + 1) * fixed_swap_leg_coverage * fixed_swap_leg_coverage * powf(y, -cms_maturity_in_years - 2) / forward
		- 2 * cms_maturity_in_years * fixed_swap_leg_coverage * powf(y, -cms_maturity_in_years - 1) / (forward * forward)
		+ 2 * (1 - powf(y, -cms_maturity_in_years)) / (forward * forward * forward);
}


float cmsReplicationIntegralFunction(float strike, float fixed_swap_leg_coverage, float cms_maturity_in_years) {
	float cash_level = cashLevel(strike, fixed_swap_leg_coverage, cms_maturity_in_years);
	return (
		2 * dCashLevel(strike, fixed_swap_leg_coverage, cms_maturity_in_years) * (strike / cash_level - 1)
		- strike * d2CashLevel(strike, fixed_swap_leg_coverage, cms_maturity_in_years)
		) / (cash_level * cash_level);
}


float cmsReplicationForwardPremium(float forward, float *normal_vols, float *strikes, size_t n_strikes, float time_to_payment_in_years, float fixed_swap_leg_coverage, float cms_maturity_in_years) {
	float cash_level = cashLevel(forward, fixed_swap_leg_coverage, cms_maturity_in_years);
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

	float last_call_price = optionForwardPremiumNormalVol(
		*(this_vol + number_of_puts),
		forward,
		*(this_strike + number_of_puts),
		time_to_payment_in_years,
		OptionType::CALL
	);
	float last_put_price = optionForwardPremiumNormalVol(
		*this_vol++,
		forward,
		*this_strike++,
		time_to_payment_in_years,
		OptionType::PUT
		);

	for (size_t i = 1; i < n_strikes; ++i) {
		if (i < number_of_puts) {
			float strike_distance = (*this_strike - last_strike);
			float put_price = optionForwardPremiumNormalVol(
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
			float call_price = optionForwardPremiumNormalVol(
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



