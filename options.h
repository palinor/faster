/**
 * (C) Copyright 2024 Aion Feehan. All Rights Reserved.
 *
 * This software is provided 'as-is', without any express or implied warranty. In no event will the authors
 * be held liable for any damages arising from the use of this software.
 *
 */

#include <math.h>
#include "piecewise_functions.h"

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




void InterpolateShiftedSabrParametersFloat(void *target, float expiry, float underlying_length, ShiftedSabrParamGrid *param_grid);
float OptionForwardPremiumNormalVol(float normal_vol, float forward, float strike, float time_to_maturity_in_years, OptionType option_type);
float OptionForwardPremiumLognormalVol(float lognormal_vol, float forward, float strike, float time_to_maturity_in_years, OptionType option_type);
float OptionForwardImpliedLognormalVolBrent(float option_forward_premium, float forward, float strike, float time_to_maturity_in_years, OptionType option_type, float tolerance);
float OptionForwardImpliedNormalVolBrent(float option_forward_premium, float forward, float strike, float time_to_maturity_in_years, OptionType option_type, float tolerance);
/*
* Implied Normal Volatility, Peter Jackel, 2017
*/
float OptionForwardImpliedNormalVolJackel(float premium, float forward, float strike, float time_to_maturity_in_years, OptionType option_type);
float ShiftedSabrImpliedNormalVol(float forward, float strike, float time_to_maturity_in_years, float sigma_0, float alpha, float beta, float rho, float zeta);
void ShiftedSabrNormalSmile(float *output_vols, float forward, float *strikes, uint n_strikes, float time_to_maturity_in_years, OptionModelParametersShiftedSabrFloat *shifted_sabr_params);
float CashLevel(float forward, float fixed_swap_leg_coverage, float cms_maturity_in_years);
float DCashLevel(float forward, float fixed_swap_leg_coverage, float cms_maturity_in_years);
float D2CashLevel(float forward, float fixed_swap_leg_coverage, float cms_maturity_in_years);
float CmsReplicationIntegralFunction(float strike, float fixed_swap_leg_coverage, float cms_maturity_in_years);
float CmsReplicationForwardPremium(float forward, float *normal_vols, float *strikes, uint n_strikes, float time_to_payment_in_years, float fixed_swap_leg_coverage, float cms_maturity_in_years);

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

void LGM1FSetConstantLambda(LGM1FTermStructureFloat *result, float lambda);
void LGM1FSetLambdaTermStructure(LGM1FTermStructureFloat *result, float *lambdas);
void LGM1FCapitalLambda(PiecewiseFunctionFloat *result, PiecewiseFunctionFloat *lambda_term_structure, float T_end);
void LGM1FDiscountFactorVolProcess(PiecewiseFunctionFloat *result, PiecewiseFunctionFloat *lambda_term_structure, PiecewiseFunctionFloat *sigma_term_structure, float T_end);
float LGM1FDiscountFactorDrift(LGM1FTermStructureFloat *lgm_term_structure, float start_time, float end_time);
float LGM1FDiscountFactorVariance(LGM1FTermStructureFloat *lgm_term_structure, float start_time, float end_time);
float LGM1FDiscountFactorTForwardChangeOfMeasure(LGM1FTermStructureFloat *lgm_term_structure, float start_time, float end_time, float t_forward);
float LGM1FDiscountFactorValue(LGM1FTermStructureFloat *lgm_term_structure, float start_time, float end_time);
float LGM1FGetTermCapletLogVariance(LGM1FTermStructureFloat *lgm_term_structure, float time_to_expiry, float underlying_length);
float LGM1FGetCompoundedCapletLognormalVol(LGM1FTermStructureFloat *lgm_term_structure, float time_to_expiry, float underlying_length, uint number_of_compoundings);
void LGM1FSetMuTermStructureFromDiscountFactors(LGM1FTermStructureFloat *lgm_term_structure, float *log_discount_factors, float *d_log_discount_factors, float *d2_log_discount_factors);
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
 );
float NormalVolToLognormalVolFloat(float normal_vol, float forward, float time_to_expiry);

/**
 * Set up to apply the Jamshidian trick to get the swaption vol
 */
float LGM1FGetSwaptionLognormalVol(LGM1FTermStructureFloat *lgm_term_structure, float time_to_expiry, float underlying_length, uint fixed_pay_frequency);
