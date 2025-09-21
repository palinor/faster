//
//  yield_curve.h
//  example_apple_metal
//
//  Created by Aion Feehan on 1/18/24.
//

#include "date.h"

enum class Precision {
	FLOAT,
	DOUBLE,
	COUNT
};
enum class RefRate {
	REF_RATE_ERROR,
	USSOFR,
	USSTERM,
	USLIBOR,
	USCPI
};
enum class PayFreq { PAY_FREQ_ERROR, M, Q, S, A };
enum class Currency { CURRENCY_ERROR, USD, EUR };
struct SwapHeader {
	Precision precision;
};
struct SwapScheduleHeader {
	Precision precision;
};

struct Swap {
	SwapHeader header;
	long id;
	Date start_date;
	Date end_date;
	Date trade_date;
	Time trade_time;
	RefRate ref_rate;
	PayFreq fixed_pay_frequency;
	PayFreq floating_pay_frequency;
	Currency currency;
	enum { ACTION_ERROR, NEW, CANCEL, CORRECT, MODIFY } action_type;
	enum { TRANSACTION_ERROR, TRADE, AMENDMENT, TERMINATION } transaction_type;
	enum { BLOCK_TRADE_ERROR, Y, N } is_block_trade;
	enum { VENUE_ERROR, ON, OFF } venue;
	void *schedule;
};

struct SwapScheduleF {
	SwapScheduleHeader header;
	float notional;
	float fixed_rate;
	float *fixed_payment_times = nullptr;
	float *floating_payment_times = nullptr;
	uint number_of_fixed_payments;
	uint number_of_floating_payments;
};

struct SwapScheduleD {
	SwapScheduleHeader header;
	double notional;
	double fixed_rate;
	double *fixed_payment_times = nullptr;
	double *floating_payment_times = nullptr;
	uint number_of_fixed_payments;
	uint number_of_floating_payments;
};

enum class YieldCurveType {
	SHORT_RATE,
	OVERNIGHT_FORWARD,
	COUNT
};


struct YieldCurveHeader {
	YieldCurveType curve_type;
	Precision precision;
};


struct YieldCurveShortRateFloat {
	YieldCurveHeader header;
	float *interpolation_times = nullptr;
	float *short_rates = nullptr;
	uint number_of_points;
};

struct YieldCurveShortRateDouble {
	YieldCurveHeader header;
	double *interpolation_times = nullptr;
	double *short_rates = nullptr;
	uint number_of_points;
};

struct YieldCurveOvernightForwardFloat {
	YieldCurveHeader header;
	float *interpolation_times = nullptr;
	float *overnight_forwards = nullptr;
	uint number_of_points;
};

struct YieldCurveOvernightForwardDouble {
	YieldCurveHeader header;
	double *interpolation_times = nullptr;
	double *overnight_forwards = nullptr;
	uint number_of_points;
};

void YieldCurveInit(
	void *yield_curve, 
	YieldCurveType curve_type, 
	Precision precision, 
	uint number_of_interpolation_points,
	Arena *arena);

int YieldCurveGetShortRate(void *output_float_or_double, void *yield_curve_input, void *maturity_in_years_float_or_double);
int YieldCurveGetOvernightForward(void *output_float_or_double, void *yield_curve_input, void *maturity_in_years_float_or_double);
int DiscountFactorSpotFromYieldCurve(void *discount_factor_float_or_double, void *yield_curve_input, void *time_to_maturity_in_years);
int DiscountFactorForwardFromYieldCurve(void *discount_factor_float_or_double, void *yield_curve_input, void *time_to_start_in_years, void *time_to_end_in_years);
int YieldCurveShortRateStripFras(
	void *yield_curve_input, 
	void *forward_rates, 
	void *expiry_times_in_years, 
	uint number_of_instruments
	);
int YieldCurveOvernightForwardStripFras(
	void *yield_curve_input, 
	void *forward_rates, 
	void *expiry_times_in_years, 
	uint number_of_instruments
	);
int SwapScheduleInitF(
	void *input, 
	uint number_of_fixed_payments, 
	uint number_of_floating_payments, 
	float swap_rate
	);
int SwapScheduleDInit(
	void *input, 
	uint number_of_fixed_payments, 
	uint number_of_floating_payments, 
	double swap_rate
	);
int YieldCurveShortRateStripSwaps(
	void *yield_curve_input, 
	void *swap_list, 
	void *swap_rates, 
	uint number_of_swaps, 
	uint number_of_already_stripped_fras
	);
int YieldCurveOvernightForwardStripSwaps(
	void *yield_curve_input, 
	void *swap_list, 
	void *swap_rates, 
	uint number_of_swaps, 
	uint number_of_already_stripped_fras
	);
float LevelFromYieldCurveF(void *yield_curve, float time_to_start, float level_maturity, float frequency);
float SwapRateYieldCurveF(void *yield_curve, float time_to_start, float swap_maturity, float fixed_frequency);
