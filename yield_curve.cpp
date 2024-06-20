#include <math.h>
#include <stdlib.h>

#include "date.cpp"

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
	size_t number_of_fixed_payments;
	size_t number_of_floating_payments;
};

struct SwapScheduleD {
	SwapScheduleHeader header;
	double notional;
	double fixed_rate;
	double *fixed_payment_times = nullptr;
	double *floating_payment_times = nullptr;
	size_t number_of_fixed_payments;
	size_t number_of_floating_payments;
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
	size_t number_of_points;
};

struct YieldCurveShortRateDouble {
	YieldCurveHeader header;
	double *interpolation_times = nullptr;
	double *short_rates = nullptr;
	size_t number_of_points;
};

struct YieldCurveOvernightForwardFloat {
	YieldCurveHeader header;
	float *interpolation_times = nullptr;
	float *overnight_forwards = nullptr;
	size_t number_of_points;
};

struct YieldCurveOvernightForwardDouble {
	YieldCurveHeader header;
	double *interpolation_times = nullptr;
	double *overnight_forwards = nullptr;
	size_t number_of_points;
};

void YieldCurveInit(void *yield_curve, YieldCurveType curve_type, Precision precision) {
	YieldCurveHeader *header = (YieldCurveHeader *)yield_curve;
	header->curve_type = curve_type;
	header->precision = precision;
}

int YieldCurveGetShortRate(void *output_float_or_double, void *yield_curve_input, void *maturity_in_years_float_or_double) {
	YieldCurveHeader *header = (YieldCurveHeader *)yield_curve_input;
	switch (header->precision) {
	case Precision::FLOAT: {
		float *output = (float *)output_float_or_double;
		float *maturity_in_years = (float *)maturity_in_years_float_or_double;
		switch (header->curve_type) {
		case YieldCurveType::SHORT_RATE: {
			YieldCurveShortRateFloat *yield_curve = (YieldCurveShortRateFloat *)yield_curve_input;
			float *this_maturity_time = yield_curve->interpolation_times;
			float *this_short_rate = yield_curve->short_rates;
			for (size_t i = 0; i < yield_curve->number_of_points; ++i) {
				if (*this_maturity_time++ < *maturity_in_years) {
					++this_short_rate;
				}
				else {
					*output = *this_short_rate;
					return 0;
				}
			}
		} break;
		case YieldCurveType::OVERNIGHT_FORWARD: {
			YieldCurveOvernightForwardFloat *yield_curve = (YieldCurveOvernightForwardFloat *)yield_curve_input;
			float *this_maturity_time = yield_curve->interpolation_times;
			float *this_overnight_forward = yield_curve->overnight_forwards;
			for (size_t i = 0; i < yield_curve->number_of_points; ++i) {
				if (*this_maturity_time++ < *maturity_in_years) {
					++this_overnight_forward;
				}
				else {
					// todo(AION) there should be coverage here, but for now all days are the same
					*output = logf(1 + *this_overnight_forward);
					return 0;
				}
			}
		} break;
		default: {
			return -1;
		}
		}
	} break;
	case Precision::DOUBLE: {
		double *output = (double *)output_float_or_double;
		double *maturity_in_years = (double *)maturity_in_years_float_or_double;
		switch (header->curve_type) {
		case YieldCurveType::SHORT_RATE: {
			YieldCurveShortRateDouble *yield_curve = (YieldCurveShortRateDouble *)yield_curve_input;
			double *this_maturity_time = yield_curve->interpolation_times;
			double *this_short_rate = yield_curve->short_rates;
			for (size_t i = 0; i < yield_curve->number_of_points; ++i) {
				if (*this_maturity_time++ < *maturity_in_years) {
					++this_short_rate;
				}
				else {
					*output = *this_short_rate;
					return 0;
				}
			}
		} break;
		case YieldCurveType::OVERNIGHT_FORWARD: {
			YieldCurveOvernightForwardDouble *yield_curve = (YieldCurveOvernightForwardDouble *)yield_curve_input;
			double *this_maturity_time = yield_curve->interpolation_times;
			double *this_overnight_forward = yield_curve->overnight_forwards;
			for (size_t i = 0; i < yield_curve->number_of_points; ++i) {
				if (*this_maturity_time++ < *maturity_in_years) {
					++this_overnight_forward;
				}
				else {
					// todo(AION) there should be coverage here, but for now all days are the same
					*output = log(1 + *this_overnight_forward);
					return 0;
				}
			}
		} break;
		default: {
			return -1;
		}
		}
	} break;
	default: {
		return -1;
	}
	}
	return 0;
}

int YieldCurveGetOvernightForward(void *output_float_or_double, void *yield_curve_input, void *maturity_in_years_float_or_double) {
	YieldCurveHeader *header = (YieldCurveHeader *)yield_curve_input;
	switch (header->precision) {
	case Precision::FLOAT: {
		float *output = (float *)output_float_or_double;
		float *maturity_in_years = (float *)maturity_in_years_float_or_double;
		switch (header->curve_type) {
		case YieldCurveType::SHORT_RATE: {
			YieldCurveShortRateFloat *yield_curve = (YieldCurveShortRateFloat *)yield_curve_input;
			float *this_maturity_time = yield_curve->interpolation_times;
			float *this_short_rate = yield_curve->short_rates;
			for (size_t i = 0; i < yield_curve->number_of_points; ++i) {
				if (*this_maturity_time++ < *maturity_in_years) {
					++this_short_rate;
				}
				else {
					*output = expf(*this_short_rate) - 1;
					return 0;
				}
			}
		} break;
		case YieldCurveType::OVERNIGHT_FORWARD: {
			YieldCurveOvernightForwardFloat *yield_curve = (YieldCurveOvernightForwardFloat *)yield_curve_input;
			float *this_maturity_time = yield_curve->interpolation_times;
			float *this_overnight_forward = yield_curve->overnight_forwards;
			for (size_t i = 0; i < yield_curve->number_of_points; ++i) {
				if (*this_maturity_time++ < *maturity_in_years) {
					++this_overnight_forward;
				}
				else {
					// todo(AION) there should be coverage here, but for now all days are the same
					*output = *this_overnight_forward;
					return 0;
				}
			}
		} break;
		default: {
			return -1;
		}
		}
	} break;
	case Precision::DOUBLE: {
		double *output = (double *)output_float_or_double;
		double *maturity_in_years = (double *)maturity_in_years_float_or_double;
		switch (header->curve_type) {
		case YieldCurveType::SHORT_RATE: {
			YieldCurveShortRateDouble *yield_curve = (YieldCurveShortRateDouble *)yield_curve_input;
			double *this_maturity_time = yield_curve->interpolation_times;
			double *this_short_rate = yield_curve->short_rates;
			for (size_t i = 0; i < yield_curve->number_of_points; ++i) {
				if (*this_maturity_time++ < *maturity_in_years) {
					++this_short_rate;
				}
				else {
					*output = exp(*this_short_rate) - 1;
					return 0;
				}
			}
		} break;
		case YieldCurveType::OVERNIGHT_FORWARD: {
			YieldCurveOvernightForwardDouble *yield_curve = (YieldCurveOvernightForwardDouble *)yield_curve_input;
			double *this_maturity_time = yield_curve->interpolation_times;
			double *this_overnight_forward = yield_curve->overnight_forwards;
			for (size_t i = 0; i < yield_curve->number_of_points; ++i) {
				if (*this_maturity_time++ < *maturity_in_years) {
					++this_overnight_forward;
				}
				else {
					// todo(AION) there should be coverage here, but for now all days are the same
					*output = *this_overnight_forward;
					return 0;
				}
			}
		} break;
		default: {
			return -1;
		}
		}
	} break;
	default: {
		return -1;
	}
	}
	return 0;
}

int DiscountFactorSpotFromYieldCurve(void *discount_factor_float_or_double, void *yield_curve_input, void *time_to_maturity_in_years) {
	YieldCurveHeader *header = (YieldCurveHeader *)yield_curve_input;
	switch (header->precision) {
	case Precision::FLOAT: {
		float *output = (float *)discount_factor_float_or_double;
		float *time_to_maturity_in_years_float = (float *)time_to_maturity_in_years;
		switch (header->curve_type) {
		case YieldCurveType::SHORT_RATE: {
			YieldCurveShortRateFloat *yield_curve = (YieldCurveShortRateFloat *)yield_curve_input;
			float *this_maturity_time = yield_curve->interpolation_times;
			float *this_short_rate = yield_curve->short_rates;
			float log_discount_factor = 0;
			for (size_t i = 0; i < yield_curve->number_of_points; ++i) {
				if (*(this_maturity_time + 1) < *time_to_maturity_in_years_float) {
					log_discount_factor += *this_short_rate++ * (*(this_maturity_time + 1) - *this_maturity_time++);
				}
				else {
					log_discount_factor += *this_short_rate * (*time_to_maturity_in_years_float - *this_maturity_time);
					break;
				}
			}
			*output = expf(-log_discount_factor);
			return 0;
		} break;
		case YieldCurveType::OVERNIGHT_FORWARD: {
			YieldCurveOvernightForwardFloat *yield_curve = (YieldCurveOvernightForwardFloat *)yield_curve_input;
			float *this_maturity_time = yield_curve->interpolation_times;
			float *this_overnight_forward = yield_curve->overnight_forwards;
			float discount_factor = 1;
			int time_to_maturity_in_days = (int)(*time_to_maturity_in_years_float * 365);
			for (int i = 0; i < time_to_maturity_in_days; ++i) {
				if ((float)(i) / 365 > *this_maturity_time) {
					++this_maturity_time;
					++this_overnight_forward;
				}
				// todo(AION): this really needs a calendar
				discount_factor /= (1 + *this_overnight_forward / 365);
			}
			*output = discount_factor;
			return 0;
		} break;
		default: return -1;
		}
	}
	case Precision::DOUBLE: {
		double *output = (double *)discount_factor_float_or_double;
		double *time_to_maturity_in_years_float = (double *)time_to_maturity_in_years;
		switch (header->curve_type) {
		case YieldCurveType::SHORT_RATE: {
			YieldCurveShortRateDouble *yield_curve = (YieldCurveShortRateDouble *)yield_curve_input;
			double *this_maturity_time = yield_curve->interpolation_times;
			double *this_short_rate = yield_curve->short_rates;
			double log_discount_factor = 0;
			for (size_t i = 0; i < yield_curve->number_of_points; ++i) {
				if (*(this_maturity_time + 1) < *time_to_maturity_in_years_float) {
					log_discount_factor += *this_short_rate++ * (*(this_maturity_time + 1) - *this_maturity_time++);
				}
				else {
					log_discount_factor += *this_short_rate * (*time_to_maturity_in_years_float - *this_maturity_time);
					break;
				}
			}
			*output = exp(-log_discount_factor);
			return 0;
		} break;
		case YieldCurveType::OVERNIGHT_FORWARD: {
			YieldCurveOvernightForwardDouble *yield_curve = (YieldCurveOvernightForwardDouble *)yield_curve_input;
			double *this_maturity_time = yield_curve->interpolation_times;
			double *this_overnight_forward = yield_curve->overnight_forwards;
			double discount_factor = 1;
			int time_to_maturity_in_days = (int)(*time_to_maturity_in_years_float * 365);
			for (int i = 0; i < time_to_maturity_in_days; ++i) {
				if ((double)(i) / 365 > *this_maturity_time) {
					++this_maturity_time;
					++this_overnight_forward;
				}
				// todo(AION): this really needs a calendar
				discount_factor /= (1 + *this_overnight_forward / 365);
			}
			*output = discount_factor;
			return 0;
		} break;
		default: return -1;
		}
	}
	default: return -1;
	}
	return 0;
}

int DiscountFactorForwardFromYieldCurve(void *discount_factor_float_or_double, void *yield_curve_input, void *time_to_start_in_years, void *time_to_end_in_years) {
	YieldCurveHeader *header = (YieldCurveHeader *)yield_curve_input;
	switch (header->precision) {
	case Precision::FLOAT: {
		float *output = (float *)discount_factor_float_or_double;
		float *time_to_start = (float *)time_to_start_in_years;
		float *time_to_end = (float *)time_to_end_in_years;
		float discount_factor_start, discount_factor_end;
		int error = DiscountFactorSpotFromYieldCurve(&discount_factor_start, yield_curve_input, time_to_start);
		if (error) return -1;
		error = DiscountFactorSpotFromYieldCurve(&discount_factor_end, yield_curve_input, time_to_end);
		if (error) return -1;
		*output = discount_factor_end / discount_factor_start;
		return 0;
	} break;
	case Precision::DOUBLE: {
		double *output = (double *)discount_factor_float_or_double;
		double *time_to_start = (double *)time_to_start_in_years;
		double *time_to_end = (double *)time_to_end_in_years;
		double discount_factor_start;
		double discount_factor_end;
		int error = DiscountFactorSpotFromYieldCurve(&discount_factor_start, yield_curve_input, time_to_start);
		if (error) return -1;
		error = DiscountFactorSpotFromYieldCurve(&discount_factor_end, yield_curve_input, time_to_end);
		if (error) return -1;
		*output = discount_factor_end / discount_factor_start;
		return 0;

	} break;
	default: return -1;
	}
	return -1;
}


int YieldCurveShortRateStripFras(void *yield_curve_input, void *forward_rates, void *expiry_times_in_years, size_t number_of_instruments) {
	YieldCurveHeader *header = (YieldCurveHeader *)yield_curve_input;
	header->curve_type = YieldCurveType::SHORT_RATE;
	switch (header->precision) {
	case Precision::FLOAT: {
		YieldCurveShortRateFloat *yield_curve = (YieldCurveShortRateFloat *)yield_curve_input;
		float first_expiry_time = 0;
		float second_expiry_time = 0;
		if ((yield_curve->short_rates) && (yield_curve->number_of_points != number_of_instruments)) {
			void *new_block = realloc(yield_curve->short_rates, number_of_instruments * sizeof(float));
			if (!new_block) return -1;
			yield_curve->short_rates = (float *)new_block;
		}
		if ((yield_curve->interpolation_times) && (yield_curve->number_of_points != number_of_instruments)) {
			void *new_block = realloc(yield_curve->interpolation_times, number_of_instruments * sizeof(float));
			if (!new_block) return -1;
			yield_curve->interpolation_times = (float *)new_block;
		}
		yield_curve->number_of_points = number_of_instruments;
		if (!yield_curve->short_rates) yield_curve->short_rates = (float *)malloc(number_of_instruments * sizeof(float));
		if (!yield_curve->short_rates) return -1;
		if (!yield_curve->interpolation_times) yield_curve->interpolation_times = (float *)malloc(number_of_instruments * sizeof(float));
		if (!yield_curve->interpolation_times) return -1;

		float *this_forward_rate = (float *)forward_rates;
		float *this_expiry_time = (float *)expiry_times_in_years;
		float *yield_curve_short_rate = yield_curve->short_rates;
		float *yield_curve_expiry_time = yield_curve->interpolation_times;
		for (size_t i = 0; i < number_of_instruments; ++i) {
			second_expiry_time = *this_expiry_time;
			*yield_curve_expiry_time++ = *this_expiry_time++;
			float coverage = second_expiry_time - first_expiry_time;
			*yield_curve_short_rate++ = logf(1 + coverage * (*this_forward_rate++)) / coverage;
			first_expiry_time = second_expiry_time;
		}
	} break;
	case Precision::DOUBLE: {
		YieldCurveShortRateDouble *yield_curve = (YieldCurveShortRateDouble *)yield_curve_input;
		double first_expiry_time = 0;
		double second_expiry_time = 0;
		if ((yield_curve->short_rates) && (yield_curve->number_of_points != number_of_instruments)) {
			void *new_block = realloc(yield_curve->short_rates, number_of_instruments * sizeof(double));
			if (!new_block) return -1;
			yield_curve->short_rates = (double *)new_block;
		}
		if ((yield_curve->interpolation_times) && (yield_curve->number_of_points != number_of_instruments)) {
			void *new_block = realloc(yield_curve->interpolation_times, number_of_instruments * sizeof(double));
			if (!new_block) return -1;
			yield_curve->interpolation_times = (double *)new_block;
		}
		yield_curve->number_of_points = number_of_instruments;
		if (!yield_curve->short_rates) yield_curve->short_rates = (double *)malloc(number_of_instruments * sizeof(double));
		if (!yield_curve->short_rates) return -1;
		if (!yield_curve->interpolation_times) yield_curve->interpolation_times = (double *)malloc(number_of_instruments * sizeof(double));
		if (!yield_curve->interpolation_times) return -1;

		double *this_forward_rate = (double *)forward_rates;
		double *this_expiry_time = (double *)expiry_times_in_years;
		double *yield_curve_short_rate = yield_curve->short_rates;
		double *yield_curve_expiry_time = yield_curve->interpolation_times;
		for (size_t i = 0; i < number_of_instruments; ++i) {
			second_expiry_time = *this_expiry_time;
			*yield_curve_expiry_time++ = *this_expiry_time++;
			double coverage = second_expiry_time - first_expiry_time;
			*yield_curve_short_rate++ = log(1 + coverage * (*this_forward_rate++)) / coverage;
			first_expiry_time = second_expiry_time;
		}
	} break;
	default: return -1;
	}
	return 0;
}

int YieldCurveOvernightForwardStripFras(void *yield_curve_input, void *forward_rates, void *expiry_times_in_years, size_t number_of_instruments) {
	YieldCurveHeader *header = (YieldCurveHeader *)yield_curve_input;
	header->curve_type = YieldCurveType::OVERNIGHT_FORWARD;
	switch (header->precision) {
	case Precision::FLOAT: {
		YieldCurveOvernightForwardFloat *yield_curve = (YieldCurveOvernightForwardFloat *)yield_curve_input;
		float first_expiry_time = 0;
		float second_expiry_time = 0;
		if ((yield_curve->overnight_forwards) && (yield_curve->number_of_points != number_of_instruments)) {
			yield_curve->overnight_forwards = (float *)realloc(yield_curve->overnight_forwards, number_of_instruments * sizeof(float));
			if (!yield_curve->overnight_forwards) return -1;
		}
		if ((yield_curve->interpolation_times) && (yield_curve->number_of_points != number_of_instruments)) {
			yield_curve->interpolation_times = (float *)realloc(yield_curve->interpolation_times, number_of_instruments * sizeof(float));
			if (!yield_curve->interpolation_times) return -1;
		}
		yield_curve->number_of_points = number_of_instruments;
		if (!yield_curve->overnight_forwards) yield_curve->overnight_forwards = (float *)malloc(number_of_instruments * sizeof(float));
		if (!yield_curve->overnight_forwards) return -1;
		if (!yield_curve->interpolation_times) yield_curve->interpolation_times = (float *)malloc(number_of_instruments * sizeof(float));
		if (!yield_curve->interpolation_times) return -1;

		float *this_forward_rate = (float *)forward_rates;
		float *this_expiry_time = (float *)expiry_times_in_years;
		float *yield_curve_overnight_forward = yield_curve->overnight_forwards;
		float *yield_curve_expiry_time = yield_curve->interpolation_times;
		for (size_t i = 0; i < number_of_instruments; ++i) {
			second_expiry_time = *this_expiry_time;
			*yield_curve_expiry_time++ = *this_expiry_time++;
			float coverage = second_expiry_time - first_expiry_time;
			float coverage_in_days = coverage * 365;
			*yield_curve_overnight_forward++ = 365 * (expf(logf(1 + coverage * (*this_forward_rate++)) / coverage_in_days) - 1);
			first_expiry_time = second_expiry_time;
		}
	} break;
	case Precision::DOUBLE: {
		YieldCurveOvernightForwardDouble *yield_curve = (YieldCurveOvernightForwardDouble *)yield_curve_input;
		double first_expiry_time = 0;
		double second_expiry_time = 0;
		if ((yield_curve->overnight_forwards) && (yield_curve->number_of_points != number_of_instruments)) {
			void *new_block = realloc(yield_curve->overnight_forwards, number_of_instruments * sizeof(double));
			if (!new_block) return -1;
			yield_curve->overnight_forwards = (double *)new_block;
		}
		if ((yield_curve->interpolation_times) && (yield_curve->number_of_points != number_of_instruments)) {
			void *new_block = realloc(yield_curve->interpolation_times, number_of_instruments * sizeof(double));
			if (!new_block) return -1;
			yield_curve->interpolation_times = (double *)new_block;
		}
		yield_curve->number_of_points = number_of_instruments;
		if (!yield_curve->overnight_forwards) yield_curve->overnight_forwards = (double *)malloc(number_of_instruments * sizeof(double));
		if (!yield_curve->interpolation_times) yield_curve->interpolation_times = (double *)malloc(number_of_instruments * sizeof(double));

		double *this_forward_rate = (double *)forward_rates;
		double *this_expiry_time = (double *)expiry_times_in_years;
		double *yield_curve_overnight_forward = yield_curve->overnight_forwards;
		double *yield_curve_expiry_time = yield_curve->interpolation_times;
		for (size_t i = 0; i < number_of_instruments; ++i) {
			second_expiry_time = *this_expiry_time;
			*yield_curve_expiry_time++ = *this_expiry_time++;
			double coverage = second_expiry_time - first_expiry_time;
			double coverage_in_days = coverage * 365;
			*yield_curve_overnight_forward++ = 365 * (exp(log(1 + coverage * (*this_forward_rate++)) / coverage_in_days) - 1);
			first_expiry_time = second_expiry_time;
		}
	} break;
	default: return -1;
	}
	return 0;
}


int SwapScheduleInitF(SwapScheduleF *swap_schedule, size_t number_of_fixed_payments, size_t number_of_floating_payments, float swap_rate) {
	swap_schedule->header.precision = Precision::FLOAT;
	swap_schedule->fixed_rate = swap_rate;
	if ((swap_schedule->fixed_payment_times) && (swap_schedule->number_of_fixed_payments != number_of_fixed_payments)) {
		swap_schedule->fixed_payment_times = (float *)realloc(swap_schedule->fixed_payment_times, number_of_fixed_payments * sizeof(float));
	}
	if (!swap_schedule->fixed_payment_times) swap_schedule->fixed_payment_times = (float *)malloc(number_of_fixed_payments * sizeof(float));
	if (!swap_schedule->fixed_payment_times) return -1;

	if ((swap_schedule->floating_payment_times) && (swap_schedule->number_of_floating_payments != number_of_floating_payments)) {
		swap_schedule->floating_payment_times = (float *)realloc(swap_schedule->floating_payment_times, number_of_floating_payments * sizeof(float));
	}
	if (!swap_schedule->floating_payment_times) swap_schedule->floating_payment_times = (float *)malloc(number_of_floating_payments * sizeof(float));
	if (!swap_schedule->floating_payment_times) return -1;
	swap_schedule->number_of_fixed_payments = number_of_fixed_payments;
	swap_schedule->number_of_floating_payments = number_of_floating_payments;
	float *fixed_payment_time = swap_schedule->fixed_payment_times;
	float *floating_payment_time = swap_schedule->floating_payment_times;
	for (size_t i = 0; i < number_of_fixed_payments; ++i) {
		*fixed_payment_time++ = (float)(i + 1);
	}
	for (size_t i = 0; i < number_of_floating_payments; ++i) {
		*floating_payment_time++ = (float)(i + 1);
	}
	return 0;
}

int SwapScheduleDInit(SwapScheduleD *swap_schedule, size_t number_of_fixed_payments, size_t number_of_floating_payments, double swap_rate) {
	swap_schedule->header.precision = Precision::DOUBLE;
	swap_schedule->fixed_rate = swap_rate;
	if ((swap_schedule->fixed_payment_times) && (swap_schedule->number_of_fixed_payments != number_of_fixed_payments)) {
		swap_schedule->fixed_payment_times = (double *)realloc(swap_schedule->fixed_payment_times, number_of_fixed_payments * sizeof(double));
	}
	if (!swap_schedule->fixed_payment_times) swap_schedule->fixed_payment_times = (double *)malloc(number_of_fixed_payments * sizeof(double));
	if (!swap_schedule->fixed_payment_times) return -1;

	if ((swap_schedule->floating_payment_times) && (swap_schedule->number_of_floating_payments != number_of_floating_payments)) {
		swap_schedule->floating_payment_times = (double *)realloc(swap_schedule->floating_payment_times, number_of_floating_payments * sizeof(double));
	}
	if (!swap_schedule->floating_payment_times) swap_schedule->floating_payment_times = (double *)malloc(number_of_floating_payments * sizeof(double));
	if (!swap_schedule->floating_payment_times) return -1;

	swap_schedule->number_of_fixed_payments = number_of_fixed_payments;
	swap_schedule->number_of_floating_payments = number_of_floating_payments;
	double *fixed_payment_time = swap_schedule->fixed_payment_times;
	double *floating_payment_time = swap_schedule->floating_payment_times;
	for (size_t i = 0; i < number_of_fixed_payments; ++i) {
		*fixed_payment_time++ = (double)(i + 1);
	}
	for (size_t i = 0; i < number_of_floating_payments; ++i) {
		*floating_payment_time++ = (double)(i + 1);
	}
	return 0;
}

inline void SwapScheduleFFree(SwapScheduleF *swap) {
	free(swap->fixed_payment_times);
	free(swap->floating_payment_times);
}

inline void SwapScheduleDFree(SwapScheduleD *swap) {
	free(swap->fixed_payment_times);
	free(swap->floating_payment_times);
}

int YieldCurveShortRateStripSwaps(void *yield_curve_input, void *swap_list, void *swap_rates, size_t number_of_swaps, size_t number_of_already_stripped_fras) {
	YieldCurveHeader *header = (YieldCurveHeader *)yield_curve_input;
	switch (header->precision) {
	case Precision::FLOAT: {
		YieldCurveShortRateFloat *yield_curve = (YieldCurveShortRateFloat *)yield_curve_input;
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
		float *this_swap_rate = (float *)swap_rates;
		float *previous_swap_rate = (float *)swap_rates;
		float *yield_curve_short_rate = yield_curve->short_rates + number_of_already_stripped_fras;
		float *yield_curve_expiry_time = yield_curve->interpolation_times + number_of_already_stripped_fras;
		float fixed_leg_sum = 0;
		size_t number_of_known_swap_discount_factors = 0;
		SwapScheduleF *first_swap = (SwapScheduleF *)swap_list;
		SwapScheduleF *last_swap = first_swap + number_of_swaps - 1;
		float last_fra_payment_time = yield_curve->interpolation_times[number_of_already_stripped_fras - 1];
		float *this_swap_payment_time = first_swap->fixed_payment_times;
		for (size_t i = 0; i < first_swap->number_of_fixed_payments; ++i) {
			if (last_fra_payment_time > *this_swap_payment_time++) ++number_of_known_swap_discount_factors;
		}
		this_swap_payment_time = first_swap->fixed_payment_times;
		float last_swap_payment_time = 0;
		float this_discount_factor = 1;
		for (size_t i = 0; i < number_of_known_swap_discount_factors; ++i) {
			float coverage = (*this_swap_payment_time - last_swap_payment_time);
			int err = DiscountFactorSpotFromYieldCurve(&this_discount_factor, yield_curve, this_swap_payment_time);
			if (err) return -1;
			fixed_leg_sum += coverage * this_discount_factor;
			last_swap_payment_time = *this_swap_payment_time;
			++this_swap_payment_time;
		}
		float last_discount_factor = this_discount_factor;
		float last_swap_rate = 1;
		SwapScheduleF *this_swap = first_swap;
		for (size_t i = 0; i < number_of_swaps; ++i) {
			this_swap_payment_time = this_swap->fixed_payment_times + this_swap->number_of_fixed_payments - 1;
			fixed_leg_sum *= *this_swap_rate / last_swap_rate;
			float coverage = (*this_swap_payment_time - last_swap_payment_time);
			this_discount_factor = (1 - fixed_leg_sum) / (1 + coverage * (*this_swap_rate));
			*yield_curve_expiry_time++ = *this_swap_payment_time;
			*yield_curve_short_rate++ = logf(last_discount_factor / this_discount_factor) / coverage;
			fixed_leg_sum += coverage * this_discount_factor * (*this_swap_rate);
			last_swap_rate = *this_swap_rate++;
			last_swap_payment_time = *this_swap_payment_time;
			++this_swap;
		}

	} break;
	case Precision::DOUBLE: {
		YieldCurveShortRateDouble *yield_curve = (YieldCurveShortRateDouble *)yield_curve_input;
		double first_expiry_time = 0;
		double second_expiry_time = 0;
		size_t number_of_total_instruments = number_of_swaps + number_of_already_stripped_fras;
		if ((yield_curve->short_rates) && (yield_curve->number_of_points != number_of_total_instruments)) {
			yield_curve->short_rates = (double *)realloc(yield_curve->short_rates, number_of_total_instruments * sizeof(double));
		}
		if ((yield_curve->interpolation_times) && (yield_curve->number_of_points != number_of_total_instruments)) {
			yield_curve->interpolation_times = (double *)realloc(yield_curve->interpolation_times, number_of_total_instruments * sizeof(double));
		}
		yield_curve->number_of_points = number_of_total_instruments;
		if (!yield_curve->short_rates) {
			if (number_of_already_stripped_fras > 0) {
				return -1;
			}
			else {
				yield_curve->short_rates = (double *)malloc(number_of_swaps * sizeof(double));
			}
		}
		if (!yield_curve->interpolation_times) {
			if (number_of_already_stripped_fras > 0) {
				return -1;
			}
			else {
				yield_curve->interpolation_times = (double *)malloc(number_of_swaps * sizeof(double));
			}
		}
		double *this_swap_rate = (double *)swap_rates;
		double *previous_swap_rate = (double *)swap_rates;
		double *yield_curve_short_rate = yield_curve->short_rates + number_of_already_stripped_fras;
		double *yield_curve_expiry_time = yield_curve->interpolation_times + number_of_already_stripped_fras;
		double fixed_leg_sum = 0;
		size_t number_of_known_swap_discount_factors = 0;
		SwapScheduleD *first_swap = (SwapScheduleD *)swap_list;
		SwapScheduleD *last_swap = first_swap + number_of_swaps - 1;
		double last_fra_payment_time = yield_curve->interpolation_times[number_of_already_stripped_fras - 1];
		double *this_swap_payment_time = first_swap->fixed_payment_times;
		for (size_t i = 0; i < first_swap->number_of_fixed_payments; ++i) {
			if (last_fra_payment_time > *this_swap_payment_time++) ++number_of_known_swap_discount_factors;
		}
		this_swap_payment_time = first_swap->fixed_payment_times;
		double last_swap_payment_time = 0;
		double this_discount_factor = 1;
		for (size_t i = 0; i < number_of_known_swap_discount_factors; ++i) {
			double coverage = (*this_swap_payment_time - last_swap_payment_time);
			int err = DiscountFactorSpotFromYieldCurve(&this_discount_factor, yield_curve, this_swap_payment_time);
			if (err) return -1;
			fixed_leg_sum += coverage * this_discount_factor;
			last_swap_payment_time = *this_swap_payment_time;
			++this_swap_payment_time;
		}
		double last_discount_factor = this_discount_factor;
		double last_swap_rate = 1;
		SwapScheduleD *this_swap = first_swap;
		for (size_t i = 0; i < number_of_swaps; ++i) {
			this_swap_payment_time = this_swap->fixed_payment_times + this_swap->number_of_fixed_payments - 1;
			fixed_leg_sum *= *this_swap_rate / last_swap_rate;
			double coverage = (*this_swap_payment_time - last_swap_payment_time);
			this_discount_factor = (1 - fixed_leg_sum) / (1 + coverage * (*this_swap_rate));
			*yield_curve_expiry_time++ = *this_swap_payment_time;
			*yield_curve_short_rate++ = log(last_discount_factor / this_discount_factor) / coverage;
			fixed_leg_sum += coverage * this_discount_factor * (*this_swap_rate);
			last_swap_rate = *this_swap_rate++;
			last_swap_payment_time = *this_swap_payment_time;
			++this_swap;
		}
	} break;
	default: return -1;
	}
	return 0;
}

int YieldCurveOvernightForwardStripSwaps(void *yield_curve_input, void *swap_list, void *swap_rates, size_t number_of_swaps, size_t number_of_already_stripped_fras) {
	YieldCurveHeader *header = (YieldCurveHeader *)yield_curve_input;
	switch (header->precision) {
	case Precision::FLOAT: {
		YieldCurveOvernightForwardFloat *yield_curve = (YieldCurveOvernightForwardFloat *)yield_curve_input;
		float first_expiry_time = 0;
		float second_expiry_time = 0;
		size_t number_of_total_instruments = number_of_swaps + number_of_already_stripped_fras;
		if ((yield_curve->overnight_forwards) && (yield_curve->number_of_points != number_of_total_instruments)) {
			yield_curve->overnight_forwards = (float *)realloc(yield_curve->overnight_forwards, number_of_total_instruments * sizeof(float));
		}
		if ((yield_curve->interpolation_times) && (yield_curve->number_of_points != number_of_total_instruments)) {
			yield_curve->interpolation_times = (float *)realloc(yield_curve->interpolation_times, number_of_total_instruments * sizeof(float));
		}
		yield_curve->number_of_points = number_of_total_instruments;
		if (!yield_curve->overnight_forwards) {
			if (number_of_already_stripped_fras > 0) {
				return -1;
			}
			else {
				yield_curve->overnight_forwards = (float *)malloc(number_of_swaps * sizeof(float));
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
		float *this_swap_rate = (float *)swap_rates;
		float *previous_swap_rate = (float *)swap_rates;
		float *yield_curve_overnight_forward = yield_curve->overnight_forwards + number_of_already_stripped_fras;
		float *yield_curve_expiry_time = yield_curve->interpolation_times + number_of_already_stripped_fras;
		float fixed_leg_sum = 0;
		size_t number_of_known_swap_discount_factors = 0;
		SwapScheduleF *first_swap = (SwapScheduleF *)swap_list;
		SwapScheduleF *last_swap = first_swap + number_of_swaps - 1;
		float last_fra_payment_time = yield_curve->interpolation_times[number_of_already_stripped_fras - 1];
		float *this_swap_payment_time = first_swap->fixed_payment_times;
		for (size_t i = 0; i < first_swap->number_of_fixed_payments; ++i) {
			if (last_fra_payment_time > *this_swap_payment_time++) ++number_of_known_swap_discount_factors;
		}
		this_swap_payment_time = first_swap->fixed_payment_times;
		float last_swap_payment_time = 0;
		float this_discount_factor = 1;
		for (size_t i = 0; i < number_of_known_swap_discount_factors; ++i) {
			float coverage = (*this_swap_payment_time - last_swap_payment_time);
			int err = DiscountFactorSpotFromYieldCurve(&this_discount_factor, yield_curve, this_swap_payment_time);
			if (err) return -1;
			fixed_leg_sum += coverage * this_discount_factor;
			last_swap_payment_time = *this_swap_payment_time;
			++this_swap_payment_time;
		}
		float last_discount_factor = this_discount_factor;
		float last_swap_rate = 1;
		SwapScheduleF *this_swap = first_swap;
		for (size_t i = 0; i < number_of_swaps; ++i) {
			this_swap_payment_time = this_swap->fixed_payment_times + this_swap->number_of_fixed_payments - 1;
			fixed_leg_sum *= *this_swap_rate / last_swap_rate;
			float coverage = (*this_swap_payment_time - last_swap_payment_time);
			this_discount_factor = (1 - fixed_leg_sum) / (1 + coverage * (*this_swap_rate));
			*yield_curve_expiry_time++ = *this_swap_payment_time;
			float coverage_in_days = coverage * 365;
			*yield_curve_overnight_forward++ = 365 * (expf(logf(last_discount_factor / this_discount_factor) / coverage_in_days) - 1);
			fixed_leg_sum += coverage * this_discount_factor * (*this_swap_rate);
			last_swap_rate = *this_swap_rate++;
			last_swap_payment_time = *this_swap_payment_time;
			++this_swap;
		}

	} break;
	case Precision::DOUBLE: {
		YieldCurveOvernightForwardDouble *yield_curve = (YieldCurveOvernightForwardDouble *)yield_curve_input;
		double first_expiry_time = 0;
		double second_expiry_time = 0;
		size_t number_of_total_instruments = number_of_swaps + number_of_already_stripped_fras;
		if ((yield_curve->overnight_forwards) && (yield_curve->number_of_points != number_of_total_instruments)) {
			yield_curve->overnight_forwards = (double *)realloc(yield_curve->overnight_forwards, number_of_total_instruments * sizeof(double));
		}
		if ((yield_curve->interpolation_times) && (yield_curve->number_of_points != number_of_total_instruments)) {
			yield_curve->interpolation_times = (double *)realloc(yield_curve->interpolation_times, number_of_total_instruments * sizeof(double));
		}
		yield_curve->number_of_points = number_of_total_instruments;
		if (!yield_curve->overnight_forwards) {
			if (number_of_already_stripped_fras > 0) {
				return -1;
			}
			else {
				yield_curve->overnight_forwards = (double *)malloc(number_of_swaps * sizeof(double));
			}
		}
		if (!yield_curve->interpolation_times) {
			if (number_of_already_stripped_fras > 0) {
				return -1;
			}
			else {
				yield_curve->interpolation_times = (double *)malloc(number_of_swaps * sizeof(double));
			}
		}
		double *this_swap_rate = (double *)swap_rates;
		double *previous_swap_rate = (double *)swap_rates;
		double *yield_curve_overnight_forward = yield_curve->overnight_forwards + number_of_already_stripped_fras;
		double *yield_curve_expiry_time = yield_curve->interpolation_times + number_of_already_stripped_fras;
		double fixed_leg_sum = 0;
		size_t number_of_known_swap_discount_factors = 0;
		SwapScheduleD *first_swap = (SwapScheduleD *)swap_list;
		SwapScheduleD *last_swap = first_swap + number_of_swaps - 1;
		double last_fra_payment_time = yield_curve->interpolation_times[number_of_already_stripped_fras - 1];
		double *this_swap_payment_time = first_swap->fixed_payment_times;
		for (size_t i = 0; i < first_swap->number_of_fixed_payments; ++i) {
			if (last_fra_payment_time > *this_swap_payment_time++) ++number_of_known_swap_discount_factors;
		}
		this_swap_payment_time = first_swap->fixed_payment_times;
		double last_swap_payment_time = 0;
		double this_discount_factor = 1;
		for (size_t i = 0; i < number_of_known_swap_discount_factors; ++i) {
			double coverage = (*this_swap_payment_time - last_swap_payment_time);
			int err = DiscountFactorSpotFromYieldCurve(&this_discount_factor, yield_curve, this_swap_payment_time);
			if (err) return -1;
			fixed_leg_sum += coverage * this_discount_factor;
			last_swap_payment_time = *this_swap_payment_time;
			++this_swap_payment_time;
		}
		double last_discount_factor = this_discount_factor;
		double last_swap_rate = 1;
		SwapScheduleD *this_swap = first_swap;
		for (size_t i = 0; i < number_of_swaps; ++i) {
			this_swap_payment_time = this_swap->fixed_payment_times + this_swap->number_of_fixed_payments - 1;
			fixed_leg_sum *= *this_swap_rate / last_swap_rate;
			double coverage = (*this_swap_payment_time - last_swap_payment_time);
			this_discount_factor = (1 - fixed_leg_sum) / (1 + coverage * (*this_swap_rate));
			*yield_curve_expiry_time++ = *this_swap_payment_time;
			double coverage_in_days = coverage * 365;
			*yield_curve_overnight_forward++ = 365 * (exp(log(last_discount_factor / this_discount_factor) / coverage_in_days) - 1);
			fixed_leg_sum += coverage * this_discount_factor * (*this_swap_rate);
			last_swap_rate = *this_swap_rate++;
			last_swap_payment_time = *this_swap_payment_time;
			++this_swap;
		}
	} break;
	default: return -1;
	}
	return 0;
}

float LevelFromYieldCurveF(void *yield_curve, float time_to_start, float level_maturity, float frequency) {
	float result = 0;
	for (float discount_factor_length = 1 / frequency; discount_factor_length <= level_maturity; discount_factor_length += 1 / frequency) {
		float df;
		float time_to_end = time_to_start + discount_factor_length;
		DiscountFactorForwardFromYieldCurve(&df, yield_curve, &time_to_start, &time_to_end);
		result += df / frequency;
	}
	return result;
}

float SwapRateYieldCurveF(void *yield_curve, float time_to_start, float swap_maturity, float floating_frequency, float fixed_frequency) {
	float level = LevelFromYieldCurveF(yield_curve, time_to_start, swap_maturity, fixed_frequency);
	float df_start, df_end;
	float time_to_end = time_to_start + swap_maturity;
	DiscountFactorSpotFromYieldCurve(&df_start, yield_curve, &time_to_start);
	DiscountFactorSpotFromYieldCurve(&df_end, yield_curve, &time_to_end);
	float floating_leg_pv = df_start - df_end;
	return floating_leg_pv / level;
}
