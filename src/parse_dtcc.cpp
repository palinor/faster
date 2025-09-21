#include <cerrno>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define CURL_STATICLIB
#include "curl/curl.h"
#include "yield_curve.cpp"

#pragma comment (lib, "libcurl_a.lib")
#pragma comment (lib, "Ws2_32.lib")
#pragma comment (lib, "Wldap32.lib")
#pragma comment (lib, "Normaliz.lib")
#pragma comment (lib, "Crypt32.lib")
#pragma comment (lib, "advapi32.lib")

typedef unsigned char uchar;
typedef unsigned int uint;

struct StringBuffer {
	size_t capacity;
	size_t length;
	char *string = nullptr;
};


// string buffer
#define MAX_STRING_SIZE 4096

static const size_t date_string_size = 10;

inline size_t MinSizeT(size_t a, size_t b) {
	return a > b ? b : a;
}

inline size_t MaxSizeT(size_t a, size_t b) {
	return a > b ? a : b;
}


int StringBufferInit(StringBuffer *buffer) {
	buffer->capacity = 32;
	buffer->length = 0;
	if (buffer->string) free(buffer->string);
	buffer->string = (char *)malloc(32);
	if (!buffer->string) {
		return 1;
	}
	buffer->string[0] = '\0';
	return 0;
}

int StringBufferResize(StringBuffer *buffer) {
	size_t new_capacity = MinSizeT(2 * (buffer->capacity), MAX_STRING_SIZE);
	buffer->capacity = new_capacity;
	char *new_block = (char *)realloc(buffer->string, new_capacity);
	if (!new_block) {
		free(buffer->string);
		return 1;
	}
	buffer->string = new_block;
	return 0;
}

int StringBufferAppendChar(StringBuffer *buffer, const char *new_char) {
	if (buffer->capacity + 1 > MAX_STRING_SIZE) return 1;
	if (buffer->length + 1 > buffer->capacity) StringBufferResize(buffer);
	buffer->string[buffer->length++] = *new_char;
	buffer->string[buffer->length] = '\0';
	return 0;
}

int StringBufferAppendStringWithSize(StringBuffer *buffer, const char *new_string, size_t new_string_size) {
	if (buffer->capacity + new_string_size > MAX_STRING_SIZE) return 1;
	while (buffer->length + new_string_size > buffer->capacity) StringBufferResize(buffer);
	memcpy(buffer->string + buffer->length, new_string, new_string_size);
	buffer->length += new_string_size;
	buffer->string[buffer->length] = '\0';
	return 0;
}

int StringBufferAppendString(StringBuffer *buffer, const char *new_string) {
	size_t new_string_size = strlen(new_string);
	if (buffer->capacity + new_string_size > MAX_STRING_SIZE) return 1;
	while (buffer->length + new_string_size > buffer->capacity) StringBufferResize(buffer);
	memcpy(buffer->string + buffer->length, new_string, new_string_size);
	buffer->length += new_string_size;
	buffer->string[buffer->length] = '\0';
	return 0;
}

int StringBufferAppendStringBuffer(StringBuffer *buffer, StringBuffer *new_string_buffer) {
	if (buffer->capacity + new_string_buffer->length > MAX_STRING_SIZE) return 1;
	while (buffer->length + new_string_buffer->length > buffer->capacity) StringBufferResize(buffer);
	memcpy(buffer->string + buffer->length, new_string_buffer->string, new_string_buffer->length);
	buffer->length += new_string_buffer->length;
	buffer->string[buffer->length] = '\0';
	return 0;
}

void StringBufferClip(StringBuffer *buffer, size_t elements_to_remove) {
	buffer->length -= elements_to_remove;
	buffer->string[buffer->length] = '\0';
}

inline void StringBufferReset(StringBuffer *buffer) {
	buffer->length = 0;
}

void StringBufferClear(StringBuffer *buffer) {
	if (buffer->length > 0) free(buffer->string);
}

/*** Date parsing ***/
#define DATE_STR_LEN 11      // strlen("2022-09-10") + '\0'
#define DATETIME_STR_LEN 21  // strlen("2022-09-10T20:15:56") + '\0'
#define DATETIME_MONTH_OFFSET 5
#define DATETIME_DAY_OFFSET 8
#define DATETIME_HOUR_OFFSET 11
#define DATETIME_MINUTE_OFFSET 14
#define DATETIME_SECOND_OFFSET 17


int DateToString(char *output_string_buffer, Date *input_date) {
	snprintf(output_string_buffer, 5, "%i", input_date->year);
	output_string_buffer[4] = '_';
	if (input_date->month < 10) {
		snprintf(output_string_buffer + 5, 2, "0");
		snprintf(output_string_buffer + 6, 2, "%i", input_date->month);
	}
	else {
		snprintf(output_string_buffer + 5, 3, "%i", input_date->month);
	}
	output_string_buffer[7] = '_';
	if (input_date->day < 10) {
		snprintf(output_string_buffer + 8, 2, "0");
		snprintf(output_string_buffer + 9, 2, "%i", input_date->day);
	}
	else {
		snprintf(output_string_buffer + 8, 3, "%i", input_date->day);
	}
	output_string_buffer[10] = '\0';
	return 0;
}



void StripCommas(char *output, const char *input, size_t input_size) {
	int write_idx = 0;
	int read_idx = 0;
	char this_char = 0;
	for (size_t i = 0; i < input_size; i++) {
		if ((this_char = input[read_idx++]) != ',') output[write_idx++] = this_char;
	}
}

long HandleStrtol(const char *input) {
	// strip out commas
	char buffer[32];
	size_t input_len = strlen(input);
	StripCommas(buffer, input, input_len + 1);
	char *endptr = NULL;
	long res = strtol(buffer, &endptr, 10);
	if (res == 0) {
		if (errno != 0) {
			char err_msg[64];
			sprintf(err_msg, 
			    "Encountered error in strtol, errno : %d, input: %s",
				errno, input);
			// todo(AION) figure out what to do with error messages
			// log(err_msg);
		}
		else if (buffer == endptr) {
			char err_msg[64];
			sprintf(err_msg,
				"Encountered error in strtol, no characters read from input: %s",
				input);
			// log(err_msg);
		}
	}
	return res;
}

float HandleStrtof(const char *input) {
	char buffer[32];
	size_t input_len = strlen(input);
	StripCommas(buffer, input, input_len + 1);
	char *endptr = NULL;
	float res = strtof(buffer, &endptr);
	if (res == 0) {
		if (errno != 0) {
			char err_msg[64];
			sprintf(err_msg,
				"Encountered error in strtof, errno : %d, input: %s",
				errno, input);
			// log(err_msg);
		}
		else if (buffer == endptr) {
			char err_msg[64];
			sprintf(err_msg,
				"Encountered error in strtof, no characters read from input: %s",
				input);
			// log(err_msg);
		}
	}
	return res;
}

double HandleStrtod(const char *input) {
	char buffer[32];
	size_t input_len = strlen(input);
	StripCommas(buffer, input, input_len + 1);
	char *endptr = NULL;
	float res = strtod(buffer, &endptr);
	if (res == 0) {
		if (errno != 0) {
			char err_msg[64];
			sprintf(err_msg,
				"Encountered error in strtod, errno : %d, input: %s",
				errno, input);
			// log(err_msg);
		}
		else if (buffer == endptr) {
			char err_msg[64];
			sprintf(err_msg,
				"Encountered error in strtod, no characters read from input: %s",
				input);
			// log(err_msg);
		}
	}
	return res;
}

// assume date is in format YYYY-MM-DD
int ParseDate(const char *date_string, Date *result_date) {
	size_t expected_len = strlen("YYYY-MM-DD");
	if (strlen(date_string) != expected_len) return 1;
	char year_c[5], month_c[3], day_c[3];
	strncpy(year_c, date_string, 4);
	strncpy(month_c, date_string + DATETIME_MONTH_OFFSET, 2);
	strncpy(day_c, date_string + DATETIME_DAY_OFFSET, 2);
	year_c[4] = '\0';
	month_c[2] = '\0';
	day_c[2] = '\0';
	result_date->year = HandleStrtol(year_c);
	result_date->month = HandleStrtol(month_c);
	result_date->day = HandleStrtol(day_c);
	return 0;
}

// handle how case where we want 1 to be cast to the string "01"
int CastIntTo2Chars(char *output, int input) {
	if ((input < 0) || (input > 99)) {
		return 1;
	}
	else if (input < 10) {
		output[0] = '0';
		sprintf(output + 1, "%d", input);
	}
	else {
		sprintf(output, "%d", input);
	}
	return 0;
}

// take struct tm and turn it back into a string like YYYY-MM-DD
int DateStringFromDate(char *output_string, Date *input_date) {
	int err = 0;
	sprintf(output_string, "%d", input_date->year);
	err = CastIntTo2Chars(output_string + 5, input_date->month);
	if (err) return err;
	err = CastIntTo2Chars(output_string + 8, input_date->day);
	if (err) return err;
	output_string[DATETIME_MONTH_OFFSET - 1] = '-';
	output_string[DATETIME_DAY_OFFSET - 1] = '-';
	output_string[DATE_STR_LEN - 1] = '\0';
	return 0;
}

// assume datetime is in format YYYY-MM-DDTHH:MM:SS
int ParseDatetime(const char *datetime, Date *result_date, Time *result_time) {
	size_t expected_len = strlen("YYYY-MM-DDTHH:MM:SS");
	if (strlen(datetime) != expected_len) return 1;
	char year_c[5], month_c[3], day_c[3], hour_c[3], minute_c[3], second_c[3];
	strncpy(year_c, datetime, 4);
	strncpy(month_c,datetime + DATETIME_MINUTE_OFFSET, 2);
	strncpy(day_c, datetime + DATETIME_DAY_OFFSET, 2);
	strncpy(hour_c, datetime + DATETIME_HOUR_OFFSET, 2);
	strncpy(minute_c,datetime + DATETIME_MINUTE_OFFSET, 2);
	strncpy(second_c, datetime + DATETIME_SECOND_OFFSET, 2);
	year_c[4] = '\0';
	month_c[2] = '\0';
	day_c[2] = '\0';
	hour_c[2] = '\0';
	minute_c[2] = '\0';
	second_c[2] = '\0';
	result_date->year = HandleStrtol(year_c);
	result_date->month = HandleStrtol(month_c);
	result_date->day = HandleStrtol(day_c);
	result_time->hour = HandleStrtol(hour_c);
	result_time->minute = HandleStrtol(minute_c);
	result_time->second = HandleStrtol(second_c);
	result_time->millisecond = 0;
	return 0;
}

// convert struct tm back to format like YYYY-MM-DD HH:MM:SS
int DatetimeStringFromDateAndTime(char *output_string, Date *input_date, Time *input_time) {
	sprintf(output_string, "%d", input_date->year);
	int err = CastIntTo2Chars(output_string + DATETIME_MONTH_OFFSET, input_date->month);
	if (err) return err;
	sprintf(output_string + DATETIME_DAY_OFFSET, "%d", input_date->day);
	sprintf(output_string + DATETIME_HOUR_OFFSET, "%d", input_time->hour);
	sprintf(output_string + DATETIME_MINUTE_OFFSET, "%d", input_time->minute);
	sprintf(output_string + DATETIME_SECOND_OFFSET, "%d", input_time->second);
	output_string[DATETIME_MONTH_OFFSET - 1] = '-';
	output_string[DATETIME_DAY_OFFSET - 1] = '-';
	output_string[DATETIME_HOUR_OFFSET - 1] = ' ';
	output_string[DATETIME_MINUTE_OFFSET - 1] = ':';
	output_string[DATETIME_SECOND_OFFSET - 1] = ':';
	output_string[DATETIME_STR_LEN] = '\0';
	return 0;
}

/*** Data structures ***/
#define DTCC_ID_COL "DISSEMINATION IDENTIFIER"
#define DTCC_START_COL "EFFECTIVE DATE"
#define DTCC_END_COL "EXPIRATION DATE"
#define DTCC_TRADE_TIME_COL "EXECUTION TIMESTAMP"
#define DTCC_FIXED_RATE_COL_1 "FIXED RATE-LEG 1"
#define DTCC_FIXED_RATE_COL_2 "FIXED RATE-LEG 2"
#define DTCC_NOTIONAL_COL_1 "NOTIONAL AMOUNT-LEG 1"
#define DTCC_NOTIONAL_COL_2 "NOTIONAL AMOUNT-LEG 2"
#define DTCC_ACTION_TYPE_COL "ACTION TYPE"
#define DTCC_TRANSACTION_TYPE_COL "TRANSACTION TYPE"
#define DTCC_BLOCK_TRADE_COL "BLOCK TRADE ELECTION INDICATOR"
#define DTCC_CLEARING_COL "CLEARED"
#define DTCC_VENUE_COL "EXECUTION VENUE TYPE"
#define DTCC_REF_RATE_COL_1 "UNDERLIER ID-LEG 1"
#define DTCC_REF_RATE_COL_2 "UNDERLIER ID-LEG 2"
#define DTCC_CURRENCY_COL_1 "NOTIONAL CURRENCY-LEG 1"
#define DTCC_CURRENCY_COL_2 "NOTIONAL CURRENCY-LEG 2"
#define DTCC_FIXED_PAYMENT_FREQUENCY_COL_1 "FIXED RATE PAYMENT FREQUENCY PERIOD-LEG 1"
#define DTCC_FLOATING_PAYMENT_FREQUENCY_COL_1 "FLOATING RATE PAYMENT FREQUENCY PERIOD-LEG 1"
#define DTCC_FIXED_PAYMENT_FREQUENCY_COL_2 "FIXED RATE PAYMENT FREQUENCY PERIOD-LEG 2"
#define DTCC_FLOATING_PAYMENT_FREQUENCY_COL_2 "FLOATING RATE PAYMENT FREQUENCY PERIOD-LEG 2"

void SwapAttributeLine(StringBuffer *buffer, const char *attr_name, char *attr_value) {
	char line[32];
	sprintf(line, "%s: %s \n", attr_name, attr_value);
	StringBufferAppendString(buffer, line);
}


void LogSwapDetails(Swap *swap) {
	char msg_buffer[256];
	StringBuffer buffer;
	StringBufferInit(&buffer);
	// ID
	sprintf(msg_buffer, "%ld", swap->id);
	SwapAttributeLine(&buffer, "ID", msg_buffer);
	// Notional
	if (swap->header.precision == Precision::FLOAT) {
		SwapScheduleF *schedule = (SwapScheduleF *)(swap->schedule);
		sprintf(msg_buffer, "%lf", schedule->notional);
	}
	else {
		SwapScheduleD *schedule = (SwapScheduleD *)(swap->schedule);
		sprintf(msg_buffer, "%ld", schedule->notional);

	}
	SwapAttributeLine(&buffer, "Notional", msg_buffer);
	// Start
	DateStringFromDate(msg_buffer, &swap->start_date);
	SwapAttributeLine(&buffer, "Start date", msg_buffer);
	// End
	DateStringFromDate(msg_buffer, &swap->end_date);
	SwapAttributeLine(&buffer, "End date", msg_buffer);
	// Trade time
	DatetimeStringFromDateAndTime(msg_buffer, &swap->trade_date, &swap->trade_time);
	SwapAttributeLine(&buffer, "Traded at", msg_buffer);
	// Fixed rate
	if (swap->header.precision == Precision::FLOAT) {
		SwapScheduleF *schedule = (SwapScheduleF *)(swap->schedule);
		sprintf(msg_buffer, "%f", schedule->fixed_rate);
	}
	else {
		SwapScheduleD *schedule = (SwapScheduleD *)(swap->schedule);
		sprintf(msg_buffer, "%d", schedule->fixed_rate);
	}
	SwapAttributeLine(&buffer, "Fixed rate", msg_buffer);
	printf(buffer.string);
	StringBufferClear(&buffer);
}

enum class SwapAttribute {
	ID,
	START_DATE,
	END_DATE,
	TRADE_TIME,
	FIXED_RATE,
	NOTIONAL,
	ACTION_TYPE,
	TRANSACTION_TYPE,
	IS_BLOCK_TRADE,
	VENUE,
	REF_RATE_IN_COL_1,
	REF_RATE_IN_COL_2,
	CURRENCY,
	PAY_FREQ_1,
	PAY_FREQ_2,
	FLOAT_PAY_FREQ,
	FIXED_PAY_FREQ,
	REF_RATE,
	PARSE_ERROR,
};

/*** CSV parsing ***/
// Read a line into a buffer
int ReadLine(FILE *file_handle, char *buffer, size_t max_buffer_size) {
	size_t buff_size = 0;
	char this_char = 0;
	while ((this_char = (char)getc(file_handle)) != '\n') {
		if (buff_size == max_buffer_size) return -1;
		else if (this_char == EOF) {
			return (int)buff_size;
		}
		buffer[buff_size++] = this_char;
	}
	buffer[buff_size + 1] = '\0';
	return (int)buff_size;
}

void StringBufferUppercase(StringBuffer *target_string) {
	char *this_char = target_string->string;
	for (size_t char_idx = 0; char_idx < target_string->length; ++char_idx) {
		*this_char++ = (char)toupper((int)*this_char);
	}
}

int StringBufferContains(StringBuffer *target_string, StringBuffer *match_string) {
	if (target_string->length < match_string->length) return 0;
	char *last_chunk_start = target_string->string + target_string->length - match_string->length;
	for (char *chunk_start = target_string->string; chunk_start < last_chunk_start; ++chunk_start) {
		if (!strncmp(chunk_start, match_string->string, match_string->length)) return 1;
	}
	return 0;
}

void StringUppercase(char *target_string) {
	size_t target_string_size = strlen(target_string);
	for (size_t char_idx = 0; char_idx < target_string_size; char_idx++) {
		target_string[char_idx] = (char)toupper((int)target_string[char_idx]);
	}
}

int StringContains(const char *target_string, const char *match_string) {
	size_t target_string_len = strlen(target_string);
	size_t match_string_len = strlen(match_string);
	if (match_string_len > target_string_len) {
		return 0;
	}
	for (char *chunk_start = const_cast<char *>(target_string);
		chunk_start < target_string + target_string_len - match_string_len; chunk_start++) {
		if (strncmp(chunk_start, match_string, match_string_len) == 0) {
			return 1;
		}
	}
	return 0;
}

SwapAttribute EvaluateColname(const char *col_name) {
	StringBuffer upper_case_col_name;
	StringBufferInit(&upper_case_col_name);
	StringBufferAppendString(&upper_case_col_name, col_name);
	StringBufferUppercase(&upper_case_col_name);
	SwapAttribute res;
	if (strcmp(upper_case_col_name.string, DTCC_ID_COL) == 0) {
		res = SwapAttribute::ID;
	}
	else if (strcmp(upper_case_col_name.string, DTCC_START_COL) == 0) {
		res = SwapAttribute::START_DATE;
	}
	else if (strcmp(upper_case_col_name.string, DTCC_END_COL) == 0) {
		res = SwapAttribute::END_DATE;
	}
	else if (strcmp(upper_case_col_name.string, DTCC_TRADE_TIME_COL) == 0) {
		res = SwapAttribute::TRADE_TIME;
	}
	else if ((strcmp(upper_case_col_name.string, DTCC_FIXED_RATE_COL_1) == 0) ||
		(strcmp(upper_case_col_name.string, DTCC_FIXED_RATE_COL_2) == 0)) {
		res = SwapAttribute::FIXED_RATE;
	}
	else if ((strcmp(upper_case_col_name.string, DTCC_NOTIONAL_COL_1) == 0) ||
		(strcmp(upper_case_col_name.string, DTCC_NOTIONAL_COL_2) == 0)) {
		res = SwapAttribute::NOTIONAL;
	}
	else if (strcmp(upper_case_col_name.string, DTCC_ACTION_TYPE_COL) == 0) {
		res = SwapAttribute::ACTION_TYPE;
	}
	else if (strcmp(upper_case_col_name.string, DTCC_TRANSACTION_TYPE_COL) == 0) {
		res = SwapAttribute::TRANSACTION_TYPE;
	}
	else if (strcmp(upper_case_col_name.string, DTCC_BLOCK_TRADE_COL) == 0) {
		res = SwapAttribute::IS_BLOCK_TRADE;
	}
	else if (strcmp(upper_case_col_name.string, DTCC_VENUE_COL) == 0) {
		res = SwapAttribute::VENUE;
	}
	else if (strcmp(upper_case_col_name.string, DTCC_REF_RATE_COL_1) == 0) {
		res = SwapAttribute::REF_RATE_IN_COL_1;
	}
	else if ((strcmp(upper_case_col_name.string, DTCC_REF_RATE_COL_2) == 0)) {
		res = SwapAttribute::REF_RATE_IN_COL_2;
	}
	else if ((strcmp(upper_case_col_name.string, DTCC_CURRENCY_COL_1) == 0) ||
		(strcmp(upper_case_col_name.string, DTCC_CURRENCY_COL_2) == 0)) {
		res = SwapAttribute::CURRENCY;
	}
	else if ((strcmp(upper_case_col_name.string, DTCC_FLOATING_PAYMENT_FREQUENCY_COL_1) == 0) ||
		(strcmp(upper_case_col_name.string, DTCC_FLOATING_PAYMENT_FREQUENCY_COL_2) == 0)) {
		res = SwapAttribute::FLOAT_PAY_FREQ;
	}
	else if ((strcmp(upper_case_col_name.string, DTCC_FIXED_PAYMENT_FREQUENCY_COL_1) == 0) ||
		(strcmp(upper_case_col_name.string, DTCC_FIXED_PAYMENT_FREQUENCY_COL_2) == 0)) {
		res = SwapAttribute::FIXED_PAY_FREQ;
	}
	else if ((strcmp(upper_case_col_name.string, DTCC_REF_RATE_COL_1) == 0) ||
		(strcmp(upper_case_col_name.string, DTCC_REF_RATE_COL_2) == 0)) {
		res = SwapAttribute::REF_RATE;
	}
	else {
		res = SwapAttribute::PARSE_ERROR;
	}
	free(upper_case_col_name.string);
	return res;
}


RefRate ParseRefRate(char *ref_rate) {
	enum RefRate res;
	StringUppercase(ref_rate);
	const char sofr[] = "SOFR";
	const char compound[] = "COMPOUND";
	const char term[] = "TERM";
	const char cpi[] = "CPI";
	const char libor[] = "LIBOR";
	if (StringContains(ref_rate, sofr) &&
		StringContains(ref_rate, compound)) {
		res = RefRate::USSOFR;
	}
	else if (StringContains(ref_rate, sofr) &&
		StringContains(ref_rate, term)) {
		res = RefRate::USSTERM;
	}
	else if (StringContains(ref_rate, cpi)) {
		res = RefRate::USCPI;
	}
	else if (StringContains(ref_rate, libor)) {
		res = RefRate::USLIBOR;
	}
	else {
		res = RefRate::REF_RATE_ERROR;
	}
	return res;
}

PayFreq ParsePayFreq(const char *input) {
	PayFreq res;
	if (strcmp(input, "1M") == 0) {
		res = PayFreq::M;
	}
	else if (strcmp(input, "3M") == 0) {
		res = PayFreq::Q;
	}
	else if (strcmp(input, "6M") == 0) {
		res = PayFreq::S;
	}
	else if (strcmp(input, "1Y") == 0) {
		res = PayFreq::A;
	}
	else {
		res = PayFreq::PAY_FREQ_ERROR;
	}
	return res;
}

// Reads the line of .csv into elem_array and returns the number of found
// elements
int ParseDTCCLine(
	char *elem_array,  // target pointer to update (pre-allocated memory of size
	// max_colname_len * max_n_cols)
	const char *input_line,  // line of .csv text
	int n_chars,             // size of input_line
	size_t max_colname_len   // size allocated for each title in elem_array
) {
	size_t elem_idx = 0;
	int max_buff_size = 64;
	char *buffer = (char *)malloc(max_buff_size * sizeof(char));
	if (!buffer) return 1;
	int buffer_size = 0;
	int char_idx = 0;
	char this_char = 0;
	int open_quote = 0;
	while (char_idx <= n_chars) {
		this_char = input_line[char_idx];
		// we assume special cases (containing "," in their value) are passed in
		// between " "
		if ((this_char == ',') && (open_quote == 0)) {
			buffer[buffer_size] = '\0';
			strcpy(&elem_array[elem_idx], buffer);
			elem_idx = elem_idx + max_colname_len;
			buffer_size = 0;
		}
		else if ((this_char == '\n') || (this_char == EOF)) {
			break;
		}
		else if (this_char == '"') {
			if (open_quote == 1) {
				open_quote = 0;
			}
			else {
				open_quote = 1;
			}
		}
		else {
			if (buffer_size + 1 == max_buff_size) {
				printf("Buffer contains %s\n", buffer);
				free(buffer);
				return -1;
			}
			buffer[buffer_size++] = this_char;
		}
		char_idx++;
	}
	free(buffer);
	return (int)elem_idx / (int)max_colname_len;
}

void SwapSetValue(Swap *swap, SwapAttribute attr_name, char *attr_value) {
	SwapHeader *header = (SwapHeader *)swap;
	switch (attr_name) {
	case SwapAttribute::ID:
		swap->id = HandleStrtol(attr_value);
		break;
	case SwapAttribute::START_DATE:
		ParseDate(attr_value, &swap->start_date);
		break;
	case SwapAttribute::END_DATE:
		ParseDate(attr_value, &swap->end_date);
		break;
	case SwapAttribute::TRADE_TIME:
		ParseDatetime(attr_value, &swap->trade_date, &swap->trade_time);
		break;
	case SwapAttribute::FIXED_RATE:
		if (header->precision == Precision::FLOAT) {
			SwapScheduleF *schedule = (SwapScheduleF *)(swap->schedule);
			if (strlen(attr_value) > 0) schedule->fixed_rate = HandleStrtof(attr_value);
		}
		else {
			SwapScheduleD *schedule = (SwapScheduleD *)(swap->schedule);
			if (strlen(attr_value) > 0) schedule->fixed_rate = HandleStrtod(attr_value);
		}
		break;
	case SwapAttribute::NOTIONAL:
		if (header->precision == Precision::FLOAT) {
			SwapScheduleF *schedule = (SwapScheduleF *)(swap->schedule);
			schedule->notional = HandleStrtof(attr_value);
		}
		else {
			SwapScheduleD *schedule = (SwapScheduleD *)(swap->schedule);
			schedule->notional = HandleStrtod(attr_value);
		}
		break;
	case SwapAttribute::ACTION_TYPE:
		if (strcmp(attr_value, "NEWT") == 0) {
			swap->action_type = Swap::NEW;
		}
		else if (strcmp(attr_value, "CORR") == 0) {
			swap->action_type = Swap::CORRECT;
		}
		else if (strcmp(attr_value, "MODI") == 0) {
			swap->action_type = Swap::MODIFY;
		}
		else if (strcmp(attr_value, "CANCEL") == 0) {
			swap->action_type = Swap::CANCEL;
		}
		else {
			swap->action_type = Swap::ACTION_ERROR;
		}
		break;
	case SwapAttribute::TRANSACTION_TYPE:
		if (strcmp(attr_value, "Trade") == 0) {
			swap->transaction_type = Swap::TRADE;
		}
		else if (strcmp(attr_value, "Amendment") == 0) {
			swap->transaction_type = Swap::AMENDMENT;
		}
		else if (strcmp(attr_value, "Termination") == 0) {
			swap->transaction_type = Swap::TERMINATION;
		}
		else {
			swap->transaction_type = Swap::TRANSACTION_ERROR;
		}
		break;
	case SwapAttribute::IS_BLOCK_TRADE:
		if (strcmp(attr_value, "Y") == 0) {
			swap->is_block_trade = Swap::Y;
		}
		else if (strcmp(attr_value, "N") == 0) {
			swap->is_block_trade = Swap::N;
		}
		else {
			swap->is_block_trade = Swap::BLOCK_TRADE_ERROR;
		}
		break;
	case SwapAttribute::VENUE:
		if (strcmp(attr_value, "OFF") == 0) {
			swap->venue = Swap::OFF;
		}
		else if (strcmp(attr_value, "ON") == 0) {
			swap->venue = Swap::ON;
		}
		else {
			swap->venue = Swap::VENUE_ERROR;
		}
		break;
	case SwapAttribute::CURRENCY:
		if (strcmp(attr_value, "USD") == 0) {
			swap->currency = Currency::USD;
		}
		else if (strcmp(attr_value, "EUR") == 0) {
			swap->currency = Currency::EUR;
		}
		else {
			swap->currency = Currency::CURRENCY_ERROR;
		}
		break;
	case SwapAttribute::FLOAT_PAY_FREQ:
		swap->floating_pay_frequency = ParsePayFreq(attr_value);
		break;
	case SwapAttribute::FIXED_PAY_FREQ:
		swap->fixed_pay_frequency = ParsePayFreq(attr_value);
		break;
	case SwapAttribute::REF_RATE:
		swap->ref_rate = ParseRefRate(attr_value);
		break;
	case SwapAttribute::PARSE_ERROR:
		break;
	default:
		break;
	}
}

// Reads data from one line of the csv file into a Swap structure
int SwapFromCSVLine(
	Swap *result_swap, // output
	const char *input_line,  // input line of .csv text
	size_t line_size,           // size in input_line
	const char *data_cols,   // column names, assumed to have
	// been Parsed from ParseLine
	size_t max_colname_len  // size of each column name element
) {
	StringBuffer buffer, pay_frequency_1, pay_frequency_2;
	StringBufferInit(&buffer);
	bool has_pay_frequency_1 = false;
	bool has_pay_frequency_2 = false;

	int char_idx = 0;
	size_t col_idx = 0;
	int col_1_is_float = 0;
	// Use this as a flag for parsing entries enclosed with " "
	int open_quotes = 0;
	while (char_idx < line_size) {
		const char this_char = input_line[char_idx];
		// if we reach a ',' and we have not found open quotes, we check out column
		// name and Parse our buffer into the corresponding attribute of swap
		if ((this_char == ',') && (!open_quotes)) {
			if (buffer.length > 0) {
				enum SwapAttribute attr_name = EvaluateColname(data_cols + col_idx);
				if (attr_name == SwapAttribute::REF_RATE_IN_COL_1) {
					col_1_is_float = 1;
					result_swap->ref_rate = ParseRefRate(buffer.string);
				}
				else if (attr_name == SwapAttribute::REF_RATE_IN_COL_2) {
					col_1_is_float = 0;
					result_swap->ref_rate = ParseRefRate(buffer.string);
				}
				else if (attr_name == SwapAttribute::PAY_FREQ_1) {
					StringBufferInit(&pay_frequency_1);
					has_pay_frequency_1 = true;
					StringBufferAppendStringBuffer(&pay_frequency_1, &buffer);
				}
				else if (attr_name == SwapAttribute::PAY_FREQ_2) {
					StringBufferInit(&pay_frequency_2);
					has_pay_frequency_2 = true;
					StringBufferAppendStringBuffer(&pay_frequency_2, &buffer);
				}
				else {
					SwapSetValue(result_swap, attr_name, buffer.string);
				}
			}
			StringBufferReset(&buffer);
			col_idx = col_idx + max_colname_len;
		}
		else if (this_char == '"') {
			// if we run into a quote, we switch the flag and continue parsing,
			// ignoring the quote
			if (!open_quotes) {
				open_quotes = 1;
			}
			else {
				open_quotes = 0;
			}
		}
		else {
			StringBufferAppendChar(&buffer, &this_char);
		}
		char_idx++;
	}
	if (col_1_is_float) {
		result_swap->floating_pay_frequency = has_pay_frequency_1 ? ParsePayFreq(pay_frequency_1.string) : PayFreq::PAY_FREQ_ERROR;
		result_swap->fixed_pay_frequency = has_pay_frequency_2 ? ParsePayFreq(pay_frequency_2.string) : PayFreq::PAY_FREQ_ERROR;
	}
	else {
		result_swap->floating_pay_frequency = has_pay_frequency_1 ? ParsePayFreq(pay_frequency_2.string) : PayFreq::PAY_FREQ_ERROR;
		result_swap->fixed_pay_frequency = has_pay_frequency_2 ? ParsePayFreq(pay_frequency_1.string) : PayFreq::PAY_FREQ_ERROR;
	}
	free(buffer.string);
	if (has_pay_frequency_1) free(pay_frequency_1.string);
	if (has_pay_frequency_2) free(pay_frequency_2.string);
	return 0;
}

int ListStringToSwap(Swap *swap, const char *input, size_t input_size)
{
	// Expect a list to be passed in like "Colname:Value;"
	const size_t max_buff_size = 64;
	char buff[max_buff_size], attr_name[max_buff_size], attr_value[max_buff_size];
	size_t buff_idx = 0, attr_name_idx = 0, attr_value_idx = 0;
	size_t current_input_idx = 0;
	char this_char = 0;
	char reading_attr_name = 1;
	while (current_input_idx < input_size)
	{
		buff_idx = 0, attr_name_idx = 0, attr_value_idx = 0;
		while (((this_char = input[current_input_idx++]) != ';') &&
			(current_input_idx < input_size))
		{
			buff[buff_idx++] = this_char;
			if (buff_idx > max_buff_size) return 1;
			if ((reading_attr_name == 1) && (this_char != ':'))
			{
				attr_name[attr_name_idx++] = this_char;
			}
			else if (reading_attr_name == 1)
			{
				attr_name[attr_name_idx] = '\0';
				reading_attr_name = 0;
			}
			else
			{
				attr_value[attr_value_idx++] = this_char;
			}
		}
		attr_value[attr_value_idx] = '\0';
		buff[buff_idx] = '\0';
		enum SwapAttribute Parsed_attr_name = EvaluateColname(attr_name);
		SwapSetValue(swap, Parsed_attr_name, attr_value);
	}
	return 0;
}

int SwapFromInputLine(Swap *swap, const char *input_line)
{
	size_t begin_idx = 0;
	size_t separator_idx = 0;
	size_t end_idx = 0;
	char attribute_buffer[64];
	char value_buffer[64];
	SwapAttribute this_attribute;
	while (end_idx++ < strlen(input_line))
	{
		if (input_line[end_idx] == ';')
		{
			separator_idx = begin_idx;
			while (input_line[separator_idx] != ':')
			{
				separator_idx++;
			}
			strncpy(attribute_buffer, input_line + begin_idx,
				separator_idx - begin_idx);
			strncpy(value_buffer, input_line + separator_idx,
				end_idx - separator_idx);
			begin_idx = end_idx + 1;
			this_attribute = EvaluateColname(attribute_buffer);
			SwapSetValue(swap, this_attribute, value_buffer);
		}
	}
	return 0;
}

int getCFTCDownloadUrl(char *output_url_buffer, Date *input_date) {
	char date_string[date_string_size + 1];
	DateToString(date_string, input_date);
	const char url_prefix[] = "https://kgc0418-tdw-data-0.s3.amazonaws.com/cftc/eod/CFTC_CUMULATIVE_RATES_";
	strncpy(output_url_buffer, url_prefix, strlen(url_prefix));
	strncpy(output_url_buffer + strlen(url_prefix), date_string, date_string_size);
	snprintf(output_url_buffer + strlen(url_prefix) + date_string_size, 5, ".zip");
	return 0;
}


struct ColumnNames {
	size_t max_column_name_length;
	size_t number_of_columns;
	char *contents;

	~ColumnNames() { if (contents) free(contents); }
};

int ParseLine(
	char *elem_array,  // target pointer to update (pre-allocated memory of size
	// max_colname_len * max_n_cols)
	const char *input_line,  // line of .csv text
	int n_chars,             // size of input_line
	size_t max_colname_len   // size allocated for each title in elem_array
) {
	size_t elem_idx = 0;
	int max_buffer_size = 64;
	char *buffer = (char *)malloc(max_buffer_size * sizeof(char));
	if (!buffer) return -1;
	int buffer_size = 0;
	int char_idx = 0;
	char this_char = 0;
	int open_quote = 0;
	while (char_idx <= n_chars) {
		this_char = input_line[char_idx];
		// we assume special cases (containing "," in their value) are passed in
		// between " "
		if ((this_char == ',') && (open_quote == 0)) {
			buffer[buffer_size] = '\0';
			memcpy(elem_array + elem_idx, buffer, buffer_size + 1);
			elem_idx = elem_idx + max_colname_len;
			buffer_size = 0;
		}
		else if ((this_char == '\n') || (this_char == EOF)) {
			break;
		}
		else if (this_char == '"') {
			if (open_quote == 1) {
				open_quote = 0;
			}
			else {
				open_quote = 1;
			}
		}
		else {
			if (buffer_size + 1 == max_buffer_size) {
				printf("Buffer contains %s\n", buffer);
				return -1;
			}
			buffer[buffer_size++] = this_char;
		}
		char_idx++;
	}
	free(buffer);
	return (int)(elem_idx / max_colname_len);
}



int GetDTCCData(Swap *destination_swap_list, size_t destination_swap_list_size, Date *date) {
	CURL *curl;
	FILE *file_handle;
	CURLcode res;
	int n_loaded_swaps = 0;
	char url[128];
	getCFTCDownloadUrl(url, date);

	StringBuffer outfile_name;
	StringBufferInit(&outfile_name);
	char date_string[date_string_size + 1];
	DateToString(date_string, date);
	StringBufferAppendStringWithSize(&outfile_name, date_string, date_string_size);
	StringBufferAppendString(&outfile_name, "rates.zip");

	curl_version_info_data *vinfo = curl_version_info(CURLVERSION_NOW);

	if (vinfo->features & CURL_VERSION_SSL) {
		printf("CURL: SSL enabled\n");
	}
	else {
		printf("CURL: SSL not enabled\n");
	}
	curl = curl_easy_init();
	if (curl) {

		/* Setup the https:// verification options. Note we   */
		/* do this on all requests as there may be a redirect */
		/* from http to https and we still want to verify     */
		curl_easy_setopt(curl, CURLOPT_URL, url);
		curl_easy_setopt(curl, CURLOPT_CAINFO, "./ca-bundle.crt");
		curl_easy_setopt(curl, CURLOPT_SSL_VERIFYPEER, false);
		curl_easy_setopt(curl, CURLOPT_SSL_VERIFYHOST, false);

		curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, fwrite);

		file_handle = fopen(outfile_name.string, "wb");
		curl_easy_setopt(curl, CURLOPT_WRITEDATA, file_handle);
		res = curl_easy_perform(curl);
		curl_easy_cleanup(curl);
		/*
		int error_code = fclose(file_handle);
		if (error_code != 0) {
			//free(outfile_name.string);
			return error_code;
		}
		*/
		// unzip the file
		StringBuffer system_command;
		StringBufferInit(&system_command);
		StringBufferAppendString(&system_command, "tar -xf ");
		StringBufferAppendStringBuffer(&system_command, &outfile_name);
		system(system_command.string);
		free(system_command.string);
		free(outfile_name.string);

		StringBuffer loadfile_name;
		StringBufferInit(&loadfile_name);
		StringBufferAppendString(&loadfile_name, "CFTC_CUMULATIVE_RATES_");
		StringBufferAppendString(&loadfile_name, date_string);
		StringBufferAppendString(&loadfile_name, ".csv");

		file_handle = fopen(loadfile_name.string, "rb");
		free(loadfile_name.string);
		size_t max_line_buffer_size = 4096;
		char *line_buffer = (char *)malloc(max_line_buffer_size * sizeof(char));
		if (file_handle) {
			// get column names
			size_t max_column_name_length = 64;
			size_t max_number_of_columns = 128;
			char *column_names = (char *)malloc(max_column_name_length * max_number_of_columns * sizeof(char));
			int line_size = ReadLine(file_handle, line_buffer, max_line_buffer_size);
			ParseLine(column_names, line_buffer, line_size, max_column_name_length);
			// read swaps into array
			destination_swap_list_size = destination_swap_list_size > 2048 ? 2048 : destination_swap_list_size;
			Swap *current_output_swap = destination_swap_list;
			for (int i = 0; i < destination_swap_list_size; ++i)
			{
				line_size = ReadLine(file_handle, line_buffer, max_line_buffer_size);
				if (line_size == 0)
				{
					break;
				}
				SwapFromCSVLine(current_output_swap++, line_buffer, line_size, column_names, max_column_name_length);
				++n_loaded_swaps;
			}
			printf("%d swaps loaded\n", n_loaded_swaps);
			free(line_buffer);
			free(column_names);
			fclose(file_handle);
		}
		else {
			return 0;
		}
	}
	else {
		return 0;
	}
	return n_loaded_swaps;
}

/*
StartupContext LoadSwapsFromFile(const char *filename, int max_n_cols,
	int chunk_size, int max_colname_len,
	int max_n_loaded_swaps)
{
	Colnames colnames = { 0 };
	SwapList swap_list = { 0 };
	swap_list.contents = malloc(max_n_loaded_swaps * sizeof(Swap));
	colnames.contents = malloc(max_n_cols * max_colname_len);
	colnames.max_colname_len = max_colname_len;
	char line_buffer[chunk_size];
	FILE *handler = fopen(filename, "r");
	if (handler)
	{
		// get column names
		int line_size = ReadLine(handler, line_buffer, chunk_size);
		colnames.n_colnames =
			ParseLine(colnames.contents, line_buffer, line_size, max_colname_len);
		// read swaps into array
		int n_loaded_swaps = 0;
		for (int i = 0; i < max_n_loaded_swaps; i++)
		{
			line_size = ReadLine(handler, line_buffer, chunk_size);
			if (line_size == 0)
			{
				break;
			}
			n_loaded_swaps = i + 1;
			swap_list.contents[i] = SwapFromCSVLine(
				line_buffer, line_size, colnames.contents, max_colname_len);
		}
		printf("%d swaps loaded\n", n_loaded_swaps);
		swap_list.size = n_loaded_swaps;
	}
	fclose(handler);
	StartupContext startup_context;
	startup_context.colnames = colnames;
	startup_context.swap_list = swap_list;
	return startup_context;
}

*/
