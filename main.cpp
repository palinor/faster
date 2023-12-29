#define CURL_STATICLIB
#include <errno.h>
#include <stdio.h>
#include <curl/curl.h>
#include <curl/easy.h>
#include <string.h>
#include <stdlib.h>

/*** Date parsing ***/
#define DATE_STR_LEN 11      // strlen("2022-09-10") + '\0'
#define DATETIME_STR_LEN 21  // strlen("2022-09-10T20:15:56") + '\0'
#define DATETIME_MONTH_OFFSET 5
#define DATETIME_DAY_OFFSET 8
#define DATETIME_HOUR_OFFSET 11
#define DATETIME_MINUTE_OFFSET 14
#define DATETIME_SECOND_OFFSET 17

size_t write_data(void *ptr, size_t size, size_t nmemb, FILE *stream) {
    size_t written;
    written = fwrite(ptr, size, nmemb, stream);
    return written;
}

struct Date {
    int year;
    int month;
    int day;
};

int date_to_string(char *outputString, size_t *outputSize, Date *inputDate) {
    *outputSize = 11;
    snprintf(outputString, 5, "%i", inputDate->year);
    outputString[4] = '_';
    if (inputDate->month < 10) {
        snprintf(outputString + 5, 2, "0");
        snprintf(outputString + 6, 2, "%i", inputDate->month);
    }
    else {
        snprintf(outputString + 5, 3, "%i", inputDate->month);
    }
    outputString[7] = '_';
    if (inputDate->day < 10) {
        snprintf(outputString + 8, 2, "0");
        snprintf(outputString + 9, 2, "%i", inputDate->day);
    }
    else {
        snprintf(outputString + 8, 3, "%i", inputDate->day);
    }
    outputString[10] = '\0';
    return 0;
}

int get_url(char *outputUrl, Date *inputDate) {
    char *dateString = (char *)malloc(16);
    size_t stringSize = 16;
    date_to_string(dateString, &stringSize, inputDate);
    const char urlPrefix[] = "https://kgc0418-tdw-data-0.s3.amazonaws.com/cftc/eod/CFTC_CUMULATIVE_RATES_";
    strncpy(outputUrl, urlPrefix, strlen(urlPrefix));
    strncpy(outputUrl + strlen(urlPrefix), dateString, stringSize);
    snprintf(outputUrl + strlen(urlPrefix) + stringSize - 1, 5, ".zip");
    free(dateString);
    return 0;
}

int get_dtcc_data(void) {
    CURL *curl;
    FILE *fp;
    CURLcode res;

    char *url = (char *)malloc(128);
    Date today;
    today.year = 2023;
    today.month = 12;
    today.day = 1;
    get_url(url, &today);

    const char outfilename[FILENAME_MAX] = "./json.zip";

    curl_version_info_data * vinfo = curl_version_info(CURLVERSION_NOW);

    if(vinfo->features & CURL_VERSION_SSL){
        printf("CURL: SSL enabled\n");
    }else{
        printf("CURL: SSL not enabled\n");
    }

    curl = curl_easy_init();
    if (curl) {
        fp = fopen(outfilename,"wb");        

        /* Setup the https:// verification options. Note we   */
        /* do this on all requests as there may be a redirect */
        /* from http to https and we still want to verify     */
        curl_easy_setopt(curl, CURLOPT_URL, url);
        curl_easy_setopt(curl, CURLOPT_CAINFO, "./ca-bundle.crt");
        curl_easy_setopt(curl, CURLOPT_SSL_VERIFYPEER, false);
        curl_easy_setopt(curl, CURLOPT_SSL_VERIFYHOST, false);

        curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, write_data);

        curl_easy_setopt(curl, CURLOPT_WRITEDATA, fp);
        res = curl_easy_perform(curl);
        curl_easy_cleanup(curl);

        int i=fclose(fp);
        if( i==0)
        system("unzip -j json.zip");
    }
    return 0;
}


void StripCommas(const char *input, char *output, size_t input_size) {
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
  StripCommas(input, buffer, input_len + 1);
  char *endptr = NULL;
  long res = strtol(buffer, &endptr, 10);
  if (res == 0) {
    if (errno != 0) {
      char err_msg[64];
      sprintf(err_msg, "Encountered error in strtol, errno : %d, input: %s",
              errno, input);
      log(err_msg);
    } else if (buffer == endptr) {
      char err_msg[64];
      sprintf(err_msg,
              "Encountered error in strtol, no characters read from input: %s",
              input);
      log(err_msg);
    }
  }
  return res;
}

float HandleStrtof(const char *input) {
  char buffer[32];
  size_t input_len = strlen(input);
  StripCommas(input, buffer, input_len + 1);
  char *endptr = NULL;
  float res = strtof(buffer, &endptr);
  if (res == 0) {
    if (errno != 0) {
      char err_msg[64];
      sprintf(err_msg, "Encountered error in strtof, errno : %d, input: %s",
              errno, input);
      log(err_msg);
    } else if (buffer == endptr) {
      char err_msg[64];
      sprintf(err_msg,
              "Encountered error in strtof, no characters read from input: %s",
              input);
      log(err_msg);
    }
  }
  return res;
}

// assume date is in format YYYY-MM-DD
void ParseDate(const char *date, int *year, int *month, int *day) {
  size_t expected_len = strlen("YYYY-MM-DD");
  if (strlen(date) != expected_len) Die("Unsupported date format - wrong size");
  char year_c[5], month_c[3], day_c[3];
  strncpy(year_c, date, 4);
  strncpy(month_c, date + DATETIME_MONTH_OFFSET, 2);
  strncpy(day_c, date + DATETIME_DAY_OFFSET, 2);
  year_c[4] = '\0';
  month_c[2] = '\0';
  day_c[2] = '\0';
  *year = HandleStrtol(year_c) - 1900;
  *month = HandleStrtol(month_c);
  *day = HandleStrtol(day_c);
  return;
}

// handle how case where we want 1 to be cast to the string "01"
void CastIntTo2Chars(char *output, int input) {
  if ((input < 0) || (input > 99)) {
    Die("Invalid number passed to CastIntTo2Chars");
  } else if (input < 10) {
    output[0] = '0';
    sprintf(output + 1, "%d", input);
  } else {
    sprintf(output, "%d", input);
  }
}

// take struct tm and turn it back into a string like YYYY-MM-DD
void DateFromTm(char *output_string, const struct tm input_tm) {
  sprintf(output_string, "%d", input_tm.tm_year + 1900);
  CastIntTo2Chars(output_string + 5, input_tm.tm_mon);
  CastIntTo2Chars(output_string + 8, input_tm.tm_mday);
  output_string[DATETIME_MONTH_OFFSET - 1] = '-';
  output_string[DATETIME_DAY_OFFSET - 1] = '-';
  output_string[DATE_STR_LEN - 1] = '\0';
  return;
}

// assume datetime is in format YYYY-MM-DDTHH:MM:SS
void ParseDatetime(const char *datetime, int *year, int *month, int *day,
                   int *hour, int *minute, int *second) {
  size_t expected_len = strlen("YYYY-MM-DDTHH:MM:SS");
  if (strlen(datetime) != expected_len)
    Die("Unsupported date format - wrong size");
  char year_c[5], month_c[3], day_c[3], hour_c[3], minute_c[3], second_c[3];
  strncpy(year_c, datetime, 4);
  strncpy(month_c, datetime + DATETIME_MINUTE_OFFSET, 2);
  strncpy(day_c, datetime + DATETIME_DAY_OFFSET, 2);
  strncpy(hour_c, datetime + DATETIME_HOUR_OFFSET, 2);
  strncpy(minute_c, datetime + DATETIME_MINUTE_OFFSET, 2);
  strncpy(second_c, datetime + DATETIME_SECOND_OFFSET, 2);
  year_c[4] = '\0';
  month_c[2] = '\0';
  day_c[2] = '\0';
  hour_c[2] = '\0';
  minute_c[2] = '\0';
  second_c[2] = '\0';
  *year = HandleStrtol(year_c) - 1900;
  *month = HandleStrtol(month_c);
  *day = HandleStrtol(day_c);
  *hour = HandleStrtol(hour_c);
  *minute = HandleStrtol(minute_c);
  *second = HandleStrtol(second_c);
  return;
}

// convert struct tm back to format like YYYY-MM-DD HH:MM:SS
void DatetimeFromTm(char *output_string, const struct tm input_tm) {
  sprintf(output_string, "%d", input_tm.tm_year);
  CastIntTo2Chars(output_string + DATETIME_MONTH_OFFSET, input_tm.tm_mon);
  sprintf(output_string + DATETIME_DAY_OFFSET, "%d", input_tm.tm_mday);
  sprintf(output_string + DATETIME_HOUR_OFFSET, "%d", input_tm.tm_hour);
  sprintf(output_string + DATETIME_MINUTE_OFFSET, "%d", input_tm.tm_min);
  sprintf(output_string + DATETIME_SECOND_OFFSET, "%d", input_tm.tm_sec);
  output_string[DATETIME_MONTH_OFFSET - 1] = '-';
  output_string[DATETIME_DAY_OFFSET - 1] = '-';
  output_string[DATETIME_HOUR_OFFSET - 1] = ' ';
  output_string[DATETIME_MINUTE_OFFSET - 1] = ':';
  output_string[DATETIME_SECOND_OFFSET - 1] = ':';
  output_string[DATETIME_STR_LEN] = '\0';
  return;
}

/*** Data structures ***/
#define ID_COL "Dissemination ID"
#define START_COL "Effective Date"
#define END_COL "Expiration Date"
#define TRADE_TIME_COL "Execution Timestamp"
#define FIXED_RATE_COL_1 "Fixed Rate 1"
#define FIXED_RATE_COL_2 "Fixed Rate 2"
#define NOTIONAL_COL "Notional Amount 1"
#define ACTION_TYPE_COL "Action"
#define TRANSACTION_TYPE_COL "Transaction Type"
#define BLOCK_TRADE_COL "Block Trade Election Indicator"
#define VENUE_COL "Execution Venue Type"
#define REF_RATE_COL_1 "Leg 1 - Floating Rate Index"
#define REF_RATE_COL_2 "Leg 2 - Floating Rate Index"
#define CURRENCY_COL "Notional Currency 1"
#define PAY_FREQ_COL_1 "Payment Frequency Period 1"
#define PAY_FREQ_COL_2 "Payment Frequency Period 2"
#define PAY_FREQ_COL "Payment Frequency"
#define FLOAT_PAY_FREQ_COL "Floating Payment Frequency"
#define FIXED_PAY_FREQ_COL "Fixed Payment Frequency"
#define REF_RATE_COL "Ref Rate"

typedef enum RefRate {
  REF_RATE_ERROR,
  USSOFR,
  USSTERM,
  USLIBOR,
  USCPI
} RefRate;
typedef enum PayFreq { PAY_FREQ_ERROR, M, Q, S, A } PayFreq;
typedef enum Currency { CURRENCY_ERROR, USD, EUR } Currency;

typedef struct Swap {
  long id;
  struct tm start_date;
  struct tm end_date;
  struct tm trade_time;
  float fixed_rate;
  float notional;
  RefRate ref_rate;
  PayFreq fixed_pay_freq;
  PayFreq float_pay_freq;
  Currency currency;
  enum { ACTION_ERROR, NEW, CANCEL, CORRECT } action_type;
  enum { TRANSACTION_ERROR, TRADE, AMENDMENT, TERMINATION } transaction_type;
  enum { BLOCK_TRADE_ERROR, Y, N } is_block_trade;
  enum { VENUE_ERROR, ON, OFF } venue;
} Swap;

void SwapAttributeLine(StringBuffer *buff, char *attr_name, char *attr_value) {
  char line[32];
  sprintf(line, "%s: %s \n", attr_name, attr_value);
  StringAppend(buff, line);
}

void LogSwapDetails(Swap *swap) {
  char msg_buffer[256];
  StringBuffer buff;
  StringInit(&buff);
  // ID
  sprintf(msg_buffer, "%ld", swap->id);
  SwapAttributeLine(&buff, "ID", msg_buffer);
  // Notional
  sprintf(msg_buffer, "%lf", swap->notional);
  SwapAttributeLine(&buff, "Notional", msg_buffer);
  // Start
  DateFromTm(msg_buffer, swap->start_date);
  SwapAttributeLine(&buff, "Start date", msg_buffer);
  // End
  DateFromTm(msg_buffer, swap->end_date);
  SwapAttributeLine(&buff, "End date", msg_buffer);
  // Trade time
  DatetimeFromTm(msg_buffer, swap->trade_time);
  SwapAttributeLine(&buff, "Traded at", msg_buffer);
  // Fixed rate
  sprintf(msg_buffer, "%f", swap->fixed_rate);
  SwapAttributeLine(&buff, "Fixed rate", msg_buffer);
  log(buff.string);
  StringClear(&buff);
}

typedef enum AttrToParse {
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
} AttrToParse;

/*** CSV parsing ***/
// Read a line into a buffer
int ReadLine(FILE *handler, char *buffer, int max_size) {
  int buff_size = 0;
  int this_char = 0;
  while ((this_char = getc(handler)) != '\n') {
    if (buff_size == max_size)
      Die("ReadLine - buffer too small");
    else if (this_char == EOF) {
      return buff_size;
    }
    buffer[buff_size++] = this_char;
  }
  buffer[buff_size + 1] = '\0';
  return buff_size;
}

enum AttrToParse EvaluateColname(const char *col_name) {
  enum AttrToParse res;
  if (strcmp(col_name, ID_COL) == 0) {
    res = ID;
  } else if (strcmp(col_name, START_COL) == 0) {
    res = START_DATE;
  } else if (strcmp(col_name, END_COL) == 0) {
    res = END_DATE;
  } else if (strcmp(col_name, TRADE_TIME_COL) == 0) {
    res = TRADE_TIME;
  } else if ((strcmp(col_name, FIXED_RATE_COL_1) == 0) |
             (strcmp(col_name, FIXED_RATE_COL_2) == 0)) {
    res = FIXED_RATE;
  } else if (strcmp(col_name, NOTIONAL_COL) == 0) {
    res = NOTIONAL;
  } else if (strcmp(col_name, ACTION_TYPE_COL) == 0) {
    res = ACTION_TYPE;
  } else if (strcmp(col_name, TRANSACTION_TYPE_COL) == 0) {
    res = TRANSACTION_TYPE;
  } else if (strcmp(col_name, BLOCK_TRADE_COL) == 0) {
    res = IS_BLOCK_TRADE;
  } else if (strcmp(col_name, VENUE_COL) == 0) {
    res = VENUE;
  } else if (strcmp(col_name, REF_RATE_COL_1) == 0) {
    res = REF_RATE_IN_COL_1;
  } else if ((strcmp(col_name, REF_RATE_COL_2) == 0)) {
    res = REF_RATE_IN_COL_2;
  } else if (strcmp(col_name, CURRENCY_COL) == 0) {
    res = CURRENCY;
  } else if (strcmp(col_name, PAY_FREQ_COL_1) == 0) {
    res = PAY_FREQ_1;
  } else if (strcmp(col_name, PAY_FREQ_COL_2) == 0) {
    res = PAY_FREQ_2;
  } else if (strcmp(col_name, FLOAT_PAY_FREQ_COL) == 0) {
    res = FLOAT_PAY_FREQ;
  } else if (strcmp(col_name, FIXED_PAY_FREQ_COL) == 0) {
    res = FIXED_PAY_FREQ;
  } else if (strcmp(col_name, REF_RATE_COL) == 0) {
    res = REF_RATE;
  } else {
    res = PARSE_ERROR;
  }
  return res;
}

void StringUpper(char *target_string) {
  size_t target_string_size = strlen(target_string);
  for (size_t char_idx = 0; char_idx < target_string_size; char_idx++) {
    target_string[char_idx] = toupper(target_string[char_idx]);
  }
}

int StringContains(const char *target_string, const char *match_string) {
  size_t target_string_len = strlen(target_string);
  size_t match_string_len = strlen(match_string);
  if (match_string_len > target_string_len) {
    return 0;
  }
  char chunk[match_string_len];
  for (size_t chunk_start = 0;
       chunk_start < target_string_len - match_string_len; chunk_start++) {
    strncpy(chunk, target_string + chunk_start, match_string_len);
    if (strcmp(chunk, match_string) == 0) {
      return 1;
    }
  }
  return 0;
}

enum RefRate ParseRefRate(char *ref_rate) {
  enum RefRate res;
  StringUpper(ref_rate);
  if (StringContains(ref_rate, "SOFR") &&
      StringContains(ref_rate, "COMPOUND")) {
    res = USSOFR;
  } else if (StringContains(ref_rate, "SOFR") &&
             StringContains(ref_rate, "TERM")) {
    res = USSTERM;
  } else if (StringContains(ref_rate, "CPI")) {
    res = USCPI;
  } else if (StringContains(ref_rate, "LIBOR")) {
    res = USLIBOR;
  } else {
    res = REF_RATE_ERROR;
  }
  return res;
}

enum PayFreq ParsePayFreq(const char *input) {
  enum PayFreq res;
  if (strcmp(input, "1M") == 0) {
    res = M;
  } else if (strcmp(input, "3M") == 0) {
    res = Q;
  } else if (strcmp(input, "6M") == 0) {
    res = S;
  } else if (strcmp(input, "1Y") == 0) {
    res = A;
  } else {
    res = PAY_FREQ_ERROR;
  }
  return res;
}

// Reads the line of .csv into elem_array and returns the number of found
// elements
int ParseLine(
    char *elem_array,  // target pointer to update (pre-allocated memory of size
                       // max_colname_len * max_n_cols)
    const char *input_line,  // line of .csv text
    int n_chars,             // size of input_line
    size_t max_colname_len   // size allocated for each title in elem_array
) {
  int elem_idx = 0;
  int max_buff_size = 64;
  char buffer[max_buff_size];
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
    } else if ((this_char == '\n') | (this_char == EOF)) {
      break;
    } else if (this_char == '"') {
      if (open_quote == 1) {
        open_quote = 0;
      } else {
        open_quote = 1;
      }
    } else {
      if (buffer_size + 1 == max_buff_size) {
        printf("Buffer contains %s\n", buffer);
        Die("Buffer has reached max size");
      }
      buffer[buffer_size++] = this_char;
    }
    char_idx++;
  }
  return elem_idx / max_colname_len;
}

void AssignSwapValue(Swap *swap_p, enum AttrToParse attr_name,
                     char *attr_value) {
  switch (attr_name) {
    case ID:
      swap_p->id = HandleStrtol(attr_value);
      break;
    case START_DATE:
      ParseDate(attr_value, &(swap_p->start_date.tm_year),
                &(swap_p->start_date.tm_mon), &(swap_p->start_date.tm_mday));
      break;
    case END_DATE:
      ParseDate(attr_value, &(swap_p->end_date.tm_year),
                &(swap_p->end_date.tm_mon), &(swap_p->end_date.tm_mday));
      break;
    case TRADE_TIME:
      ParseDatetime(attr_value, &(swap_p->trade_time.tm_year),
                    &(swap_p->trade_time.tm_mon), &(swap_p->trade_time.tm_mday),
                    &(swap_p->trade_time.tm_hour), &(swap_p->trade_time.tm_min),
                    &(swap_p->trade_time.tm_sec));
      break;
    case FIXED_RATE:
      if (strlen(attr_value) > 0) swap_p->fixed_rate = HandleStrtof(attr_value);
      break;
    case NOTIONAL:
      swap_p->notional = HandleStrtof(attr_value);
      break;
    case ACTION_TYPE:
      if (strcmp(attr_value, "NEW") == 0) {
        swap_p->action_type = NEW;
      } else if (strcmp(attr_value, "CORRECT") == 0) {
        swap_p->action_type = CORRECT;
      } else if (strcmp(attr_value, "CANCEL") == 0) {
        swap_p->action_type = CANCEL;
      } else {
        swap_p->action_type = ACTION_ERROR;
      }
      break;
    case TRANSACTION_TYPE:
      if (strcmp(attr_value, "Trade") == 0) {
        swap_p->transaction_type = TRADE;
      } else if (strcmp(attr_value, "Amendment") == 0) {
        swap_p->transaction_type = AMENDMENT;
      } else if (strcmp(attr_value, "Termination") == 0) {
        swap_p->transaction_type = TERMINATION;
      } else {
        swap_p->transaction_type = TRANSACTION_ERROR;
      }
      break;
    case IS_BLOCK_TRADE:
      if (strcmp(attr_value, "Y") == 0) {
        swap_p->is_block_trade = Y;
      } else if (strcmp(attr_value, "N") == 0) {
        swap_p->is_block_trade = N;
      } else {
        swap_p->is_block_trade = BLOCK_TRADE_ERROR;
      }
      break;
    case VENUE:
      if (strcmp(attr_value, "OFF") == 0) {
        swap_p->venue = OFF;
      } else if (strcmp(attr_value, "ON") == 0) {
        swap_p->venue = ON;
      } else {
        swap_p->venue = VENUE_ERROR;
      }
      break;
    case CURRENCY:
      if (strcmp(attr_value, "USD") == 0) {
        swap_p->currency = USD;
      } else if (strcmp(attr_value, "EUR") == 0) {
        swap_p->currency = EUR;
      } else {
        swap_p->currency = CURRENCY_ERROR;
      }
      break;
    case FLOAT_PAY_FREQ:
      swap_p->float_pay_freq = ParsePayFreq(attr_value);
      break;
    case FIXED_PAY_FREQ:
      swap_p->fixed_pay_freq = ParsePayFreq(attr_value);
    case REF_RATE:
      swap_p->ref_rate = ParseRefRate(attr_value);
    case PARSE_ERROR:
      break;
    default:
      break;
  }
}

// Reads data from one line of the csv file into a Swap structure
Swap SwapFromCSVLine(const char *input_line,  // input line of .csv text
                     int line_size,           // size in input_line
                     const char *data_cols,   // column names, assumed to have
                                              // been parsed from ParseLine
                     size_t max_colname_len  // size of each column name element
) {
  Swap swap;
  int max_buffer_size = 32;
  char buffer[max_buffer_size];
  int char_idx = 0;
  int col_idx = 0;
  int buffer_size = 0;
  char pay_freq_1[max_buffer_size], pay_freq_2[max_buffer_size];
  int col_1_is_float = 0;
  // Use this as a flag for parsing entries enclosed with " "
  int open_quotes = 0;
  while (char_idx < line_size) {
    const char this_char = input_line[char_idx];
    // if we reach a ',' and we have not found open quotes, we check out column
    // name and parse our buffer into the corresponding attribute of swap
    if ((this_char == ',') && (open_quotes == 0)) {
      buffer[buffer_size] = '\0';
      if (buffer_size > 0) {
        enum AttrToParse attr_name = EvaluateColname(&data_cols[col_idx]);
        if (attr_name == REF_RATE_IN_COL_1) {
          col_1_is_float = 1;
          swap.ref_rate = ParseRefRate(buffer);
        } else if (attr_name == REF_RATE_IN_COL_2) {
          col_1_is_float = 0;
          swap.ref_rate = ParseRefRate(buffer);
        } else if (attr_name == PAY_FREQ_1) {
          strncpy(pay_freq_1, buffer, buffer_size);
        } else if (attr_name == PAY_FREQ_2) {
          strncpy(pay_freq_2, buffer, buffer_size);
        } else {
          AssignSwapValue(&swap, attr_name, buffer);
        }
      }
      buffer_size = 0;
      col_idx = col_idx + max_colname_len;
    } else if (this_char == '"') {
      // if we run into a quote, we switch the flag and continue parsing,
      // ignoring the quote
      if (open_quotes == 0) {
        open_quotes = 1;
      } else {
        open_quotes = 0;
      }
    } else {
      buffer[buffer_size++] = this_char;
    }
    char_idx++;
  }
  if (col_1_is_float) {
    swap.float_pay_freq = ParsePayFreq(pay_freq_1);
    swap.fixed_pay_freq = ParsePayFreq(pay_freq_2);
  } else {
    swap.float_pay_freq = ParsePayFreq(pay_freq_2);
    swap.fixed_pay_freq = ParsePayFreq(pay_freq_1);
  }
  return swap;
}

Swap ListStringToSwap(const char *input, size_t input_size)
{
  // Expect a list to be passed in like "Colname:Value;"
  Swap swap = {0};
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
      if (buff_idx > max_buff_size)
        Die("Max buffer size reached in ListStringToSwap");
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
    enum AttrToParse parsed_attr_name = EvaluateColname(attr_name);
    AssignSwapValue(&swap, parsed_attr_name, attr_value);
  }
  return swap;
}

Swap SwapFromInputLine(const char *input_line)
{
  Swap swap = {0};
  size_t begin_idx = 0;
  size_t separator_idx = 0;
  size_t end_idx = 0;
  char attribute_buffer[64];
  char value_buffer[64];
  AttrToParse this_attribute;
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
      AssignSwapValue(&swap, this_attribute, value_buffer);
    }
  }
  return swap;
}


StartupContext LoadSwapsFromFile(const char *filename, int max_n_cols,
                                 int chunk_size, int max_colname_len,
                                 int max_n_loaded_swaps)
{
  Colnames colnames = {0};
  SwapList swap_list = {0};
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
