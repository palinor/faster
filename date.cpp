typedef unsigned char uchar;
typedef unsigned int uint;

const uchar days_in_month[12] = {
	31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31
};

enum Month {
	BLANK, JANUARY, FEBRUARY, MARCH, APRIL, MAY, JUNE, JULY, AUGUST, SEPTEMBER, OCTOBER, NOVEMBER, DECEMBER
};

bool IsLeapYear(uint year) {
	if (!(year % 400)) {
		return true;
	} else if (!(year % 100)) {
		return false;
	} else {
		return !(year % 4);
	}
}

uchar DaysInMonth(uchar month, uint year) {
	if (month == FEBRUARY) {
		return IsLeapYear(year) ? 29 : 28;
	} else {
		return days_in_month[month];
	}
}


enum DayOfTheWeek {
	MONDAY, TUESDAY, WEDNESDAY, THURSDAY, FRIDAY, SATURDAY, SUNDAY
};

enum class Calendar {
	NONE, WEEKDAYS, SIFMA
};

bool IsWeekday(DayOfTheWeek day) {
	return day < SATURDAY;
}

struct Date {
	uchar day;
	uchar month;
	uint year;
};

struct Time {
	uchar hour;
	uchar minute;
	uchar second;
	uint millisecond;
};

bool operator>(Date &date_1, Date &date_2) {
	if (date_1.year != date_2.year) {
		return date_1.year > date_2.year;
	}
	if (date_1.month != date_2.month) {
		return date_1.month > date_2.month;
	}
	return date_1.day > date_2.day;
}

bool operator>=(Date &date_1, Date &date_2) {
	if (date_1.year != date_2.year) {
		return date_1.year > date_2.year;
	}
	if (date_2.month != date_2.month) {
		return date_1.month > date_2.month;
	}
	return date_1.day >= date_2.day;
}

// note(AION) after the fuck up with pointers, maybe we don't want to combine these with operator overloading
bool operator<(Date &date_1, Date &date_2) {
	return date_2 > date_1;
}

bool operator<=(Date &date_1, Date &date_2) {
	return date_2 >= date_1;
}

bool operator==(Date &date_1, Date &date_2) {
	return (date_1.year == date_2.year && date_1.month == date_2.month && date_1.day == date_2.day);
}


void IncrementDay(Date *date) {
	if (date->day == DaysInMonth(date->month, date->year)) {
		if (date->month == 12) {
			++date->year;
			date->month = JANUARY;
			date->day = 1;
			return;
		} else {
			++date->month;
			date->day = 1;
			return;
		}
	} else {
		++date->day;
		return;
	}
}

void IncrementMonth(Date *date) {
	if (date->month == DECEMBER) {
		++date->year;
		date->month = JANUARY;
	} else {
		++date->month;
	}
}

void IncrementWeek(Date *date) {
	date->day += 7;
	if (date->day > DaysInMonth(date->month, date->year)) {
		date->day -= DaysInMonth(date->month, date->year);
		IncrementMonth(date);
	}
}

void AddCalendarDaysToDate(Date *date, int days) {
	int days_left_to_add = days;
	while (days_left_to_add > 366) {
		bool will_cross_leap_year = date->month > FEBRUARY ? IsLeapYear(date->year + 1) : IsLeapYear(date->year);
		++date->year;
		days_left_to_add -= will_cross_leap_year ? 366 : 365;
	}
	while (days_left_to_add > 31) {
		days_left_to_add -= DaysInMonth(date->month, date->year);
		IncrementMonth(date);
	}
	if (days_left_to_add > DaysInMonth(date->month, date->year) - date->day) {
		IncrementMonth(date);
		days_left_to_add -= DaysInMonth(date->month, date->year) - date->day;
		date->day = 1;
	}
	while (days_left_to_add) {
		IncrementDay(date);
		--days_left_to_add;
	}
}

int DistanceInDays(Date *left_date, Date *right_date) {
	if (*left_date > *right_date) return -1 * DistanceInDays(right_date, left_date);
	int result = 0;
	uint left_year = left_date->year;
	uchar left_month = left_date->month;
	uchar left_day = left_date->day;
	if (left_year < right_date->year) {
		result += DaysInMonth(left_month, left_year) - left_day;
		left_day = 1;
		++left_month;
		while (left_month < 12) {
			result += DaysInMonth(left_month++, left_year);
		}
		left_month = 1;
		++left_year;
	}
	while (left_year < right_date->year) {
		result += IsLeapYear(left_year++) ? 366 : 365;
	}
	while (left_month < right_date->month) {
		result += DaysInMonth(left_month++, left_year);
	}
	result += right_date->day - left_day;
	return result;
}

int DistanceInWeekdays(Date *start_date, Date *end_date) {
	if (*start_date > *end_date) {
		return -1 * DistanceInWeekdays(end_date, start_date);
	}
	Date known_monday = {6, MAY, 2024};
	int weekday = DistanceInDays(&known_monday, start_date) % 7;
	int total_distance = DistanceInDays(start_date, end_date);
	if (weekday < 0) {
		weekday += 7;
	}
	int number_of_weeks = total_distance / 7;
	int remaining_days = total_distance % 7;
	int remaining_weekdays = 0;
	for (int day = 0; day < remaining_days; ++day) {
		++weekday;
		if (weekday < 5) {
			++remaining_weekdays;
		}
		weekday %= 7;
	}
	return number_of_weeks * 5 + remaining_weekdays;
}

DayOfTheWeek GetDayOfTheWeek(Date *date) {
	Date known_monday = {6, MAY, 2024};
	return (DayOfTheWeek)(DistanceInDays(&known_monday, date) % 7);
}

inline int NewYearsCount(Date *start_date, Date *end_date) {
	if (*start_date > *end_date) {
		return 0;
	}
	if ((end_date->year == start_date->year) && ((start_date->month > 1) || (start_date->day > 2))) {
		return 0;
	}
	Date moving_date;
	moving_date.year = start_date->year;
	moving_date.month = start_date->month;
	moving_date.day = start_date->day;
	int result = 0;
	Date known_monday = {6, MAY, 2024};
	Date next_new_years;
	while (moving_date < *end_date) {
		uint next_year_to_check = ((moving_date.month == 1) && (moving_date.day < 3)) ? moving_date.year : moving_date.year + 1;
		next_new_years.year = next_year_to_check;
		next_new_years.month = JANUARY;
		next_new_years.day = 1;
		DayOfTheWeek next_new_years_day_of_the_week = (DayOfTheWeek)(DistanceInDays(&known_monday, &next_new_years) % 7);
		if (next_new_years_day_of_the_week != SATURDAY) {
			if (next_new_years_day_of_the_week == SUNDAY) {
				next_new_years.day = 2;
			}
			if (moving_date <= next_new_years) {
				++result;
			}
		}
		++moving_date.year;
	}
	return result;
}

inline int MlkDayCount(Date *start_date, Date *end_date) {
	if (*end_date < *start_date) {
		return 0;
	}
	if ((start_date->year == end_date->year) && (start_date->month > JANUARY)) {
		return 0;
	}
	Date moving_date;
	moving_date.year = start_date->year;
	moving_date.month = start_date->month;
	moving_date.day = start_date->day;
	int result = 0;
	Date known_monday = {6, MAY, 2024};
	Date mlk_day;
	while (moving_date < *end_date) {
		mlk_day.year = moving_date.year;
		mlk_day.month = JANUARY;
		mlk_day.day = 1;
		int day_of_the_week = DistanceInDays(&known_monday, &mlk_day) % 7;
		mlk_day.day = (7 - day_of_the_week) % 7 + 15;
		if ((moving_date <= mlk_day) && (mlk_day <= *end_date)) {
			++result;
		}	
		++moving_date.year;
	}
	return result;
}

inline int PresidentsDayCount(Date *start_date, Date *end_date) {
	if (*end_date < *start_date) {
		return 0;
	}
	if ((start_date->year == end_date->year) && (start_date->month > 2)) {
		return 0;
	}
	Date moving_date;
	moving_date.year = start_date->year;
	moving_date.month = start_date->month;
	moving_date.day = start_date->day;
	int result = 0;
	Date known_monday = {6, MAY, 2024};
	Date presidents_day;
	while (moving_date < *end_date) {
		presidents_day.year = moving_date.year;
		presidents_day.month = FEBRUARY;
		presidents_day.day = 1;
		int day_of_the_week = DistanceInDays(&known_monday, &presidents_day) % 7;
		presidents_day.day = (7 - day_of_the_week) % 7 + 15;
		if ((moving_date <= presidents_day) && (presidents_day <= *end_date)) {
			++result;
		}	
		++moving_date.year;
	}
	return result;
}

inline int MemorialDayCount(Date *start_date, Date *end_date) {
	if (*end_date < *start_date) {
		return 0;
	}
	if ((start_date->year == end_date->year) && (start_date->month > 5)) {
		return 0;
	}
	Date moving_date;
	moving_date.year = start_date->year;
	moving_date.month = start_date->month;
	moving_date.day = start_date->day;
	int result = 0;
	Date known_monday = {6, MAY, 2024};
	Date memorial_day;
	while (moving_date < *end_date) {
		memorial_day.year = moving_date.year;
		memorial_day.month = MAY;
		memorial_day.day = 1;
		int day_of_the_week = DistanceInDays(&known_monday, &memorial_day) % 7;
		memorial_day.day = (7 - day_of_the_week) % 7 + 1;
		while (memorial_day.day + 7 < DaysInMonth(MAY, memorial_day.year)) {
			memorial_day.day += 7;
		}
		if ((moving_date <= memorial_day) && (memorial_day <= *end_date)) {
			++result;
		}	
		++moving_date.year;
	}
	return result;
}

inline int JuneteenthCount(Date *start_date, Date *end_date) {
	if (*start_date > *end_date) {
		return 0;
	}
	if ((end_date->year == start_date->year) && (start_date->month > 7)) {
		return 0;
	}
	Date moving_date;
	moving_date.year = start_date->year;
	moving_date.month = start_date->month;
	moving_date.day = start_date->day;
	int result = 0;
	Date known_monday = {6, MAY, 2024};
	Date juneteenth;
	while (moving_date < *end_date) {
		juneteenth.year = moving_date.year;
		juneteenth.month = JUNE;
		juneteenth.day = 19;
		DayOfTheWeek juneteenth_day = (DayOfTheWeek)(DistanceInDays(&known_monday, &juneteenth) % 7);
		if (juneteenth_day == SATURDAY) {
			--juneteenth.day;
		}
		if (juneteenth_day == SUNDAY) {
			++juneteenth.day;
		}
		if ((moving_date <= juneteenth) && (juneteenth <= *end_date)) {
			++result;
		}
		++moving_date.year;
	}
	return result;
}

inline int July4Count(Date *start_date, Date *end_date) {
	if (*start_date > *end_date) {
		return 0;
	}
	if ((end_date->year == start_date->year) && (start_date->month > 7)) {
		return 0;
	}
	Date moving_date;
	moving_date.year = start_date->year;
	moving_date.month = start_date->month;
	moving_date.day = start_date->day;
	int result = 0;
	Date known_monday = {6, MAY, 2024};
	Date july_fourth;
	while (moving_date < *end_date) {
		july_fourth.year = moving_date.year;
		july_fourth.month = JULY;
		july_fourth.day = 4;
		DayOfTheWeek july_fourth_day = (DayOfTheWeek)(DistanceInDays(&known_monday, &july_fourth) % 7);
		if (july_fourth_day == SATURDAY) {
			--july_fourth.day;
		}
		if (july_fourth_day == SUNDAY) {
			++july_fourth.day;
		}
		if ((moving_date <= july_fourth) && (july_fourth <= *end_date)) {
			++result;
		}
		++moving_date.year;
	}
	return result;
}

inline int LaborDayCount(Date *start_date, Date *end_date) {
	if (*end_date < *start_date) {
		return 0;
	}
	if ((start_date->year == end_date->year) && (start_date->month > 5)) {
		return 0;
	}
	Date moving_date;
	moving_date.year = start_date->year;
	moving_date.month = start_date->month;
	moving_date.day = start_date->day;
	int result = 0;
	Date known_monday = {6, MAY, 2024};
	Date labor_day;
	while (moving_date < *end_date) {
		labor_day.year = moving_date.year;
		labor_day.month = SEPTEMBER;
		labor_day.day = 1;
		int day_of_the_week = DistanceInDays(&known_monday, &labor_day) % 7;
		labor_day.day = (7 - day_of_the_week) % 7 + 1;
		if ((moving_date <= labor_day) && (labor_day <= *end_date)) {
			++result;
		}	
		++moving_date.year;
	}
	return result;
}

inline int ColumbusDayCount(Date *start_date, Date *end_date) {
	if (*end_date < *start_date) {
		return 0;
	}
	if ((start_date->year == end_date->year) && (start_date->month > 5)) {
		return 0;
	}
	Date moving_date;
	moving_date.year = start_date->year;
	moving_date.month = start_date->month;
	moving_date.day = start_date->day;
	int result = 0;
	Date known_monday = {6, MAY, 2024};
	Date columbus_day;
	while (moving_date < *end_date) {
		columbus_day.year = moving_date.year;
		columbus_day.month = OCTOBER;
		columbus_day.day = 1;
		int day_of_the_week = DistanceInDays(&known_monday, &columbus_day) % 7;
		columbus_day.day = (7 - day_of_the_week) % 7 + 8;
		if ((moving_date <= columbus_day) && (columbus_day <= *end_date)) {
			++result;
		}	
		++moving_date.year;
	}
	return result;
}

inline int VeteransDayCount(Date *start_date, Date *end_date) {
	if (*start_date > *end_date) {
		return 0;
	}
	if ((end_date->year == start_date->year) && (start_date->month > 7)) {
		return 0;
	}
	Date moving_date;
	moving_date.year = start_date->year;
	moving_date.month = start_date->month;
	moving_date.day = start_date->day;
	int result = 0;
	Date known_monday = {6, MAY, 2024};
	Date veterans_day;
	while (moving_date < *end_date) {
		veterans_day.year = moving_date.year;
		veterans_day.month = JULY;
		veterans_day.day = 4;
		DayOfTheWeek veterans_day_of_the_week = (DayOfTheWeek)(DistanceInDays(&known_monday, &veterans_day) % 7);
		if (veterans_day_of_the_week == SATURDAY) {
			--veterans_day.day;
		}
		if (veterans_day_of_the_week == SUNDAY) {
			++veterans_day.day;
		}
		if ((moving_date <= veterans_day) && (veterans_day <= *end_date)) {
			++result;
		}
		++moving_date.year;
	}
	return result;
}

inline int ThanksgivingDayCount(Date *start_date, Date *end_date) {
	if (*end_date < *start_date) {
		return 0;
	}
	if ((start_date->year == end_date->year) && (start_date->month > 11)) {
		return 0;
	}
	Date moving_date;
	moving_date.year = start_date->year;
	moving_date.month = start_date->month;
	moving_date.day = start_date->day;
	int result = 0;
	Date known_monday = {6, MAY, 2024};
	Date thanksgiving;
	while (moving_date < *end_date) {
		thanksgiving.year = moving_date.year;
		thanksgiving.month = OCTOBER;
		thanksgiving.day = 1;
		int day_of_the_week = DistanceInDays(&known_monday, &thanksgiving) % 7;
		thanksgiving.day = (3 - day_of_the_week) % 7 + 15;
		if ((moving_date <= thanksgiving) && (thanksgiving <= *end_date)) {
			++result;
		}	
		++moving_date.year;
	}
	return result;
}

inline int ChristmasDayCount(Date *start_date, Date *end_date) {
	if (*start_date > *end_date) {
		return 0;
	}
	if (*start_date == *end_date) {
		return 0;
	}
	Date moving_date;
	moving_date.year = start_date->year;
	moving_date.month = start_date->month;
	moving_date.day = start_date->day;
	int result = 0;
	Date known_monday = {6, MAY, 2024};
	Date christmas_day;
	while (moving_date < *end_date) {
		christmas_day.year = moving_date.year;
		christmas_day.month = DECEMBER;
		christmas_day.day = 25;
		DayOfTheWeek christmas_day_of_the_week = (DayOfTheWeek)(DistanceInDays(&known_monday, &christmas_day) % 7);
		if (christmas_day_of_the_week == SATURDAY) {
			--christmas_day.day;
		}
		if (christmas_day_of_the_week == SUNDAY) {
			++christmas_day.day;
		}
		if ((moving_date <= christmas_day) && (christmas_day <= *end_date)) {
			++result;
		}
		++moving_date.year;
	}
	return result;
}

int CalendarDays(Date *start_date, Date *end_date, Calendar calendar) {
	if (*end_date < *start_date) {
		return -1 * CalendarDays(end_date, start_date, calendar);
	}
	switch (calendar) {
	case Calendar::NONE: {
		return DistanceInDays(start_date, end_date);
	} break;
	case Calendar::WEEKDAYS: {
		return DistanceInWeekdays(start_date, end_date);
	} break;
	case Calendar::SIFMA: {
		int number_of_weekdays = DistanceInWeekdays(start_date, end_date);
		number_of_weekdays -= NewYearsCount(start_date, end_date);
		number_of_weekdays -= MlkDayCount(start_date, end_date);
		number_of_weekdays -= PresidentsDayCount(start_date, end_date);
		number_of_weekdays -= MemorialDayCount(start_date, end_date);
		number_of_weekdays -= JuneteenthCount(start_date, end_date);
		number_of_weekdays -= July4Count(start_date, end_date);
		number_of_weekdays -= LaborDayCount(start_date, end_date);
		number_of_weekdays -= ColumbusDayCount(start_date, end_date);
		number_of_weekdays -= VeteransDayCount(start_date, end_date);
		number_of_weekdays -= ThanksgivingDayCount(start_date, end_date);
		number_of_weekdays -= ChristmasDayCount(start_date, end_date);
		return number_of_weekdays;
	} break;
	default: {
		return DistanceInDays(start_date, end_date);
	}
	}
}

inline float CoverageAct365(Date *start_date, Date *end_date) {
	return ((float)DistanceInDays(start_date, end_date)) / 365.0f;
}

inline float CoverageAct360(Date *start_date, Date *end_date) {
	return ((float)DistanceInDays(start_date, end_date)) / 360.0f;
}

inline float Coverage30360(Date *start_date, Date *end_date) {
	int number_of_years = end_date->year - start_date->year;
	int number_of_months = end_date->month - start_date->month;
	// we just assume that the day will line up if we are using coverages where months are all 30 days 
	return (float)number_of_years + (float)number_of_months * 30.0f / 360.0f;
}

inline float Coverage30365(Date *start_date, Date *end_date) {
	int number_of_years = end_date->year - start_date->year;
	int number_of_months = end_date->month - start_date->month;
	// we just assume that the day will line up if we are using coverages where months are all 30 days 
	return (float)number_of_years + (float)number_of_months * 30.0f / 365.0f;
}
