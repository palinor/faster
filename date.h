#ifndef uchar
typedef unsigned char uchar;
#endif
#ifndef uint
typedef unsigned int uint;
#endif

const uchar days_in_month[12] = {
	31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31
};

enum Month {
	BLANK, JANUARY, FEBRUARY, MARCH, APRIL, MAY, JUNE, JULY, AUGUST, SEPTEMBER, OCTOBER, NOVEMBER, DECEMBER
};

bool IsLeapYear(uint year);
uchar DaysInMonth(uchar month, uint year);
enum DayOfTheWeek {
	MONDAY, TUESDAY, WEDNESDAY, THURSDAY, FRIDAY, SATURDAY, SUNDAY
};

enum class Calendar {
	NONE, WEEKDAYS, SIFMA
};

bool IsWeekday(DayOfTheWeek day);

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

bool operator>(Date &date_1, Date &date_2);
bool operator>=(Date &date_1, Date &date_2);
// note(AION) after the fuck up with pointers, maybe we don't want to combine these with operator overloading
bool operator<(Date &date_1, Date &date_2);
bool operator<=(Date &date_1, Date &date_2);

bool operator==(Date &date_1, Date &date_2);

void IncrementDay(Date *date);
void IncrementMonth(Date *date);
void IncrementWeek(Date *date);
void AddCalendarDaysToDate(Date *date, int days);
int DistanceInDays(Date *left_date, Date *right_date);
int DistanceInWeekdays(Date *start_date, Date *end_date);

DayOfTheWeek GetDayOfTheWeek(Date *date);

int CalendarDays(Date *start_date, Date *end_date, Calendar calendar);