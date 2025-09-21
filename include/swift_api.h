#include <math.h>

typedef struct yield_curve {
	char *tenor;
	float value;
} yield_curve;

yield_curve PlaceholderYieldCurve(const char *tenor, float value);
