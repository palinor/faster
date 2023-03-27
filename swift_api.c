#include <stdlib.h>
#include <string.h>
#include "swift_api.h"

yield_curve PlaceholderYieldCurve(const char *tenor, float value) {
	yield_curve result;
	result.tenor = malloc(strlen(tenor) * sizeof(char));
	strcpy(result.tenor, tenor);
	result.value = value;
	return result;
}