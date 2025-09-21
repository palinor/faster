#include <fstream>
#include <cassert>
#include <cmath>


struct TapeNodeFloat {
	float weights[2];
	size_t parents[2];
};

struct TapeNodeDouble {
	double weights[2];
	size_t parents[2];
};

struct TapeGradientFloat {
	float *derivatives;
	size_t number_of_derivatives;
};

struct TapeGradientDouble {
	double *derivatives;
	size_t number_of_derivatives;
};

struct TapeFloat {
	size_t current_node_index;
	TapeNodeFloat *nodes;
	size_t number_of_nodes;
};

struct TapeDouble {
	size_t current_node_index;
	TapeNodeDouble *nodes;
	size_t number_of_nodes;
};

struct AadFloat {
	TapeFloat *tape;
	size_t index_on_tape;
	float value;
};

struct AadDouble {
	TapeDouble *tape;
	size_t index_on_tape;
	double value;
};

int AadFloatInit(AadFloat *result, float x, TapeFloat *tape) {
	assert(tape->current_node_index + 1 < tape->number_of_nodes);
	result->index_on_tape = tape->current_node_index;
	result->value = x;
	result->tape = tape;
	TapeNodeFloat *new_node = tape->nodes + tape->current_node_index;
	++tape->current_node_index;
	new_node->weights[0] = 0;
	new_node->weights[1] = 0;
	new_node->parents[0] = 1;
	new_node->parents[1] = 1;
	return 0;
}

int AddSingleVariableStepToTape(AadFloat *result, size_t x_index, float f_x, float f_derivative, TapeFloat *tape) {
	assert(tape->current_node_index + 1 < tape->number_of_nodes);
	result->index_on_tape = tape->current_node_index;
	result->tape = tape;
	result->value = f_x;
	TapeNodeFloat *new_node = tape->nodes + tape->current_node_index;
	++tape->current_node_index;
	new_node->weights[0] = f_derivative;
	new_node->weights[1] = 0;
	new_node->parents[0] = x_index;
	new_node->parents[1] = 1;
	return 0;
}

int AddDoubleVariableStepToTape(AadFloat *result, size_t x_index, size_t y_index, float f_x_y, float f_derivative_x, float f_derivative_y, TapeFloat *tape) {
	assert(tape->current_node_index + 1 < tape->number_of_nodes);
	result->index_on_tape = tape->current_node_index;
	result->tape = tape;
	result->value = f_x_y;
	TapeNodeFloat *new_node = tape->nodes + tape->current_node_index;
	++tape->current_node_index;
	new_node->weights[0] = f_derivative_x;
	new_node->weights[1] = f_derivative_y;
	new_node->parents[0] = x_index;
	new_node->parents[1] = y_index;
	return 0;
}


AadFloat operator+(const AadFloat &x, const AadFloat &y) {
	assert(x.tape == y.tape);
	AadFloat result;
	AddDoubleVariableStepToTape(
		&result,
		x.index_on_tape,
		y.index_on_tape,
		x.value + y.value,
		1,
		1,
		x.tape
	);
	return result;
}

AadFloat operator+(float x, AadFloat y) {
	AadFloat x_new;
	AadFloatInit(&x_new, x, y.tape);
	return x_new + y;
}

AadFloat operator+(AadFloat x, float y) {
	AadFloat y_new;
	AadFloatInit(&y_new, y, x.tape);
	return x + y_new;
}

AadFloat operator-(const AadFloat &x, const AadFloat &y) {
	assert(x.tape == y.tape);
	AadFloat result;
	AddDoubleVariableStepToTape(
		&result,
		x.index_on_tape,
		y.index_on_tape,
		x.value - y.value,
		1,
		-1,
		x.tape);
	return result;
}

AadFloat operator-(float x, AadFloat y) {
	AadFloat x_new;
	AadFloatInit(&x_new, x, y.tape);
	return x_new - y;
}

AadFloat operator-(AadFloat x, float y) {
	AadFloat y_new;
	AadFloatInit(&y_new, y, x.tape);
	return x - y_new;
}

AadFloat operator*(AadFloat x, AadFloat y) {
	assert(x.tape == y.tape);
	AadFloat result;
	AddDoubleVariableStepToTape(
		&result,
		x.index_on_tape,
		y.index_on_tape,
		x.value * y.value,
		y.value,
		x.value,
		x.tape
	);
	return result;
}

AadFloat operator*(float x, AadFloat y) {
	AadFloat x_new;
	AadFloatInit(&x_new, x, y.tape);
	return x_new * y;
}

AadFloat operator*(AadFloat x, float y) {
	AadFloat y_new;
	AadFloatInit(&y_new, y, x.tape);
	return x * y_new;
}

AadFloat operator/(AadFloat x, AadFloat y) {
	assert(x.tape == y.tape);
	float fXY = x.value / y.value;
	float dfX = 1 / y.value;
	float dfY = -x.value / (y.value * y.value);
	AadFloat result;
	AddDoubleVariableStepToTape(
		&result,
		x.index_on_tape,
		y.index_on_tape,
		fXY,
		dfX,
		dfY,
		x.tape
	);
	return result;
}

AadFloat operator/(float x, AadFloat y) {
	AadFloat x_new;
	AadFloatInit(&x_new, x, y.tape);
	return x_new / y;
}

AadFloat operator/(AadFloat x, float y) {
	AadFloat y_new;
	AadFloatInit(&y_new, y, x.tape);
	return x / y_new;
}

AadFloat exp(AadFloat x) {
	AadFloat result;
	AddSingleVariableStepToTape(
		&result,
		x.index_on_tape,
		std::expf(x.value),
		std::expf(x.value),
		x.tape
	);
	return result;
}

AadFloat log(AadFloat x) {
	AadFloat result;
	AddSingleVariableStepToTape(
		&result,
		x.index_on_tape,
		std::logf(x.value),
		1 / x.value,
		x.tape
	);
	return result;
}

AadFloat pow(AadFloat x, AadFloat y) {
	assert(x.tape == y.tape);
	return exp(y * log(x));
}

AadFloat pow(AadFloat x, float y) {
	AadFloat y_new;
	AadFloatInit(&y_new, y, x.tape);
	return pow(x, y_new);
}

AadFloat sqrt(AadFloat x) {
	AadFloat result;
	AddSingleVariableStepToTape(
		&result,
		x.index_on_tape,
		std::sqrtf(x.value),
		1 / (2 * std::sqrtf(x.value)),
		x.tape
	);
	return result;
}

float GradientWithRespectTo(TapeGradientFloat *grad, AadFloat x) {
	return grad->derivatives[x.index_on_tape];
}

int PopulateTapeGradient(TapeGradientFloat *result, AadFloat *final_value) {
	size_t number_of_derivatives = final_value ->index_on_tape + 1;
	void *new_derivative_block = calloc(number_of_derivatives, sizeof(float));
	if (!new_derivative_block) {
		return -1;
	}
	result->derivatives = reinterpret_cast<float *>(new_derivative_block);
	result->number_of_derivatives = number_of_derivatives;
	result->derivatives[number_of_derivatives - 1] = 1;
	for (size_t k = 1; k <= number_of_derivatives; k++) {
		size_t i = number_of_derivatives - k;
		TapeNodeFloat *node = final_value->tape->nodes + i;
		float derivative = result->derivatives[i];
		result->derivatives[node->parents[0]] += node->weights[0] * derivative;
		result->derivatives[node->parents[1]] += node->weights[1] * derivative;
	}
	return 0;
}

void AadTests() {
	TapeFloat t;
	const size_t maxTapeSize = 256;
	t.nodes = reinterpret_cast<TapeNodeFloat *>(malloc(maxTapeSize * sizeof(TapeNodeFloat)));
	t.current_node_index = 0;
	t.number_of_nodes = maxTapeSize;
	AadFloat x, y, z1, z2;
	AadFloatInit(&x, 9, &t);
	AadFloatInit(&y, 0.5, &t);
	z1 = x * y;
	z2 = pow(x, y);

	TapeGradientFloat grad1, grad2;
	PopulateTapeGradient(&grad1, &z1);
	PopulateTapeGradient(&grad2, &z2);


	float dz1dx = GradientWithRespectTo(&grad1, x);
	float dz1dy = GradientWithRespectTo(&grad1, y);
	float dz2dx = GradientWithRespectTo(&grad2, x);
	float dz2dy = GradientWithRespectTo(&grad2, y);

	assert(fabs(dz1dx - y.value) < 1e-15);
	assert(fabs(dz1dy - x.value) < 1e-15);
	assert(fabs(dz2dy - (log(x) * z2).value) < 1e-15);
	assert(fabs(dz2dx - 1 / (2 * z2.value)) < 1e-15);

	return;
}

struct AadFloat_shifted_sabr_params {
	AadFloat forward_;
	AadFloat sigma0_;
	AadFloat alpha_;
	AadFloat beta_;
	AadFloat rho_;
	AadFloat zeta_;
	AadFloat timeToExpiry_;
};

//todo(AION): there is a bug in here somewhere. Need to double check the formula:w
AadFloat AadFloatShiftedSabr(AadFloat &strike, AadFloat_shifted_sabr_params *params) {
	AadFloat moneyness = params->forward_ - strike;
	AadFloat shiftedForward = params->forward_ + params->zeta_;
	AadFloat shiftedStrike = strike + params->zeta_;
	AadFloat i0;
	if (fabs(moneyness.value) < 1e-8) {
		i0 = params->sigma0_ * pow(shiftedForward, params->beta_);
	}
	else {
		AadFloat z;
		if (params->beta_.value > 1 - 1e-8) {
			z = params->alpha_ * moneyness / params->sigma0_;
		}
		else {
			AadFloat betaInv = 1 - params->beta_;
			z = params->alpha_ / params->sigma0_ * (pow(shiftedForward, betaInv) - pow(shiftedStrike, betaInv)) / betaInv;
		}
		AadFloat inSqrtTerm = 1 - 2 * params->rho_ * z + z * z;
		i0 = params->alpha_ * moneyness / log((sqrt(inSqrtTerm) + z - params->rho_) / (1 - params->rho_));

	}

	AadFloat beta2 = params->beta_ * params->beta_;
	AadFloat sigma02 = params->sigma0_ * params->sigma0_;
	AadFloat rho2 = params->rho_ * params->rho_;
	AadFloat alpha2 = params->alpha_ * params->alpha_;
	AadFloat i1Term1 = ((beta2 - 2 * params->beta_) * sigma02) / (24 * pow((shiftedStrike + shiftedForward) / 2, 2 * (1 - params->beta_)));
	AadFloat i1Term2 = (params->rho_ * params->alpha_ * params->beta_ * params->sigma0_) / (4 * pow((shiftedStrike + shiftedForward) / 2, 1 - params->beta_));
	AadFloat i1Term3 = ((2 - 3 * rho2) * alpha2) / 24;
	AadFloat i1 = i1Term1 + i1Term2 + i1Term3;
	return i0 * (1 + i1 * params->timeToExpiry_);

}

void TestShiftedSabr() {
	TapeFloat t;
	const size_t maxTapeSize = 1024;
	float strikeStart = -2e-2;
	float strikeEnd = 1e-1;
	size_t nStrikes = 100;
	float strikeStep = (strikeEnd - strikeStart) / nStrikes;
	std::ofstream myFile;
	myFile.open("shifted_sabr_smile.csv");
	for (size_t i = 0; i < nStrikes; i++) {
		t.nodes = reinterpret_cast<TapeNodeFloat *>(malloc(maxTapeSize * sizeof(TapeNodeFloat)));
		t.current_node_index = 0;
		t.number_of_nodes = maxTapeSize;
		AadFloat_shifted_sabr_params params;
		AadFloatInit(&(params.forward_), 2.5e-2, &t);
		AadFloatInit(&(params.timeToExpiry_), 10, &t);
		AadFloatInit(&(params.alpha_), 0.28, &t);
		AadFloatInit(&(params.beta_), 0.36, &t);
		AadFloatInit(&(params.rho_), 0.05, &t);
		AadFloatInit(&(params.zeta_), 0.03, &t);
		AadFloatInit(&(params.sigma0_), 0.019, &t);
		AadFloat thisStrike, impliedVol;
		AadFloatInit(&thisStrike, strikeStart + i * strikeStep, &t);
		impliedVol = AadFloatShiftedSabr(thisStrike, &params);
		myFile << thisStrike.value << "," << impliedVol.value << "\n";
		free(t.nodes);
	}
	myFile.close();
}

