#include <fstream>
#include <cassert>
#include <cmath>

const size_t TAPE_SIZE = 100;
struct float_tape_node {
	double weights[2];
	size_t parents[2];
};

struct float_grad {
	float *derivs_;
	size_t derivsSize_;
};

struct float_tape {
	size_t node_index_;
	float_tape_node *nodes_;
	size_t nNodes_;
};

struct aad_float {
	float_tape *tape_;
	size_t index_;
	float value_;
};

aad_float AddFloatToTape(float x, float_tape *tape) {
	assert(tape->node_index_ + 1 < tape->nNodes_);
	aad_float result;
	result.index_ = tape->node_index_;
	result.tape_ = tape;
	result.value_ = x;
	float_tape_node *newNode = tape->nodes_ + tape->node_index_;
	tape->node_index_++;
	newNode->weights[0] = 0;
	newNode->weights[1] = 0;
	newNode->parents[0] = 1;
	newNode->parents[1] = 1;
	return result;
}

aad_float AddFunctionOfXToTape(size_t xIndex, float fX, float fDerivative, float_tape *tape) {
	assert(tape->node_index_ + 1 < tape->nNodes_);
	aad_float result;
	result.index_ = tape->node_index_;
	result.tape_ = tape;
	result.value_ = fX;
	float_tape_node *newNode = tape->nodes_ + tape->node_index_;
	tape->node_index_++;
	newNode->weights[0] = fDerivative;
	newNode->weights[1] = 0;
	newNode->parents[0] = xIndex;
	newNode->parents[1] = 1;
	return result;
}

aad_float AddFunctionOfXYToTape(size_t xIndex, size_t yIndex, float fXY, float fDerivativeX, float fDerivativeY, float_tape *tape) {
	assert(tape->node_index_ + 1 < tape->nNodes_);
	aad_float result;
	result.index_ = tape->node_index_;
	result.tape_ = tape;
	result.value_ = fXY;
	float_tape_node *newNode = tape->nodes_ + tape->node_index_;
	tape->node_index_++;
	newNode->weights[0] = fDerivativeX;
	newNode->weights[1] = fDerivativeY;
	newNode->parents[0] = xIndex;
	newNode->parents[1] = yIndex;
	return result;
}


aad_float operator+(aad_float x, aad_float y) {
	assert(x.tape_ == y.tape_);
	aad_float result = AddFunctionOfXYToTape(x.index_, y.index_, x.value_ + y.value_, 1, 1, x.tape_);
	return result;
}

aad_float operator+(float x, aad_float y) {
	aad_float xNew = AddFloatToTape(x, y.tape_);
	return xNew + y;
}

aad_float operator+(aad_float x, float y) {
	aad_float yNew = AddFloatToTape(y, x.tape_);
	return x + yNew;
}

aad_float operator-(aad_float x, aad_float y) {
	assert(x.tape_ == y.tape_);
	return AddFunctionOfXYToTape(x.index_, y.index_, x.value_ - y.value_, 1, -1, x.tape_);
}

aad_float operator-(float x, aad_float y) {
	aad_float xNew = AddFloatToTape(x, y.tape_);
	return xNew - y;
}

aad_float operator-(aad_float x, float y) {
	aad_float yNew = AddFloatToTape(y, x.tape_);
	return x - yNew;
}

aad_float operator*(aad_float x, aad_float y) {
	assert(x.tape_ == y.tape_);
	return AddFunctionOfXYToTape(x.index_, y.index_, x.value_ * y.value_, y.value_, x.value_, x.tape_);
}

aad_float operator*(float x, aad_float y) {
	aad_float xNew = AddFloatToTape(x, y.tape_);
	return xNew * y;
}

aad_float operator*(aad_float x, float y) {
	aad_float yNew = AddFloatToTape(y, x.tape_);
	return x * yNew;
}

aad_float operator/(aad_float x, aad_float y) {
	assert(x.tape_ == y.tape_);
	float fXY = x.value_ / y.value_;
	float dfX = 1 / y.value_;
	float dfY = -x.value_ / (y.value_ * y.value_);
	return AddFunctionOfXYToTape(x.index_, y.index_, fXY, dfX, dfY, x.tape_);
}

aad_float operator/(float x, aad_float y) {
	aad_float xNew = AddFloatToTape(x, y.tape_);
	return xNew / y;
}

aad_float operator/(aad_float x, float y) {
	aad_float yNew = AddFloatToTape(y, x.tape_);
	return x / yNew;
}

aad_float exp(aad_float x) {
	return AddFunctionOfXToTape(x.index_, std::exp(x.value_), std::exp(x.value_), x.tape_);
}

aad_float log(aad_float x) {
	return AddFunctionOfXToTape(x.index_, std::log(x.value_), 1 / x.value_, x.tape_);
}

aad_float pow(aad_float x, aad_float y) {
	assert(x.tape_ == y.tape_);
	return exp(y * log(x));
}

aad_float pow(aad_float x, float y) {
	aad_float yNew = AddFloatToTape(y, x.tape_);
	return pow(x, yNew);
}

aad_float sqrt(aad_float x) {
	return AddFunctionOfXToTape(x.index_, std::sqrt(x.value_), 1 / (2 * std::sqrt(x.value_)), x.tape_);
}

float GradientWithRespectTo(float_grad *grad, aad_float x) {
	return grad->derivs_[x.index_];
}

float_grad ComputeGradient(aad_float *finalValue) {
	float_grad result;
	size_t nGrads = finalValue->index_ + 1;
	result.derivs_ = reinterpret_cast<float *>(calloc(nGrads, sizeof(float)));
	result.derivsSize_ = nGrads;
	result.derivs_[nGrads - 1] = 1;
	for (size_t k = 1; k <= result.derivsSize_; k++) {
		size_t i = result.derivsSize_ - k;
		float_tape_node node = finalValue->tape_->nodes_[i];
		float deriv = result.derivs_[i];
		result.derivs_[node.parents[0]] += node.weights[0] * deriv;
		result.derivs_[node.parents[1]] += node.weights[1] * deriv;
	}
	return result;
}

void AadTests() {
	float_tape t;
	const size_t maxTapeSize = 256;
	t.nodes_ = reinterpret_cast<float_tape_node *>(malloc(maxTapeSize * sizeof(float_tape_node)));
	t.node_index_ = 0;
	t.nNodes_ = maxTapeSize;
	aad_float x = AddFloatToTape(9, &t);
	aad_float y = AddFloatToTape(0.5, &t);
	aad_float z1 = x * y;
	aad_float z2 = pow(x, y);

	float_grad grad1 = ComputeGradient(&z1);
	float_grad grad2 = ComputeGradient(&z2);
	float dz1dx = GradientWithRespectTo(&grad1, x);
	float dz1dy = GradientWithRespectTo(&grad1, y);
	float dz2dx = GradientWithRespectTo(&grad2, x);
	float dz2dy = GradientWithRespectTo(&grad2, y);

	assert(fabs(dz1dx - y.value_) < 1e-15);
	assert(fabs(dz1dy - x.value_) < 1e-15);
	assert(fabs(dz2dy - (log(x) * z2).value_) < 1e-15);
	assert(fabs(dz2dx - 1 / (2 * z2.value_)) < 1e-15);

	return;
}

struct aad_float_shifted_sabr_params {
	aad_float forward_;
	aad_float sigma0_;
	aad_float alpha_;
	aad_float beta_;
	aad_float rho_;
	aad_float zeta_;
	aad_float timeToExpiry_;
};

//todo(AION): there is a bug in here somewhere. Need to double check the formula:w
aad_float AadFloatShiftedSabr(aad_float strike, aad_float_shifted_sabr_params *params) {
	aad_float moneyness = params->forward_ - strike;
	aad_float shiftedForward = params->forward_ + params->zeta_;
	aad_float shiftedStrike = strike + params->zeta_;
	aad_float i0;
	if (fabs(moneyness.value_) < 1e-8) {
		i0 = params->sigma0_ * pow(shiftedForward, params->beta_);
	}
	else {
		aad_float z;
		if (params->beta_.value_ > 1 - 1e-8) {
			z = params->alpha_ * moneyness / params->sigma0_;
		}
		else {
			aad_float betaInv = 1 - params->beta_;
			z = params->alpha_ / params->sigma0_ * (pow(shiftedForward, betaInv) - pow(shiftedStrike, betaInv)) / betaInv;
		}
		aad_float inSqrtTerm = 1 - 2 * params->rho_ * z + z * z;
		i0 = params->alpha_ * moneyness / log((sqrt(inSqrtTerm) + z - params->rho_) / (1 - params->rho_));

	}

	aad_float beta2 = params->beta_ * params->beta_;
	aad_float sigma02 = params->sigma0_ * params->sigma0_;
	aad_float rho2 = params->rho_ * params->rho_;
	aad_float alpha2 = params->alpha_ * params->alpha_;
	aad_float i1Term1 = ((beta2 - 2 * params->beta_) * sigma02) / (24 * pow((shiftedStrike + shiftedForward) / 2, 2 * (1 - params->beta_)));
	aad_float i1Term2 = (params->rho_ * params->alpha_ * params->beta_ * params->sigma0_) / (4 * pow((shiftedStrike + shiftedForward) / 2, 1 - params->beta_));
	aad_float i1Term3 = ((2 - 3 * rho2) * alpha2) / 24;
	aad_float i1 = i1Term1 + i1Term2 + i1Term3;
	return i0 * (1 + i1 * params->timeToExpiry_);

}

void TestShiftedSabr() {
	float_tape t;
	const size_t maxTapeSize = 1024;
	float strikeStart = -2e-2;
	float strikeEnd = 1e-1;
	size_t nStrikes = 100;
	float strikeStep = (strikeEnd - strikeStart) / nStrikes;
	std::ofstream myFile;
	myFile.open("shifted_sabr_smile.csv");
	for (size_t i = 0; i < nStrikes; i++) {
		t.nodes_ = reinterpret_cast<float_tape_node *>(malloc(maxTapeSize * sizeof(float_tape_node)));
		t.node_index_ = 0;
		t.nNodes_ = maxTapeSize;
		aad_float_shifted_sabr_params params;
		params.forward_ = AddFloatToTape(2.5e-2, &t);
		params.timeToExpiry_ = AddFloatToTape(10, &t);
		params.alpha_ = AddFloatToTape(0.28, &t);
		params.beta_ = AddFloatToTape(0.36, &t);
		params.rho_ = AddFloatToTape(0.05, &t);
		params.zeta_ = AddFloatToTape(0.03, &t);
		params.sigma0_ = AddFloatToTape(0.019, &t);
		aad_float thisStrike = AddFloatToTape(strikeStart + i * strikeStep, &t);
		aad_float impliedVol = AadFloatShiftedSabr(thisStrike, &params);
		myFile << thisStrike.value_ << "," << impliedVol.value_ << "\n";
		free(t.nodes_);
	}
	myFile.close();
}

