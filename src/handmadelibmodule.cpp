#define PY_SSIZE_T_CLEAN
#include "Python.h"
#include "aad.cpp"

static PyObject *ShiftedSabr(PyObject *self, PyObject *args) {

	float forward;
	float strike;
	float time_to_expiry;
	float sigma_0;
	float alpha;
	float beta;
	float rho;
	float zeta;

	if (!PyArg_ParseTuple(args, "ffffffff",
		&forward,
		&time_to_expiry,
		&alpha,
		&beta,
		&rho,
		&zeta,
		&sigma_0,
		&strike
	)) {
		return NULL;
	}


	float shifted_forward = forward + zeta;
	float shifted_strike = strike + zeta;
	float moneyness = forward - strike;
	float i_0 = 0;
	float i_1 = 0;
	float z = 0;
	if (beta > 0.99) {
		z = alpha * moneyness / sigma_0;
	} else {
		z = alpha / sigma_0 * (powf(shifted_forward, 1 - beta) - powf(shifted_strike, 1 - beta)) / (1 - beta);
	}
	if (fabsf(moneyness) < 1e-10) {
		i_0 = sigma_0 * powf(shifted_strike, beta);
	} else {
		i_0 = alpha * moneyness / logf(
			(sqrtf(1 - 2 * rho * z + z * z) + z - rho) / (1 - rho)
		);
	}
	i_1 = (((beta - 1) * (beta - 1) - 1) * sigma_0 * sigma_0 / powf((shifted_forward + shifted_strike) / 2, 2 * (1 - beta))
		+ 6 * rho * alpha * sigma_0 * beta / powf((shifted_forward + shifted_strike) / 2, 1 - beta)
		+ 2 * alpha * alpha - 3 * rho * rho * alpha * alpha
		) / 24;
	float result = i_0 * (1 + i_1 * time_to_expiry);
	return PyFloat_FromDouble(result);
}

static PyObject *AadShiftedSabr(PyObject *self, PyObject *args) {
	float forward;
	float time_to_expiry;
	float alpha;
	float beta;
	float rho;
	float zeta;
	float sigma0;
	float strike;

	if (!PyArg_ParseTuple(args, "ffffffff",
		&forward,
		&time_to_expiry,
		&alpha,
		&beta,
		&rho,
		&zeta,
		&sigma0,
		&strike
	)) {
		return NULL;
	}

	TapeFloat tape;
	const size_t starting_tape_size = 1024;
	tape.nodes = reinterpret_cast<TapeNodeFloat *>(malloc(starting_tape_size * sizeof(TapeNodeFloat)));
	tape.current_node_index = 0;
	tape.number_of_nodes = starting_tape_size;
	AadFloat_shifted_sabr_params params;
	AadFloatInit(&(params.forward_), forward, &tape);
	AadFloatInit(&(params.timeToExpiry_), time_to_expiry, &tape);
	AadFloatInit(&(params.alpha_), alpha, &tape);
	AadFloatInit(&(params.beta_), beta, &tape);
	AadFloatInit(&(params.rho_), rho, &tape);
	AadFloatInit(&(params.zeta_), zeta, &tape);
	AadFloatInit(&(params.sigma0_), sigma0, &tape);
	AadFloat thisStrike, impliedVol;
	AadFloatInit(&thisStrike, strike, &tape);
	impliedVol = AadFloatShiftedSabr(thisStrike, &params);
	free(tape.nodes);
	return PyFloat_FromDouble(impliedVol.value);
}

static PyMethodDef handmade_methods[] = {
	{"AadShiftedSabr", AadShiftedSabr, METH_VARARGS, "Test the AAD shifted sabr"},
	{"ShiftedSabr", ShiftedSabr, METH_VARARGS, "Test the regular shifted sabr"},
	{NULL, NULL, 0, NULL}
};
static PyModuleDef module_def = {
	.m_base = PyModuleDef_HEAD_INIT,
	.m_name = "HandmadeLib",
	.m_size = 0,
	.m_methods = handmade_methods
};

PyMODINIT_FUNC
PyInit_handmadelib() {
	return PyModuleDef_Init(&module_def);
}
