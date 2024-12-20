#include <math.h>
#include <stdio.h>
#ifdef __APPLE__
#include "imgui/backends/imgui_impl_metal.h"
#include "implot/implot.h"
#include "imgui/examples/example_apple_metal/options.h"
#include "imgui/examples/example_apple_metal/yield_curve.h"
#else
#define IMGUI_DEFINE_MATH_OPERATORS
#include "imgui/imgui_draw.cpp"
#include "imgui/imgui_tables.cpp"
#include "imgui/imgui_widgets.cpp"
#include "imgui/imgui.cpp"
#include "implot/implot_items.cpp"
#include "implot/implot.cpp"
#include "yield_curve.cpp"
#include "options.cpp"
#endif


struct DisplayData {

	uint number_of_fras = 8;
	uint number_of_swaps = 5;
	uint number_of_plot_points = 200;
	void *yield_curve = {};
	Arena *arena;

	Swap *swap_list = nullptr;
	float *plot_maturities = nullptr;
	float *plot_short_rates = nullptr;
	float *plot_overnight_forwards = nullptr;
	float *plot_discount_factors = nullptr;
	float *fra_times_to_maturity = nullptr;
	float *forward_rates = nullptr;
	float *swap_fixed_payment_schedule = nullptr;
	float *swap_rates = nullptr;
	float *forward_starting_swap_rates = nullptr;

	uint number_of_cms = 5;
	uint number_of_strikes_per_smile = 5;
	float *cms_replication_vols = nullptr;
	float *cms_replication_strikes = nullptr;
	float *strikes = nullptr;
	float *normal_vols = nullptr;
	float *times_to_maturity = nullptr;
	float *cms_values = nullptr;
	float cms_maturity = 10;

	float slider_width = 200;
	enum { SHORT_RATE, OVERNIGHT_FORWARD } curve_type = SHORT_RATE;

	float f = 0.0f;
	int counter = 0;
	ShiftedSabrParamGrid shifted_sabr_params = {};
	bool *smile_is_selected_for_plotting = nullptr;
	uint expiry_idx = 0;
	uint underlying_idx = 0;
	OptionModelParametersShiftedSabrFloat *interpolated_shifted_sabr_parameters = nullptr;

	LGM1FTermStructureFloat lgm_term_structure = { 0 };

};

void displayDataInit(
	DisplayData *data,
	uint number_of_fras,
	uint number_of_swaps,
	uint number_of_plot_points,
	uint number_of_cms,
	uint number_of_strikes_per_smile,
	uint number_of_expiries,
	uint number_of_underlyings,
	float *times_to_expiry,
	float *underlying_maturity_length,
	float cms_maturity,
	float slider_width,
	Arena *arena
) {
	if (!arena) {
		arena = ArenaAllocateDefault();
	}
	data->arena = arena;
	data->number_of_fras = number_of_fras;
	data->number_of_swaps = number_of_swaps;
	data->number_of_plot_points = number_of_plot_points;
	data->number_of_strikes_per_smile = number_of_strikes_per_smile;
	data->cms_maturity = cms_maturity;
	data->slider_width = slider_width;
	data->fra_times_to_maturity = (float *)ArenaGetMemory(number_of_fras * sizeof(float), arena);
	data->forward_rates = (float *)ArenaGetMemory(number_of_fras * sizeof(float), arena);
	data->swap_fixed_payment_schedule = (float *)ArenaGetMemory(number_of_swaps * sizeof(float), arena);
	data->swap_rates = (float *)ArenaGetMemory(number_of_swaps * sizeof(float), arena);
	data->forward_starting_swap_rates = (float *)ArenaGetMemory(number_of_cms * sizeof(float), arena);
	data->plot_maturities = (float *)ArenaGetMemory(number_of_plot_points * sizeof(float), arena);
	data->plot_short_rates = (float *)ArenaGetMemory(number_of_plot_points * sizeof(float), arena);
	data->plot_discount_factors = (float *)ArenaGetMemory(number_of_plot_points * sizeof(float), arena);
	data->plot_overnight_forwards = (float *)ArenaGetMemory(number_of_plot_points * sizeof(float), arena);
	data->times_to_maturity = (float *)ArenaGetMemory(number_of_cms * sizeof(float), arena);
	data->cms_values = (float *)ArenaGetMemory(number_of_cms * sizeof(float), arena);
	data->cms_replication_vols = (float *)ArenaGetMemory(number_of_strikes_per_smile * sizeof(float), arena);
	data->cms_replication_strikes = (float *)ArenaGetMemory(number_of_strikes_per_smile * sizeof(float), arena);
	uint number_of_smiles = number_of_expiries * number_of_underlyings;
	data->strikes = (float *)ArenaGetMemory(number_of_strikes_per_smile * number_of_smiles * sizeof(float), arena);
	data->normal_vols = (float *)ArenaGetMemory(number_of_strikes_per_smile * number_of_smiles * sizeof(float), arena);
	data->smile_is_selected_for_plotting = (bool *)ArenaGetMemory(sizeof(bool) * number_of_smiles, arena);
	data->shifted_sabr_params.parameters = (OptionModelParametersShiftedSabrFloat *)ArenaGetMemory(sizeof(OptionModelParametersShiftedSabrFloat) * number_of_smiles, arena);
	data->interpolated_shifted_sabr_parameters = (OptionModelParametersShiftedSabrFloat *)ArenaGetMemory(sizeof(OptionModelParametersShiftedSabrFloat), arena);
	data->shifted_sabr_params.number_of_expiries = number_of_expiries;
	data->shifted_sabr_params.number_of_underlyings = number_of_underlyings;
	data->shifted_sabr_params.times_to_expiry = times_to_expiry;
	data->shifted_sabr_params.underlying_maturity_length = underlying_maturity_length;

	const float starting_sigma_0 = 0.01;
	const float starting_alpha = 1;
	const float starting_beta = 0.15;
	const float starting_rho = 0.2;
	const float starting_zeta = 0.03;

	OptionModelParametersShiftedSabrFloat *this_parameter_set = data->shifted_sabr_params.parameters;
	bool *this_smile_is_selected_for_plotting = data->smile_is_selected_for_plotting;
	for (uint i = 0; i < number_of_expiries; ++i) {
		for (uint j = 0; j < number_of_underlyings; ++j) {
			this_parameter_set->model_header.model_type = OptionModelType::SHIFTED_SABR;
			this_parameter_set->model_header.vol_type = OptionVolType::NORMAL;
			this_parameter_set->sigma_0 = starting_sigma_0;
			this_parameter_set->alpha = starting_alpha;
			this_parameter_set->beta = starting_beta;
			this_parameter_set->rho = starting_rho;
			this_parameter_set->zeta = starting_zeta;
			++this_parameter_set;
			*this_smile_is_selected_for_plotting++ = false;
		}
	}


	for (uint i = 0; i < number_of_fras; ++i) {
		data->fra_times_to_maturity[i] = 0.25 * (i + 1);
		data->forward_rates[i] = 0.05 + i * 0.005;
	}
	data->swap_fixed_payment_schedule[0] = 5;
	data->swap_fixed_payment_schedule[1] = 7;
	data->swap_fixed_payment_schedule[2] = 10;
	data->swap_fixed_payment_schedule[3] = 15;
	data->swap_fixed_payment_schedule[4] = 30;

	for (uint i = 0; i < number_of_swaps; ++i) {
		data->swap_rates[i] = 0.05 + i * 0.01;
	}

	for (uint i = 0; i < number_of_cms; ++i) {
		data->times_to_maturity[i] = (1 + static_cast<float>(i)) * 5;
	}

	uint number_of_stripping_instruments = number_of_swaps + number_of_fras;
	data->yield_curve = ArenaGetMemory(sizeof(YieldCurveShortRateFloat), arena);
	YieldCurveInit(
		data->yield_curve,
		YieldCurveType::SHORT_RATE,
		Precision::FLOAT,
		number_of_stripping_instruments,
		arena
	);

	data->swap_list = (Swap *)ArenaGetMemory(number_of_swaps * sizeof(Swap), arena);
	for (uint i = 0; i < number_of_swaps; ++i) {
		uint number_fixed_payments = (uint)(*(data->swap_fixed_payment_schedule + i));
		uint number_floating_payments = 10;
		SwapScheduleF *schedule = (SwapScheduleF *)ArenaGetMemory(sizeof(SwapScheduleF), arena);
		schedule->fixed_payment_times = (float *)ArenaGetMemory(number_fixed_payments * sizeof(float), arena);
		schedule->floating_payment_times = (float *)ArenaGetMemory(number_floating_payments * sizeof(float), arena);
		data->swap_list[i].schedule = schedule;
		SwapScheduleInitF(
			schedule,
			number_fixed_payments,
			number_floating_payments,
			data->swap_rates[i]
		);
	}

	uint term_structure_size = 8;
	data->lgm_term_structure.number_of_points = term_structure_size;
	data->lgm_term_structure.tenors = (float *)ArenaGetMemory(term_structure_size * sizeof(float), arena);
	data->lgm_term_structure.tenors[0] = 0;
	data->lgm_term_structure.tenors[1] = 1;
	data->lgm_term_structure.tenors[2] = 2;
	data->lgm_term_structure.tenors[3] = 5;
	data->lgm_term_structure.tenors[4] = 10;
	data->lgm_term_structure.tenors[5] = 15;
	data->lgm_term_structure.tenors[6] = 20;
	data->lgm_term_structure.tenors[7] = 30;
	data->lgm_term_structure.lambda_term_structure = (PiecewiseFunctionFloat *)ArenaGetMemory(sizeof(PiecewiseFunctionFloat), arena);
	data->lgm_term_structure.mu_term_structure = (PiecewiseFunctionFloat *)ArenaGetMemory(sizeof(PiecewiseFunctionFloat), arena);
	data->lgm_term_structure.sigma_term_structure = (PiecewiseFunctionFloat *)ArenaGetMemory(sizeof(PiecewiseFunctionFloat), arena);
	PiecewiseFunctionFloatAlloc(data->lgm_term_structure.lambda_term_structure, term_structure_size, arena);
	PiecewiseFunctionFloatAlloc(data->lgm_term_structure.mu_term_structure, term_structure_size, arena);
	PiecewiseFunctionFloatAlloc(data->lgm_term_structure.sigma_term_structure, term_structure_size, arena);
	memcpy(data->lgm_term_structure.lambda_term_structure->time_term_structure, data->lgm_term_structure.tenors, term_structure_size * sizeof(float));
	memcpy(data->lgm_term_structure.mu_term_structure->time_term_structure, data->lgm_term_structure.tenors, term_structure_size * sizeof(float));
	memcpy(data->lgm_term_structure.sigma_term_structure->time_term_structure, data->lgm_term_structure.tenors, term_structure_size * sizeof(float));
	data->lgm_term_structure.lambdas = (float **)ArenaGetMemory(term_structure_size * sizeof(float *), arena);
	data->lgm_term_structure.sigmas = (float **)ArenaGetMemory(term_structure_size * sizeof(float *), arena);
	data->lgm_term_structure.mus = (float **)ArenaGetMemory(term_structure_size * sizeof(float *), arena);
	for (uint i = 0; i < term_structure_size; ++i) {
		data->lgm_term_structure.lambdas[i] = data->lgm_term_structure.lambda_term_structure->piecewise_functions[i].polynomials[0].coefficients;
		data->lgm_term_structure.mus[i] = data->lgm_term_structure.mu_term_structure->piecewise_functions[i].polynomials[0].coefficients;
		data->lgm_term_structure.sigmas[i] = data->lgm_term_structure.sigma_term_structure->piecewise_functions[i].polynomials[0].coefficients;
	}

}

float GetStrikeValue(uint strike_idx, uint number_of_strikes_per_smile, float forward) {
	float far_left_strike = forward - 0.03f;
	float far_right_strike = forward + 0.1f;
	return far_left_strike + (float)strike_idx / (float)number_of_strikes_per_smile * far_right_strike;
}

void imGuiRenderLoop(DisplayData *data) {
	ImGuiIO &io = ImGui::GetIO();
	ImGui::Begin("Yield Curve Viewer");

	if (ImPlot::BeginPlot("Yield Curve USD")) {
		ImPlot::SetupAxesLimits(0, 35, -0.05, 1);
		ImPlot::PlotLine("Short Rate", data->plot_maturities, data->plot_short_rates, data->number_of_plot_points);
		ImPlot::PlotLine("Discount Factors", data->plot_maturities, data->plot_discount_factors, data->number_of_plot_points);
		ImPlot::PlotLine("Overnight Forwards", data->plot_maturities, data->plot_overnight_forwards, data->number_of_plot_points);
		ImPlot::EndPlot();
	}
	ImGui::Text("Application average %.3f ms/frame (%.1f FPS)", 1000.0f / io.Framerate, io.Framerate);
	ImGui::End();

	ImGui::Begin("Input yield curve data");

	static bool model_is_selected[2] = { false, false };
	if (ImGui::TreeNode("Yield Curve Model Type")) {
		if (ImGui::Selectable("Short Rate", model_is_selected[0])) {
			if (!ImGui::GetIO().KeyCtrl) {
				memset(model_is_selected, 0, sizeof(model_is_selected));
				model_is_selected[0] = true;
			}
		}

		if (ImGui::Selectable("Overnight Forward", model_is_selected[1])) {
			if (!ImGui::GetIO().KeyCtrl) {
				memset(model_is_selected, 0, sizeof(model_is_selected));
				model_is_selected[1] = true;
			}
		}
		ImGui::TreePop();
	}

	if (model_is_selected[0]) {
		data->curve_type = DisplayData::SHORT_RATE;
	} else if (model_is_selected[1]) {
		data->curve_type = DisplayData::OVERNIGHT_FORWARD;
	}

	for (uint fra_idx = 0; fra_idx < data->number_of_fras; ++fra_idx) {

		char fra_end_label[32];
		char forward_rate_label[32];
		snprintf(fra_end_label, (int)32, "Fra %i end", fra_idx + 1);
		snprintf(forward_rate_label, (int)32, "Forward %i", fra_idx + 1);

		ImGui::PushItemWidth(data->slider_width);
		ImGui::SliderFloat(fra_end_label, data->fra_times_to_maturity + fra_idx, 0, *(data->swap_fixed_payment_schedule));
		ImGui::SameLine();
		ImGui::SliderFloat(forward_rate_label, data->forward_rates + fra_idx, -0.02, 0.2);
	}
	int minimal_swap_maturity = (int)(data->fra_times_to_maturity[data->number_of_fras - 1]);
	float *this_swap_maturity = data->swap_fixed_payment_schedule;
	float *this_swap_rate = data->swap_rates;

	for (uint swap_idx = 0; swap_idx < data->number_of_swaps; ++swap_idx) {
		char swap_end_label[32];
		char swap_rate_label[32];
		snprintf(swap_end_label, (int)32, "Swap %i end", swap_idx + 1);
		snprintf(swap_rate_label, (int)32, "Swap %i", swap_idx + 1);

		ImGui::PushItemWidth(data->slider_width);
		ImGui::SliderFloat(swap_end_label, this_swap_maturity, minimal_swap_maturity, 50);
		int number_fixed_payments = (int)(*this_swap_maturity);
		ImGui::SameLine();
		ImGui::SliderFloat(swap_rate_label, this_swap_rate, -0.02, 0.2);
		SwapScheduleF *schedule = (SwapScheduleF *)data->swap_list[swap_idx].schedule;
		schedule->fixed_rate = *(data->swap_rates + swap_idx);
		minimal_swap_maturity = data->swap_fixed_payment_schedule[swap_idx];
		SwapScheduleInitF(
			data->swap_list[swap_idx].schedule,
			number_fixed_payments,
			10,
			*this_swap_rate++
		);
		++this_swap_maturity;
	}

	switch (data->curve_type) {
	case DisplayData::SHORT_RATE: {
		YieldCurveShortRateFloat *yield_curve = (YieldCurveShortRateFloat *)data->yield_curve;
		YieldCurveShortRateStripFras(
			data->yield_curve,
			data->forward_rates,
			data->fra_times_to_maturity,
			data->number_of_fras
		);
		YieldCurveShortRateStripSwaps(
			data->yield_curve,
			data->swap_list,
			data->swap_rates,
			data->number_of_swaps,
			data->number_of_fras
		);
		float *this_short_rate = yield_curve->short_rates;
		float *this_yield_curve_maturity_time = yield_curve->interpolation_times;
		float this_maturity_time = 0;
		float first_maturity_time = 0;
		float last_maturity_time = yield_curve->interpolation_times[yield_curve->number_of_points - 1];
		float plot_maturity_step = (last_maturity_time - first_maturity_time) / ((float)(data->number_of_plot_points));
		float log_discount_factor = 0;
		for (uint i = 0; i < data->number_of_plot_points; ++i) {
			if (this_maturity_time > *this_yield_curve_maturity_time) {
				++this_yield_curve_maturity_time;
				++this_short_rate;
			}
			log_discount_factor -= *this_short_rate * plot_maturity_step;
			data->plot_maturities[i] = this_maturity_time;
			data->plot_short_rates[i] = *this_short_rate;
			data->plot_discount_factors[i] = expf(log_discount_factor);
			data->plot_overnight_forwards[i] = 365 * (expf(*this_short_rate / 365) - 1);
			this_maturity_time += plot_maturity_step;
		}
	} break;
	case DisplayData::OVERNIGHT_FORWARD: {
		YieldCurveOvernightForwardFloat *yield_curve = (YieldCurveOvernightForwardFloat *)data->yield_curve;
		YieldCurveOvernightForwardStripFras(
			yield_curve,
			data->forward_rates,
			data->fra_times_to_maturity,
			data->number_of_fras
		);
		YieldCurveOvernightForwardStripSwaps(
			yield_curve,
			data->swap_list,
			data->swap_rates,
			data->number_of_swaps,
			data->number_of_fras
		);
		float *this_overnight_forward = yield_curve->overnight_forwards;
		float *this_yield_curve_maturity_time = yield_curve->interpolation_times;
		float this_maturity_time = 0;
		float first_maturity_time = 0;
		float last_maturity_time = yield_curve->interpolation_times[yield_curve->number_of_points - 1];
		float plot_maturity_step = (last_maturity_time - first_maturity_time) / ((float)(data->number_of_plot_points));
		float discount_factor = 1;
		for (uint i = 0; i < data->number_of_plot_points; ++i) {
			this_maturity_time += plot_maturity_step;
			if (this_maturity_time > *this_yield_curve_maturity_time) {
				++this_yield_curve_maturity_time;
				++this_overnight_forward;
			}
			discount_factor /= (1 + *this_overnight_forward * plot_maturity_step);
			data->plot_maturities[i] = this_maturity_time;
			data->plot_short_rates[i] = 365 * logf(1 + *this_overnight_forward / 365);
			data->plot_discount_factors[i] = discount_factor;
			data->plot_overnight_forwards[i] = *this_overnight_forward;
		}
	} break;
	}
	ImGui::End();

	ImGui::Begin("Smile input data");

	uint number_of_expiries = data->shifted_sabr_params.number_of_expiries;
	uint number_of_underlyings = data->shifted_sabr_params.number_of_underlyings;
	float *expiries = data->shifted_sabr_params.times_to_expiry;
	float *underlying_tenors = data->shifted_sabr_params.underlying_maturity_length;
	OptionModelParametersShiftedSabrFloat *model = data->shifted_sabr_params.parameters;
	for (uint i = 0; i < number_of_expiries; ++i) {
		for (uint j = 0; j < number_of_underlyings; ++j) {
			float *this_strike = data->strikes + i * data->number_of_strikes_per_smile;
			float *this_vol = data->normal_vols + i * data->number_of_strikes_per_smile;
			float frequency = underlying_tenors[j] < 1 ? 1 / underlying_tenors[j] : 1;
			float forward = SwapRateYieldCurveF(data->yield_curve, expiries[i], underlying_tenors[j], 1);
			for (uint k = 0; k < data->number_of_strikes_per_smile; ++k) {
				*this_strike = GetStrikeValue(k, data->number_of_strikes_per_smile, forward);
				*this_vol++ = ShiftedSabrImpliedNormalVol(
					forward,
					*this_strike++,
					expiries[i],
					model->sigma_0,
					model->alpha,
					model->beta,
					model->rho,
					model->zeta
				);
			}
			++model;
		}
	}

	static uint selected_expiry_idx;
	static uint selected_underlying_idx;

	static ImGuiTableFlags flags = ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg;
	char underlying_label[32];
	char expiry_label[32];
	char button_label[32];
	uint number_of_rows = data->shifted_sabr_params.number_of_expiries;
	uint number_of_columns = data->shifted_sabr_params.number_of_underlyings;

	OptionModelParametersShiftedSabrFloat *shifted_sabr_params = data->shifted_sabr_params.parameters + selected_expiry_idx * number_of_columns + selected_underlying_idx;
	if (ImGui::BeginTable("OptionGrid", data->shifted_sabr_params.number_of_underlyings, flags)) {
		for (uint col_idx = 0; col_idx < number_of_columns; ++col_idx) {
			snprintf(underlying_label, (uint)32, "%.2fY", data->shifted_sabr_params.underlying_maturity_length[col_idx]);
			ImGui::TableSetupColumn(underlying_label);
		}
		ImGui::TableHeadersRow();

		bool *smile_is_selected_for_plotting = data->smile_is_selected_for_plotting;
		for (uint row_idx = 0; row_idx < number_of_rows; ++row_idx) {
			ImGui::TableNextRow();
			for (uint col_idx = 0; col_idx < number_of_columns; ++col_idx) {
				ImGui::TableSetColumnIndex(col_idx);
				snprintf(
					button_label,
					(uint)32,
					"%.2fY%.2fY",
					data->shifted_sabr_params.times_to_expiry[row_idx],
					data->shifted_sabr_params.underlying_maturity_length[col_idx]
				);
				if (ImGui::Button(button_label, ImVec2(-FLT_MIN, 0.0f))) {
					selected_expiry_idx = row_idx;
					selected_underlying_idx = col_idx;
					*smile_is_selected_for_plotting = !*smile_is_selected_for_plotting;
				}
				++smile_is_selected_for_plotting;
			}
		}
		ImGui::EndTable();
	}


	snprintf(expiry_label, (uint)32, "%.2fY%.2fY", data->shifted_sabr_params.times_to_expiry[selected_expiry_idx], data->shifted_sabr_params.underlying_maturity_length[selected_underlying_idx]);
	ImGui::Text(expiry_label);
	char parameter_label[32];
	snprintf(parameter_label, (uint)32, "Sigma0");
	ImGui::PushItemWidth(data->slider_width);
	ImGui::SliderFloat(parameter_label, &(shifted_sabr_params->sigma_0), 0.0001, 0.1);

	snprintf(parameter_label, (uint)32, "Alpha");
	ImGui::PushItemWidth(data->slider_width);
	ImGui::SliderFloat(parameter_label, &(shifted_sabr_params->alpha), 0.0001, 1.5);

	snprintf(parameter_label, (uint)32, "Beta");
	ImGui::PushItemWidth(data->slider_width);
	ImGui::SliderFloat(parameter_label, &(shifted_sabr_params->beta), 0.0001, 0.99);

	snprintf(parameter_label, (uint)32, "Rho");
	ImGui::PushItemWidth(data->slider_width);
	ImGui::SliderFloat(parameter_label, &(shifted_sabr_params->rho), -0.8, 0.8);

	snprintf(parameter_label, (uint)32, "Zeta");
	ImGui::PushItemWidth(data->slider_width);
	ImGui::SliderFloat(parameter_label, &(shifted_sabr_params->zeta), 0.0001, 0.1);

	OptionModelParametersShiftedSabrFloat *these_parameters = data->shifted_sabr_params.parameters;
	uint number_of_vols_per_row = number_of_underlyings * data->number_of_strikes_per_smile;
	float *this_expiry = data->shifted_sabr_params.times_to_expiry;
	for (uint i = 0; i < number_of_expiries; ++i) {
		float *this_underlying_tenor = data->shifted_sabr_params.underlying_maturity_length;
		for (uint j = 0; j < number_of_underlyings; ++j) {
			float *this_smile = data->normal_vols + i * number_of_vols_per_row + j * data->number_of_strikes_per_smile;
			float *this_strike = data->strikes + i * number_of_vols_per_row + j * data->number_of_strikes_per_smile;
			float frequency = *this_underlying_tenor < 1 ? 1 / *this_underlying_tenor : 1;
			float forward = SwapRateYieldCurveF(data->yield_curve, *this_expiry, *this_underlying_tenor, frequency);
			for (uint k = 0; k < data->number_of_strikes_per_smile; ++k) {
				*this_strike = GetStrikeValue(k, data->number_of_strikes_per_smile, forward);
				*this_smile++ = ShiftedSabrImpliedNormalVol(
					forward,
					*this_strike++,
					*this_expiry,
					these_parameters->sigma_0,
					these_parameters->alpha,
					these_parameters->beta,
					these_parameters->rho,
					these_parameters->zeta
				);
			}
			++these_parameters;
			++this_underlying_tenor;
		}
		++this_expiry;
	}

	float *this_time_to_maturity = data->times_to_maturity;
	float *this_cms_value = data->cms_values;
	float *this_forward_swap_rate = data->forward_starting_swap_rates;
	OptionModelParametersShiftedSabrFloat *interpolated_parameters = (OptionModelParametersShiftedSabrFloat *)data->interpolated_shifted_sabr_parameters;
	for (uint i = 0; i < data->number_of_cms; ++i) {
		float df;
		DiscountFactorSpotFromYieldCurve(&df, data->yield_curve, this_time_to_maturity);
		float forward = SwapRateYieldCurveF(data->yield_curve, *this_time_to_maturity, data->cms_maturity, 1);
		InterpolateShiftedSabrParametersFloat(interpolated_parameters, *this_time_to_maturity, data->cms_maturity, &(data->shifted_sabr_params));
		float *this_vol = data->cms_replication_vols;
		float *this_strike = data->cms_replication_strikes;
		for (uint k = 0; k < data->number_of_strikes_per_smile; ++k) {
			*this_strike = GetStrikeValue(k, data->number_of_strikes_per_smile, forward);
			*this_vol++ = ShiftedSabrImpliedNormalVol(
				forward,
				*this_strike++,
				*this_time_to_maturity,
				interpolated_parameters->sigma_0,
				interpolated_parameters->alpha,
				interpolated_parameters->beta,
				interpolated_parameters->rho,
				interpolated_parameters->zeta
			);
		}
		*this_forward_swap_rate++ = forward;
		*this_cms_value++ = CmsReplicationForwardPremium(
			forward,
			data->cms_replication_vols,
			data->cms_replication_strikes,
			data->number_of_strikes_per_smile,
			*this_time_to_maturity++,
			1,
			data->cms_maturity
		);
	}
	ImGui::End();

	ImGui::Begin("Options model");
	if (ImPlot::BeginPlot("Smiles")) {
		char maturity_label[32];
		float *this_time_to_expiry = data->shifted_sabr_params.times_to_expiry;
		bool *this_smile_is_selected_for_plotting = data->smile_is_selected_for_plotting;
		for (uint i = 0; i < data->shifted_sabr_params.number_of_expiries; ++i) {
			float *this_underlying_tenor = data->shifted_sabr_params.underlying_maturity_length;
			for (uint j = 0; j < data->shifted_sabr_params.number_of_underlyings; ++j) {
				if (*this_smile_is_selected_for_plotting++) {
					float *strikes = data->strikes + i * data->shifted_sabr_params.number_of_underlyings * data->number_of_strikes_per_smile + j * data->number_of_strikes_per_smile;
					float *vols = data->normal_vols + i * data->shifted_sabr_params.number_of_underlyings * data->number_of_strikes_per_smile + j * data->number_of_strikes_per_smile;
					snprintf(maturity_label, (uint)32, "%.2fY%.2fY", *this_time_to_expiry, *this_underlying_tenor);
					ImPlot::PlotLine(
						maturity_label,
						strikes,
						vols,
						data->number_of_strikes_per_smile
					);
					float frequency = *this_underlying_tenor < 1 ? 1 / *this_underlying_tenor : 1;
					float forward = SwapRateYieldCurveF(data->yield_curve, *this_time_to_expiry, *this_underlying_tenor, frequency);
					InterpolateShiftedSabrParametersFloat(
						interpolated_parameters,
						*this_time_to_expiry,
						*this_underlying_tenor,
						&(data->shifted_sabr_params)
					);
					float atm_vol = ShiftedSabrImpliedNormalVol(
						forward,
						forward,
						*this_time_to_expiry,
						interpolated_parameters->sigma_0,
						interpolated_parameters->alpha,
						interpolated_parameters->beta,
						interpolated_parameters->rho,
						interpolated_parameters->zeta
					);
					ImPlot::SetNextMarkerStyle(ImPlotMarker_Cross, -1, ImVec4(0, 0, 0, -1), 10);
					ImPlot::PlotLine(maturity_label, &forward, &atm_vol, 1);
				}
				++this_underlying_tenor;
			}
			++this_time_to_expiry;
		}
		ImPlot::EndPlot();
	}
	ImGui::End();

	ImGui::Begin("Forward looking");
	char cms_label[32];
	char forward_swap_label[32];
	snprintf(cms_label, (uint)32, "%i Y CMS", (int)(data->cms_maturity));
	snprintf(forward_swap_label, (uint)32, "%i Y Forward swap", (int)(data->cms_maturity));
	if (ImPlot::BeginPlot("Forward swaps and CMS")) {
		ImPlot::PlotLine(cms_label, data->times_to_maturity, data->cms_values, data->number_of_cms);
		ImPlot::PlotLine(forward_swap_label, data->times_to_maturity, data->forward_starting_swap_rates, data->number_of_cms);
		ImPlot::EndPlot();
	}
	ImGui::End();

	ImGui::Begin("LGM Term Structure");
	ImGui::PushItemWidth(data->slider_width);
	static float constant_lambda_value;
	const float min_lambda_value = 0.001;
	if (constant_lambda_value < min_lambda_value) {
		constant_lambda_value = min_lambda_value;
	}
	ImGui::SliderFloat("Lambda", &constant_lambda_value, min_lambda_value, 5);
	LGM1FSetConstantLambda(&(data->lgm_term_structure), constant_lambda_value);
	if (ImPlot::BeginPlot("LGM Term Structure")) {
		const uint term_struct_size = 16;
		static float lambdas[term_struct_size];
		static float sigmas[term_struct_size];
		static float mus[term_struct_size];
		CapletVolParams target_caplets[term_struct_size];
		float caplet_expiries[term_struct_size] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 15, 17, 20, 25, 29 };
		uint number_of_caplets = data->lgm_term_structure.number_of_points - 2;
		assert (term_struct_size > number_of_caplets);
		float *tenors = data->lgm_term_structure.tenors + 1;
		for (uint i = 0; i < number_of_caplets; ++i) {
			float forward = SwapRateYieldCurveF(data->yield_curve, tenors[i], 1, 1);
			InterpolateShiftedSabrParametersFloat(
				interpolated_parameters,
				tenors[i],
				1,
				&(data->shifted_sabr_params)
			);
			target_caplets[i].normal_vol = ShiftedSabrImpliedNormalVol(
				forward,
				forward,
				tenors[i],
				interpolated_parameters->sigma_0,
				interpolated_parameters->alpha,
				interpolated_parameters->beta,
				interpolated_parameters->rho,
				interpolated_parameters->zeta
			);
			target_caplets[i].coverage = 1;
			target_caplets[i].forward = forward;
			target_caplets[i].time_to_expiry = tenors[i];
		}

		LGM1FFitSigmaTermStructureFromTermCaplets(
			&(data->lgm_term_structure),
			target_caplets,
			1e-5f,
			20
		);
		for (uint i = 0; i < data->lgm_term_structure.number_of_points; ++i) {
			lambdas[i] = *(data->lgm_term_structure.lambdas[i]);
			sigmas[i] = *(data->lgm_term_structure.sigmas[i]);
			mus[i] = *(data->lgm_term_structure.mus[i]);
		}
		ImPlot::PlotLine("Lambda", data->lgm_term_structure.tenors, lambdas, data->lgm_term_structure.number_of_points);
		ImPlot::PlotLine("Mu", data->lgm_term_structure.tenors, mus, data->lgm_term_structure.number_of_points);
		ImPlot::PlotLine("Sigma", data->lgm_term_structure.tenors, sigmas, data->lgm_term_structure.number_of_points);
		ImPlot::EndPlot();
	}
	ImGui::End();

}


