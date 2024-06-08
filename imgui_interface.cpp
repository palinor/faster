#include <math.h>
#include <stdio.h>
#ifdef __APPLE__
#include "imgui/backends/imgui_impl_metal.h"
#include "implot/implot.h"
#include "imgui/examples/example_apple_metal/options.h"
#include "imgui/examples/example_apple_metal/yield_curve.h"
#else
#define IMGUI_DEFINE_MATH_OPERATORS
#include "imgui_impl_win32.cpp"
#include "imgui_impl_dx12.cpp"
#include "imgui_draw.cpp"
#include "imgui_tables.cpp"
#include "imgui_widgets.cpp"
#include "imgui.cpp"
#include "implot_items.cpp"
#include "implot.cpp"
#include "yield_curve.cpp"
#include "options.cpp"
#endif


struct DisplayData {

	size_t number_of_fras = 8;
	size_t number_of_swaps = 5;
	size_t number_of_plot_points = 200;
	void *yield_curve = {};

	YieldCurveSwapF *swap_list = nullptr;
	float *plot_maturities = nullptr;
	float *plot_short_rates = nullptr;
	float *plot_overnight_forwards = nullptr;
	float *plot_discount_factors = nullptr;
	float *fra_times_to_maturity = nullptr;
	float *forward_rates = nullptr;
	float *swap_fixed_payment_schedule = nullptr;
	float *swap_rates = nullptr;

	size_t number_of_cms = 5;
	size_t number_of_strikes = 5;
	OptionModelParametersShiftedSabrFloat *shifted_sabr_params = nullptr;
	float *strikes = nullptr;
	float *normal_vols = nullptr;
	float *times_to_maturity = nullptr;
	float *cms_values = nullptr;
	float cms_maturity = 10;

	float slider_width = 200;
	enum { SHORT_RATE, OVERNIGHT_FORWARD } curve_type = SHORT_RATE;

	float f = 0.0f;
	int counter = 0;
};

void displayDataInit(
	DisplayData *data,
	size_t number_of_fras,
	size_t number_of_swaps,
	size_t number_of_plot_points,
	size_t number_of_cms,
	size_t number_of_strikes,
	float cms_maturity,
	float slider_width
) {
	data->number_of_fras = number_of_fras;
	data->number_of_swaps = number_of_swaps;
	data->number_of_plot_points = number_of_plot_points;
	data->number_of_strikes = number_of_strikes;
	data->cms_maturity = cms_maturity;
	if (!data->fra_times_to_maturity) data->fra_times_to_maturity = (float *)malloc(number_of_fras * sizeof(float));
	if (!data->forward_rates) data->forward_rates = (float *)malloc(number_of_fras * sizeof(float));
	if (!data->swap_fixed_payment_schedule) data->swap_fixed_payment_schedule = (float *)calloc(number_of_swaps, sizeof(float));
	if (!data->swap_rates) data->swap_rates = (float *)malloc(number_of_swaps * sizeof(float));
	if (!data->plot_maturities) data->plot_maturities = (float *)calloc(number_of_plot_points, sizeof(float));
	if (!data->plot_short_rates) data->plot_short_rates = (float *)calloc(number_of_plot_points, sizeof(float));
	if (!data->plot_discount_factors) data->plot_discount_factors = (float *)calloc(number_of_plot_points, sizeof(float));
	if (!data->plot_overnight_forwards) data->plot_overnight_forwards = (float *)calloc(number_of_plot_points, sizeof(float));
	if (!data->strikes) data->strikes = (float *)malloc(number_of_strikes * sizeof(float));
	if (!data->times_to_maturity) data->times_to_maturity = (float *)malloc(number_of_strikes * sizeof(float));
	if (!data->cms_values) data->cms_values = (float *)malloc(number_of_strikes * sizeof(float));
	if (!data->normal_vols) data->normal_vols = (float *)malloc(number_of_strikes * number_of_cms * sizeof(float));
	if (!data->shifted_sabr_params) data->shifted_sabr_params = (OptionModelParametersShiftedSabrFloat *)malloc(sizeof(OptionModelParametersShiftedSabrFloat));

	data->shifted_sabr_params->sigma_0 = 0.03;
	data->shifted_sabr_params->alpha = 0.22;
	data->shifted_sabr_params->beta = 0.4;
	data->shifted_sabr_params->rho = 0.01;
	data->shifted_sabr_params->zeta = 0.03;


	for (size_t i = 0; i < number_of_fras; ++i) {
		data->fra_times_to_maturity[i] = 0.25 * (i + 1);
		data->forward_rates[i] = 0.05 + i * 0.005;
	}
	data->swap_fixed_payment_schedule[0] = 5;
	data->swap_fixed_payment_schedule[1] = 7;
	data->swap_fixed_payment_schedule[2] = 10;
	data->swap_fixed_payment_schedule[3] = 15;
	data->swap_fixed_payment_schedule[4] = 30;

	for (size_t i = 0; i < number_of_swaps; ++i) {
		data->swap_rates[i] = 0.05 + i * 0.01;
	}

	for (size_t i = 0; i < number_of_cms; ++i) {
		data->times_to_maturity[i] = (1 + static_cast<float>(i)) * 5;
	}

	data->swap_list = (YieldCurveSwapF *)calloc(number_of_swaps, sizeof(YieldCurveSwapF));
	for (size_t i = 0; i < number_of_swaps; ++i) {
		size_t number_fixed_payments = (size_t)(*(data->swap_fixed_payment_schedule + i));
		size_t number_floating_payments = 10;
		YieldCurveSwapFInit(data->swap_list + i, number_fixed_payments, number_floating_payments, data->swap_rates[i]);
	}

	float *this_strike = data->strikes;
	for (size_t i = 0; i < data->number_of_strikes; ++i) {
		*this_strike++ = -0.01 + i * 0.001;
	}
}

void displayDataFree(DisplayData *data) {
	free(data->swap_list);
	free(data->plot_maturities);
	free(data->plot_short_rates);
	free(data->fra_times_to_maturity);
	free(data->forward_rates);
	free(data->swap_fixed_payment_schedule);
	free(data->swap_rates);
	free(data->plot_overnight_forwards);
	free(data->strikes);
	free(data->cms_values);
	free(data->times_to_maturity);
	free(data->normal_vols);
	free(data->shifted_sabr_params);
}


void imGuiRenderLoop(DisplayData *data) {
	ImGuiIO &io = ImGui::GetIO();
	ImGui::Begin("Yield Curve Viewer");

	if (ImPlot::BeginPlot("Yield Curve USD")) {
		ImPlot::SetupAxesLimits(0, 35, -0.05, 1);
		ImPlot::PlotLine("Short Rate", data->plot_maturities, data->plot_short_rates, data->number_of_plot_points);
		ImPlot::PlotLine("Discount Factors", data->plot_maturities, data->plot_discount_factors, data->number_of_plot_points);
		ImPlot::PlotLine("Overnight Forwards", data->plot_maturities, data->plot_overnight_forwards, data->number_of_plot_points);
		ImPlot::PlotLine("CMS Forwards", data->times_to_maturity, data->cms_values, data->number_of_cms);
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
	}
	else if (model_is_selected[1]) {
		data->curve_type = DisplayData::OVERNIGHT_FORWARD;
	}

	if (!(data->yield_curve)) {
		data->yield_curve = calloc(1, sizeof(YieldCurveShortRateFloat));
	}

	for (int fra_idx = 0; fra_idx < data->number_of_fras; ++fra_idx) {

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

	for (int swap_idx = 0; swap_idx < data->number_of_swaps; ++swap_idx) {
		char swap_end_label[32];
		char swap_rate_label[32];
		snprintf(swap_end_label, (int)32, "Swap %i end", swap_idx + 1);
		snprintf(swap_rate_label, (int)32, "Swap %i", swap_idx + 1);

		ImGui::PushItemWidth(data->slider_width);
		ImGui::SliderFloat(swap_end_label, this_swap_maturity, minimal_swap_maturity, 50);
		int number_fixed_payments = (int)(*this_swap_maturity);
		ImGui::SameLine();
		ImGui::SliderFloat(swap_rate_label, this_swap_rate, -0.02, 0.2);
		data->swap_list[swap_idx].swap_rate = *(data->swap_rates + swap_idx);
		minimal_swap_maturity = data->swap_fixed_payment_schedule[swap_idx];
		YieldCurveSwapFInit(data->swap_list + swap_idx, number_fixed_payments, 10, *this_swap_rate++);
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
		float plot_maturity_step_in_days = plot_maturity_step * 365;
		for (size_t i = 0; i < data->number_of_plot_points; ++i) {
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
		for (size_t i = 0; i < data->number_of_plot_points; ++i) {
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

	char parameter_label[32];
	snprintf(parameter_label, (size_t)32, "Sigma0");
	ImGui::PushItemWidth(data->slider_width);
	ImGui::SliderFloat(parameter_label, &(data->shifted_sabr_params->sigma_0), 0.0001, 1.5);

	snprintf(parameter_label, (size_t)32, "Alpha");
	ImGui::PushItemWidth(data->slider_width);
	ImGui::SliderFloat(parameter_label, &(data->shifted_sabr_params->alpha), 0.0001, 1.5);

	snprintf(parameter_label, (size_t)32, "Beta");
	ImGui::PushItemWidth(data->slider_width);
	ImGui::SliderFloat(parameter_label, &(data->shifted_sabr_params->beta), 0.0001, 0.99);

	snprintf(parameter_label, (size_t)32, "Rho");
	ImGui::PushItemWidth(data->slider_width);
	ImGui::SliderFloat(parameter_label, &(data->shifted_sabr_params->rho), -0.8, 0.8);

	snprintf(parameter_label, (size_t)32, "Zeta");
	ImGui::PushItemWidth(data->slider_width);
	ImGui::SliderFloat(parameter_label, &(data->shifted_sabr_params->zeta), 0.0001, 0.1);

	float *this_strike = data->strikes;
	for (size_t i = 0; i < data->number_of_strikes; ++i) {
		*this_strike++ = -0.01 + i * 0.001;
	}
	float *this_time_to_maturity = data->times_to_maturity;
	float *this_cms_value = data->cms_values;
	for (size_t i = 0; i < data->number_of_cms; ++i) {
		float df;
		DiscountFactorSpotFromYieldCurve(&df, data->yield_curve, this_time_to_maturity);
		float forward = SwapRateYieldCurveF(data->yield_curve, *this_time_to_maturity, data->cms_maturity, 1, 1);
		float *this_smile = data->normal_vols + data->number_of_strikes * i;
		float *this_strike = data->strikes;
		for (size_t k = 0; k < data->number_of_strikes; ++k) {
			*(this_smile + k) = ShiftedSabrImpliedNormalVol(
				forward,
				*this_strike++,
				*this_time_to_maturity,
				data->shifted_sabr_params->sigma_0,
				data->shifted_sabr_params->alpha,
				data->shifted_sabr_params->beta,
				data->shifted_sabr_params->rho,
				data->shifted_sabr_params->zeta
			);
		}
		*this_cms_value++ = CmsReplicationForwardPremium(forward, this_smile, data->strikes, data->number_of_strikes, *this_time_to_maturity++, 1, data->cms_maturity);
	}
	ImGui::End();

	ImGui::Begin("Options model");
	if (ImPlot::BeginPlot("Smiles")) {
		char maturity_label[32];
		for (size_t i = 0; i < data->number_of_cms; ++i) {
			snprintf(maturity_label, (size_t)32, "%i Y", static_cast<int>(data->times_to_maturity[i]));
			ImPlot::PlotLine(maturity_label, data->strikes, data->normal_vols + i * data->number_of_strikes, data->number_of_strikes);
		}
		ImPlot::EndPlot();
	}
	ImGui::End();

}


