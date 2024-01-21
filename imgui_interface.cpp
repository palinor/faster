#include <math.h>
#include <stdio.h>
#include "imgui.h"
#include "implot/implot.h"
#include "yield_curve.h"

struct DisplayData {

	size_t number_of_fras = 8;
	size_t number_of_swaps = 5;
	size_t number_of_plot_points = 200;
	YieldCurveFloat yield_curve = {0};
	
	YieldCurveSwapF *swap_list = nullptr;
	float *plot_maturities = nullptr;
	float *plot_short_rates = nullptr;
	float *plot_discount_factors = nullptr;
	float *fra_times_to_maturity = nullptr;
	float *forward_rates = nullptr;
	float *swap_fixed_payment_schedule = nullptr;
	float *swap_rates = nullptr;

	float slider_width = 200;

	float f = 0.0f;
	int counter = 0;
};

void displayDataInit(
	DisplayData *data,
	size_t number_of_fras,
	size_t number_of_swaps,
	size_t number_of_plot_points,
	float slider_width
) {
	data->number_of_fras = number_of_fras;
	data->number_of_swaps = number_of_swaps;
	data->number_of_plot_points = number_of_plot_points;
	if (!data->fra_times_to_maturity) data->fra_times_to_maturity = (float *)malloc(number_of_fras * sizeof(float));
	if (!data->forward_rates) data->forward_rates = (float *)malloc(number_of_fras * sizeof(float));
	if (!data->swap_fixed_payment_schedule) data->swap_fixed_payment_schedule = (float *)malloc(number_of_swaps * sizeof(float));
	if (!data->swap_rates) data->swap_rates = (float *)malloc(number_of_swaps * sizeof(float));
	if (!data->plot_maturities) data->plot_maturities = (float *)malloc(number_of_plot_points * sizeof(float));
	if (!data->plot_short_rates) data->plot_short_rates = (float *)malloc(number_of_plot_points * sizeof(float));
	if (!data->plot_discount_factors) data->plot_discount_factors = (float *)malloc(number_of_plot_points * sizeof(float));

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

	for (size_t i = 0; i < number_of_plot_points; ++i) {
		data->plot_maturities[i] = i;
		data->plot_short_rates[i] = 0;
		data->plot_discount_factors[i] = 0;
	}

	data->swap_list = (YieldCurveSwapF *)calloc(number_of_swaps, sizeof(YieldCurveSwapF));
	for (size_t i = 0; i < number_of_swaps; ++i) {
		size_t number_fixed_payments = (size_t)(*(data->swap_fixed_payment_schedule + i));
		size_t number_floating_payments = 10;
		yieldCurveSwapFInit(data->swap_list + i, number_fixed_payments, number_floating_payments, data->swap_rates[i]);
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
}


void imGuiRenderLoop(DisplayData *data) {
	ImGuiIO& io = ImGui::GetIO();
	ImGui::Begin("Yield Curve Viewer");

	if (ImPlot::BeginPlot("Yield Curve USD")) {
		ImPlot::PlotLine("Short Rate", data->plot_maturities, data->plot_short_rates, data->number_of_plot_points);
		ImPlot::PlotLine("Discount Factors", data->plot_maturities, data->plot_discount_factors, data->number_of_plot_points);
		ImPlot::EndPlot();
	}
	ImGui::Text("Application average %.3f ms/frame (%.1f FPS)", 1000.0f / io.Framerate, io.Framerate);
	ImGui::End();

	ImGui::Begin("Input data");

	for (size_t fra_idx = 0; fra_idx < data->number_of_fras; ++fra_idx) {

		char fra_end_label[32];
		char forward_rate_label[32];
		snprintf(fra_end_label, 32, "Fra %i end", fra_idx + 1);
		snprintf(forward_rate_label, 32, "Forward %i", fra_idx + 1);

		ImGui::PushItemWidth(data->slider_width);
		ImGui::SliderFloat(fra_end_label, data->fra_times_to_maturity + fra_idx, 0, *(data->swap_fixed_payment_schedule));
		ImGui::SameLine();
		ImGui::SliderFloat(forward_rate_label, data->forward_rates + fra_idx, -0.02, 0.2);
	}
	float minimal_swap_maturity = data->fra_times_to_maturity[data->number_of_fras - 1];

	for (size_t swap_idx = 0; swap_idx < data->number_of_swaps; ++swap_idx) {
		char swap_end_label[32];
		char swap_rate_label[32];
		snprintf(swap_end_label, 32, "Swap %i end", swap_idx + 1);
		snprintf(swap_rate_label, 32, "Swap %i", swap_idx + 1);

		ImGui::PushItemWidth(data->slider_width);
		ImGui::SliderFloat(swap_end_label, data->swap_fixed_payment_schedule + swap_idx, minimal_swap_maturity, 50);
		ImGui::SameLine();
		ImGui::SliderFloat(swap_rate_label, data->swap_rates + swap_idx, -0.02, 0.2);
		data->swap_list[swap_idx].swap_rate = *(data->swap_rates + swap_idx);
		minimal_swap_maturity = data->swap_fixed_payment_schedule[swap_idx];
	}


	yieldCurveFloatStripFras(
		&(data->yield_curve), 
		data->forward_rates, 
		data->fra_times_to_maturity, 
		data->number_of_fras
		);
	yieldCurveFloatStripSwaps(
		&(data->yield_curve), 
		data->swap_list, 
		data->swap_rates, 
		data->number_of_swaps, 
		data->number_of_fras
		);
	float *this_short_rate = data->yield_curve.short_rates;
	float *this_yield_curve_maturity_time = data->yield_curve.interpolation_times;
	float this_maturity_time = 0;
	float first_maturity_time = 0;
	float last_maturity_time = data->yield_curve.interpolation_times[data->yield_curve.number_of_points - 1];
	float plot_maturity_step = (last_maturity_time - first_maturity_time) / ((float)(data->number_of_plot_points));
	float log_discount_factor = 0;
	for (size_t i = 0; i < data->number_of_plot_points; ++i) {
		this_maturity_time += plot_maturity_step;
		if (this_maturity_time > *this_yield_curve_maturity_time) {
			++this_yield_curve_maturity_time;
			++this_short_rate;
		}
		log_discount_factor -= *this_short_rate * plot_maturity_step;
		data->plot_maturities[i] = this_maturity_time;
		data->plot_short_rates[i] = *this_short_rate;
		data->plot_discount_factors[i] = expf(log_discount_factor);
	}
	ImGui::End();

}


