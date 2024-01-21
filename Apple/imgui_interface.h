//
//  imgui_interface.h
//  example_apple_metal
//
//  Created by Aion Feehan on 1/18/24.
//  Copyright © 2024 Warren Moore. All rights reserved.
//

#ifndef imgui_interface_h
#include <stdio.h>
#include "imgui.h"
#include "../implot/implot.h"
#include "../yield_curve.h"

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
    float slider_width);
void displayDataFree(DisplayData *data);
void imGuiRenderLoop(DisplayData *data);


#define imgui_interface_h


#endif /* imgui_interface_h */
