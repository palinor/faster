// start with pricing a call option using bs dynamics / barrier option / american exercise call option
// once the scheme is implemented and some tests are run, try expanding the model

#include "arena_allocator.h"
#include "matrix.cpp"

struct PDE1DGridf32 {
	float *x_grid;
	size_t x_grid_size;
	float *t_grid;
	size_t t_grid_size;
};

enum class SchemeType {
	EXPLICIT,
	IMPLICIT,
	CRANK_NICOLSON
};

void FillTransitionMatrix(
	MatrixTridiagonalf32 *transition_matrix,
	PDE1DGridf32 *grid,
	float(*drift_function)(float, float),
	float(*variance_function)(float, float),
	float(*interest_rate_function)(float, float),
	float time_to_eval_matrix,
	float theta_delta_t_coef
) {

	// note that the matrix is implemented as (I + theta_delta_t_coef * A)
	// 
	// todo(AION) the boundary conditions go here
	// first row
	float x = grid->x_grid[0];
	float delta_x = grid->x_grid[1] - x;
	float drift = drift_function(time_to_eval_matrix, x);
	float variance = variance_function(time_to_eval_matrix, x);
	float interest_rate = interest_rate_function(time_to_eval_matrix, x);
	transition_matrix->upper_diagonal[0] = (drift / (2 * delta_x) + variance / (2 * delta_x * delta_x)) * theta_delta_t_coef;
	transition_matrix->diagonal[0] = (-interest_rate - variance / (delta_x * delta_x)) * theta_delta_t_coef + 1;

	// last row 
	size_t n = grid->x_grid_size;
	x = grid->x_grid[n - 1];
	delta_x = x - grid->x_grid[n - 2];
	drift = drift_function(time_to_eval_matrix, x);
	variance = variance_function(time_to_eval_matrix, x);
	interest_rate = interest_rate_function(time_to_eval_matrix, x);
	transition_matrix->diagonal[n - 1] = (-interest_rate - variance / (delta_x * delta_x) + 1) * theta_delta_t_coef / 2;
	transition_matrix->lower_diagonal[n - 1] = (-drift / (2 * delta_x) + variance / (2 * delta_x * delta_x)) * theta_delta_t_coef / 2 + 1;

	// fill the rest
	for (size_t x_index = 1; x_index < grid->x_grid_size - 1; ++x_index) {
		x = grid->x_grid[x_index];
		delta_x = grid->x_grid[x_index + 1] - x;
		drift = drift_function(time_to_eval_matrix, x);
		variance = variance_function(time_to_eval_matrix, x);
		interest_rate = interest_rate_function(time_to_eval_matrix, x);
		transition_matrix->upper_diagonal[x_index] = (drift / (2 * delta_x) + variance / (2 * delta_x * delta_x)) * theta_delta_t_coef / 2;
		transition_matrix->diagonal[x_index] = (-interest_rate - variance / (delta_x * delta_x)) * theta_delta_t_coef + 1;
		transition_matrix->lower_diagonal[x_index] = (-drift / (2 * delta_x) + variance / (2 * delta_x * delta_x)) * theta_delta_t_coef;
	}
}

//todo(AION) we need to make sure to include spatial boundary conditions as well
int FiniteDifferenceSolverf32(
	Vectorf32 *current_x_state,
	PDE1DGridf32 *grid,
	Vectorf32 *terminal_condition,
	float(*drift_function)(float, float),
	float(*variance_function)(float, float),
	float(*interest_rate_function)(float, float),
	Arena *arena,
	SchemeType scheme_type
) {
	MatrixTridiagonalf32 transition_matrix;
	Vectorf32 right_hand_side_vector;
	Vectorf32 scratch_vector;
	transition_matrix.dimension = grid->x_grid_size;
	transition_matrix.upper_diagonal = (float *)ArenaGetMemory(grid->x_grid_size * sizeof(float), arena);
	transition_matrix.lower_diagonal = (float *)ArenaGetMemory(grid->x_grid_size * sizeof(float), arena);
	transition_matrix.diagonal = (float *)ArenaGetMemory(grid->x_grid_size * sizeof(float), arena);
	scratch_vector.size = grid->x_grid_size;
	right_hand_side_vector.size = grid->x_grid_size;
	scratch_vector.contents = (float *)ArenaGetMemory(grid->x_grid_size * sizeof(float), arena);
	right_hand_side_vector.contents = (float *)ArenaGetMemory(grid->x_grid_size * sizeof(float), arena);

	memcpy(current_x_state->contents, terminal_condition->contents, terminal_condition->size);

	for (int time_step = grid->t_grid_size - 1; time_step > 0; --time_step) {

		// at this point, time_step = t + 1 in the equations (it is the time for the right hand side of the equation)
		float delta_t = grid->t_grid[time_step] - grid->t_grid[time_step - 1];
		if (scheme_type == SchemeType::CRANK_NICOLSON) {
			//CRANK-NICOLSON scheme: theta = 1/2, supposedly more stable. One forward pass, one matrix inversion.
			float time_to_eval_matrix = (grid->t_grid[time_step] + grid->t_grid[time_step - 1]) / 2;
			float theta_delta_t_coef = delta_t / 2;

			// fill the transition matrix for right-hand side step
			FillTransitionMatrix(
				&transition_matrix,
				grid,
				drift_function,
				variance_function,
				interest_rate_function,
				time_to_eval_matrix,
				theta_delta_t_coef
			);

			// compute right hand side of the expression
			MatrixTridiagonalf23ApplyToVector(&right_hand_side_vector, &transition_matrix, current_x_state);

			// set up transition matrix for backwards step
			//used to be (I + delta_t/2 A), convert this to (I - delta_t/2 A)
			MatrixTridiagonalf32MultiplyConstant(&transition_matrix, -1);
			for (size_t x_index = 0; x_index < transition_matrix.dimension; ++x_index) {
				transition_matrix.diagonal[x_index] += 2;
			}
			// compute the backwards step
			MatrixTridiagonalf32InvertEquation(current_x_state, &transition_matrix, &right_hand_side_vector, &scratch_vector);
		}
		else if (scheme_type == SchemeType::EXPLICIT) {
			float time_to_eval_matrix = grid->t_grid[time_step];
			float theta_delta_t_coef = delta_t;

			// fill the transition matrix for right-hand side step
			// EXPLICIT scheme: theta = 0, we evaluate fully at t + 1 to step back to t. No inversion required.
			FillTransitionMatrix(
				&transition_matrix,
				grid,
				drift_function,
				variance_function,
				interest_rate_function,
				time_to_eval_matrix,
				theta_delta_t_coef
			);
			// compute right hand side of the expression
			memcpy(right_hand_side_vector.contents, current_x_state->contents, current_x_state->size);
			MatrixTridiagonalf23ApplyToVector(current_x_state, &transition_matrix, &right_hand_side_vector);
			// and we are done
		}
		else if (scheme_type == SchemeType::IMPLICIT) {
			float time_to_eval_matrix = grid->t_grid[time_step - 1];
			float theta_delta_t_coef = -delta_t;

			// fill the transition matrix for right-hand side step
			// IMPLICIT scheme: theta = 1, we evaluate fully at t to step back to t. No right hand side calculation required, only one inversion.
			FillTransitionMatrix(
				&transition_matrix,
				grid,
				drift_function,
				variance_function,
				interest_rate_function,
				time_to_eval_matrix,
				theta_delta_t_coef
			);
			memcpy(right_hand_side_vector.contents, current_x_state->contents, current_x_state->size);
			MatrixTridiagonalf32InvertEquation(current_x_state, &transition_matrix, &right_hand_side_vector, &scratch_vector);
			// and we are done
		}
	}
	return 0;

}

