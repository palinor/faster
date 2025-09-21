#include <vector>
#include <iostream>

void ReturnNewSorted(std::vector<int>& result, std::vector<int>& input) {
	std::sort(input.begin(), input.end());
	size_t input_size = input.size();
	size_t middle = input_size / 2;
	size_t current_insertion_point = middle;

	for (size_t i = 0; i < input.size(); ++i) {
		result[current_insertion_point] = input[input_size - i];
		if (i % 2) {
			current_insertion_point = middle + i / 2;
			if (current_insertion_point > input_size) current_insertion_point = middle - i / 2;
		}
		else {
			current_insertion_point = middle - i / 2;
		}
	}
}

int main() {
	std::vector<int> input = { 2, 1, 3, 4, 5 };
	std::vector<int> output;
	output.reserve(input.size());

	ReturnNewSorted(output, input);
	std::cout << ouput;
}
