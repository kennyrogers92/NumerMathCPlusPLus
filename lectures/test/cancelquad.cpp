#include <iostream>
#include <iomanip>
#include <cmath>

int main() {
	// Roots of x^2 - 26x + 1 = 0 using single precision.
	float x0_f = 13.0f + std::sqrt(168.0f);
	float x1_f = 13.0f - std::sqrt(168.0f);
	// Roots of x^2 - 26x + 1 = 0 using double precision.
	double x0_d = 13.0 + std::sqrt(168.0);
	double x1_d = 13.0 - std::sqrt(168.0);
	// Display relative errors.
	std::cout << "Relative error of first root\t\t\t";
	std::cout << std::abs(x0_d - static_cast<double>(x0_f)) / std::abs(x0_f)
		<< std::endl;
	std::cout << "Relative error of second root\t\t\t";
	std::cout << std::abs(x1_d - static_cast<double>(x1_f)) / std::abs(x1_f)
		<< "\t\t" << std::endl;
	// Alternative way of computing the second root.
	x1_f = 1.0f / x0_f;
	x1_d = 1.0 / x0_d;
	std::cout << "Relative error of (alternative) second root\t";
	std::cout << std::abs(x1_d - static_cast<double>(x1_f)) / std::abs(x1_d)
		<< "\t" << std::endl;
	return 0;
}	
