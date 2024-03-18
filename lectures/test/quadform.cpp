//Implemetation of quadratic formula
#include <iomanip>
#include "../src/quadform.hpp"

int main() {
	// coefficients
	double a = 2.0, b = 3.0, c = 4.0;
	// pre-allocate array of solutions
	complex_d_t x[2];
	//apply quadratic formula
	quadform(a, b, c, x);
	std::cout << "The solutions of the quadratic equation f(x) = "
		<< a << "x^2 + " << b << "x + " << c << " = 0 are: " << std::endl;
	// print solutions and function values
	std::cout << std::scientific << std::setprecision(8);
	for (int k = 0; k < 2; k++) {
		std::cout << "x = " << x[k] << "\t";
		std::cout << "f(x) = " << eval(a, b, c, x[k]) << std::endl;
	}
	return 0;
}
