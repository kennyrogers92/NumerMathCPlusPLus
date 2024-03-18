// Sample implementation of Cardano's formula. 
#include <iomanip>
#include "../src/cardano.hpp"

int main() {
	// coefficients
	double a = 2.0, b = 3.0, c = 4.0, d = 5.0;
	// pre-allocate array of solutions
	complex_d_t x[3];
	// apply Cardano's formula
	cardano(a, b, c, d, x);
	std::cout << "The solutions of the cubic equation f(x) = "
		<< a << "x^3 + " << b << "x^2 + " << c << "x + " << d << " = 0 "
	       	<< "are:" << std::endl;
	// print solutions and function values
	std::cout << std::scientific << std::setprecision(8);
	for (int k = 0; k < 3; k++) {
		std::cout << "x = " << x[k] << "\t";
		std::cout << "f(x) = " << eval(a, b, c, d, x[k]) << std::endl;
	}
	return 0;
}
