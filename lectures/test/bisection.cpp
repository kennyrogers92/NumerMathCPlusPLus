// Sample implementation of bisection method.
#include "../src/rootscalar.hpp"

// Define f(x) = 0.25*cos^2(2x) - x^2
double f(double x) {
	return 0.25*std::pow(std::cos(2.0*x), 2) - std::pow(x, 2);
}

int main() {
	// print filename
	std::cout << "File: " << __FILE__ << std::endl;
	// parameter object
	root_scalar::param parameter;
	parameter.tol = 1e-15;
	parameter.maxit = 100;
	// approximate a root by bisection method
	root_scalar::Result result = bisection(f, 0.0, 1.0, parameter);
	result.print();
	return 0;
}
