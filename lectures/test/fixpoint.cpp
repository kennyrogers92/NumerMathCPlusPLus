#include "../src/rootscalar.hpp"

// Define g(x) = 0.5cos(2x)
double g(double x) {
	return 0.5*std::cos(2.0*x);
}

int main() {
	// print filename
	std::cout << "File: " << __FILE__ << std::endl;
	// parameter object
	root_scalar::param parameter;
	parameter.tol = 1e-15;
	parameter.maxit = 100;
	// approximate a root by bisection method
	root_scalar::Result result = fixpoint(g, 0.5, parameter);
	result.print();
	return 0;
}
