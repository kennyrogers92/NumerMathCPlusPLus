// Sample implementation of quartic formula.
#include <iomanip>
#include <vector>
#include "../src/quartic.hpp"

int main() {
	// coefficients
	double coeff1 [] = {2., 3., 4., 5., 6.};
	double coeff2 [] = {1., -2., 3., -4., 0.};
	double coeff3 [] = {4., 0., -9., 0., 2.};
	double* coeffmat [] = {coeff1, coeff2, coeff3};
	// pre-allocate array of solutions
	complex_d_t x_1[4], x_2[4], x_3[4];
	complex_d_t* roots [] = {x_1, x_2, x_3};
	// apply quartic formula
	for (int i = 0; i < 3; i++) {
		quarform(coeffmat[i][0], coeffmat[i][1], coeffmat[i][2],
			coeffmat[i][3], coeffmat[i][4], roots[i]);
	}
	// Print out solutions
	for (int i = 0; i < 3; i++) {
		std::cout << "The solutions of the quartic equation f(x) = "
			<< coeffmat[i][0] << "x^4 + " << coeffmat[i][1] << "x^3 + "
			<< coeffmat[i][2] << "x^2 + " << coeffmat[i][3] << "x + "
			<< coeffmat[i][4] << " = 0 are: " << std::endl;
		for (int j = 0; j < 4; j++) {
			std::cout << std::scientific << std::setprecision(16)
				<< "x = " << roots[i][j] << "\t";
			std::cout << "f(x) = " << eval(coeffmat[i][0], coeffmat[i][1],
				coeffmat[i][2], coeffmat[i][3], coeffmat[i][4], roots[i][j]) << std::endl;
		}
		std::cout << std::defaultfloat << std::endl;
	}
	
 	return 0;
}
