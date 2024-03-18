#include <iostream>
#include <iomanip>
#include <string>
#include <cmath>

double f(double x) {
	if (x == 0.0) return 1.0;
	else return (std::exp(x) - 1.0) / x;
}

double F(double x) {
	double z = std::exp(x);
	if (z == 1.0) return 1.0;
	else return (z - 1.0) / std::log(z);
}

int main() {
	std::cout << std::string(70, '-') << std::endl;
	std::cout << " x\t\t\tfl(f(x))\t\t\tfl(F(x))" << std::endl;
	std::cout << std::string(70, '-') << std::endl;
	for (int k = 5; k < 17; k++) {
		double x = std::pow(10.0, -k);
		std::cout << std::scientific << std::setprecision(1) << x << "\t\t"
			<< std::setprecision(16) << f(x) << "\t\t" << F(x) << std::endl;
	}
	std::cout << std::string(70, '-') << std::endl;
	return 0;
}
