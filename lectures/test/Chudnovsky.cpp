#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>
#include "../src/timer.hpp"

double Chudnovsky(int N) {
    double C = 426880.0 * std::sqrt(10005.0);
    double M = 1;
    double L = 13591409;
    double X = 1;
    double sum = M * L / X;
    double pi;
    for (int k = 0; k < N; k++) {
        L = L + 545140134.0;
        X = X * (-262537412640768000);
        M = M * ((12.0 * k + 2) * (12 * k + 6) * (12 * k + 10)) / std::pow(k + 1, 3);
        sum += M * L / X;
    }
    pi = C * 1 / sum;
    return pi;
}

double RelativeError(double exact, double approx) {
	return std::abs(exact - approx) / std::abs(exact);
}

int main() {
	std::cout << std::endl;
	std::cout << "Machine value of pi:\n";
	 	std::cout << "\t" << std::fixed << std::setprecision(48) << M_PI << std::endl;
	std::cout << std::endl;
	std::cout << "Computed values using the Chudnovsky:"
		<< std::endl;
	std::cout << std::string(78, '-') << std::endl;
	std::cout << "N\tCHUDNOVSKY\t\t\t\t\t\tRELATIVE ERROR" << std::endl;
	std::cout << std::string(78, '-') << std::endl;
	for (size_t N = 1; N < 25; N++) {
		double pi = Chudnovsky(N);
		double rel_err = RelativeError(M_PI, pi);
		std::cout << N << "\t" << std::fixed << std::setprecision(48)
          << pi << "\t" << std::scientific << std::setprecision(8)
          << rel_err << std::endl;
	}
	std::cout << std::string(78, '-') << std::endl;
	return 0;
}