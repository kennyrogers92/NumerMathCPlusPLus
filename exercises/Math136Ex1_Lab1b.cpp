/* Math 136 Exercise 1
 * Name:  Krizelda Claire Cayabyab
 *        Kenneth Raposas
 *        Hannah Denise Sarilla
 * Date: 25 March 2022
 */

#include <iostream>
#include <iomanip>
#include <string>
#include <cmath>

double PiArcSinCompensated(int N) {
	double s = 1, x = 0.125, e = 0;
	for (int i = 1; i < N; i++) {
		double y = x / (2*i + 1) - e;
		double t = s + y;
		e = (t - s) - y;
		s = t;
		x *= 0.125*(2*i + 1) / (i + 1);
	}
	return 3*s;
}

double RelativeError(double exact, double approx) {
	return std::abs(exact - approx) / std::abs(exact);
}

int main() {
	std::cout << std::endl;
	std::cout << "Machine value of pi:\n";
	 	std::cout << "\t" << std::fixed << std::setprecision(16)
		<< M_PI << std::endl;
	std::cout << std::endl;
	std::cout << "Computed values using the power series of"
		<< " arc tangent function:" << std::endl;
	std::cout << std::string(46, '-') << std::endl;
	std::cout << "N\tCOMPENSATED SUM\t\tRELATIVE ERROR" << std::endl;
	std::cout << std::string(46, '-') << std::endl;
	for (size_t N = 1; N < 25; N++) {
		double pi = PiArcSinCompensated(N);
		double rel_err =RelativeError(M_PI, pi);
		std::cout << N << "\t" << std::fixed << std::setprecision(16)
          << pi << "\t" << std::scientific << std::setprecision(8)
          << rel_err << std::endl;
	}
	std::cout << std::string(46, '-') << std::endl;
	return 0;
}
