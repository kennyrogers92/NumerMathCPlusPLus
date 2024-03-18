#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>
#include "../src/timer.hpp"

double PiArcTanForward(int N) {
	double s = 1.;
	for (int k = 1; k < N; k++) {
		if (k % 2 == 0) {
			s += 1. / (2*k + 1);
		}
		else {
			s -= 1. / (2*k + 1);
		}
	}
	s *= 4;
	return s;
}

double PiArcTanBackward(int N) {
	double s = 0.;
	for (int k = N - 1; k > -1; k--) {
		if (k % 2 == 0) {
			s += 1. / (2*k + 1);
		}
		else {
			s -= 1. / (2*k + 1);
		}
	}
	s *= 4;
	return s;
}

double PiArcTanForwardCompensated(int N) {
	double s = 1., e = 0;
	for (int k = 1; k < N; k++) {
		double x;
		if (k % 2 == 0) {
			x = 1. / (2*k + 1);
		}
		else {
			x = -1. / (2*k + 1);
		}
		double y = x - e;
		double t = s + y;
		e = (t - s) - y;
		s = t;
	}
	return 4*s;
}

int main() {
	timer stopwatch;
	std::cout << std::endl;
	std::cout << "Machine value of pi:\n";
	std::cout << "\t" << std::fixed << std::setprecision(16) << M_PI << std::endl;
	std::cout << std::endl;
	std::cout << "Computed values using the power series of arc tangent function:\n";
	std::cout << std::string(126, '-') << std::endl;
	std::cout << "N\tFORWARD SUM\t\tRELATIVE ERROR\tBACKWARD SUM\t\t"
		<< "RELATIVE ERROR\tCOMPENSATED SUM\t\tRELATIVE ERROR\n";
	std::cout << std::string(126, '-') << std::endl;
	stopwatch.start();
	for (size_t k = 0; k < 10; k++) {
		int N = std::pow(10, k);
		double PiForward = PiArcTanForward(N);
		double PiBackward = PiArcTanBackward(N);
		double PiForwardCompensated = PiArcTanForwardCompensated(N);
		std::cout << std::scientific << std::setprecision(0)
			<< static_cast<double>(N) << "\t"
			<< std::fixed << std::setprecision(16) << PiForward << "\t"
			<< std::scientific << std::setprecision(8)
			<< std::abs(PiForward - M_PI) / M_PI << "\t"
			<< std::fixed << std::setprecision(16) << PiBackward << "\t"
			<< std::scientific << std::setprecision(8)
			<< std::abs(PiBackward - M_PI) / M_PI << "\t"
			<< std::fixed << std::setprecision(16) << PiForwardCompensated << "\t"
			<< std::scientific << std::setprecision(8)
			<< std::abs(PiForwardCompensated - M_PI) / M_PI << std::endl;
	}
	stopwatch.stop();
	std::cout << std::string(126, '-') << std::endl;
	std::cout << std::endl;
	std::cout << "Elapsed time for computing pi using the power series of arc\n"
		<< "tangent function with 1 to 1e+9 terms is " << std::scientific
		<< stopwatch.get_elapsed_time() << "seconds.\n" << std::endl;
	return 0;
}
