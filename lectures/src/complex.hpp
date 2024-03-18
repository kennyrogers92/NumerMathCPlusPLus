#ifndef COMPLEX_HPP_INCLUDE
#define COMPLEX_HPP_INCLUDE

#include <iostream>
#include <complex>
#include <cmath>

struct complex_d_t : public std::complex<double>
{
	complex_d_t() : std::complex<double>() {}
	complex_d_t(const double& x) : std::complex<double> (x, 0.) {};
	complex_d_t(const double& x, const double& y) : std::complex<double> (x, y) {};
	complex_d_t(const std::complex<double>& z) : std::complex<double>(z) {};
};


// Output stream of a complex number.
std::ostream& operator<<(std::ostream &output, const complex_d_t &z) {
	if (z.imag() >= 0)
		output << z.real() << " + " << std::abs(z.imag()) << "j";
	else
		output << z.real() << " - " << std::abs(z.imag()) << "j";
	return output;
}

#endif
