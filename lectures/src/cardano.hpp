#ifndef CARDANO_HPP_INCLUDE
#define CARDANO_HPP_INCLUDE

#include "complex.hpp"

// Evaluates az^3 + bz^2 + cz + d using nested multiplication.
complex_d_t eval(double a, double b, double c, double d, complex_d_t &z) {
	complex_d_t val;
	val = a*z + b;
	val = val*z + c;
	val = val*z + d;
	return val;
}

// Implementation of Cardano's formula for the solutions of the cubic
// equation ax^3 + bx^2 + cx + d = 0.
void cardano(double a, double b, double c, double d, complex_d_t x[]) {
	// parameters in Cardano's formula
	complex_d_t j = complex_d_t(0.0, 1.0);
	complex_d_t q = (3.0*a*c - b*b) / (9.0*a*a);
	complex_d_t r = (9.0*a*b*c - 27.0*a*a*d - 2.0*b*b*b) / (54.0*a*a*a);
	complex_d_t D = q*q*q + r*r;
	complex_d_t u = r + std::sqrt(D);
	complex_d_t v = r - std::sqrt(D);
	complex_d_t s, t;
	if (u.real() >= 0.0)
		s = std::pow(u, 1.0/3.0);
	else
		s = -std::pow(-u, 1.0/3.0);
	if (v.real() >= 0.0)
		t = std::pow(v, 1.0/3.0);
	else
		t = -std::pow(-v, 1.0/3.0);
	// computations of solutions
	x[0] = s + t - b/(3*a);
	x[1] = - 0.5*(s + t) - b/(3*a) + std::sqrt(3.0)*j*(s-t)/2.0;
	x[2] = - 0.5*(s + t) - b/(3*a) - std::sqrt(3.0)*j*(s-t)/2.0;
}

#endif
