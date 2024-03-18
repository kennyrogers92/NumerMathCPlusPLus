#ifndef QUADFORM_HPP_INCLUDE
#define QUADFORM_HPP_INCLUDE

#include "complex.hpp"

//Evaluates az^2 + bz + c.
complex_d_t eval(complex_d_t a, complex_d_t b, complex_d_t c, complex_d_t &z) {
	return (a*z + b)*z + c;
}

//Implementation of the quadratic formula for the solutions of the equation ax^2 + bx + c = 0.
void quadform(double a, double b, double c, complex_d_t x[]) {
	complex_d_t sqrt_disc = std::sqrt(complex_d_t(b*b - 4.0*a*c));
	x[0] = (-b + sqrt_disc) / (2.0*a);
	x[1] = (-b - sqrt_disc) / (2.0*a);
	return;
}

#endif
