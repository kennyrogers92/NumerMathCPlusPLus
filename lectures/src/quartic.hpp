#ifndef QUARTIC_HPP_INCLUDED
#define QUARTIC_HPP_INCLUDED

#include "quadform.hpp"
#include "cardano.hpp"

// Evaluates az^4 + bz^3 + cz^2 + dz + e using nested multiplication.
complex_d_t eval(double a, double b, double c, double d, double e, complex_d_t &z) {
	complex_d_t val;
	val = a*z + b;
	val = val*z + c;
	val = val*z + d;
	val = val*z + e;
	return val;
}

// Implementation of quartic formula.
void quarform(double a, double b, double c, double d, double e, complex_d_t x[]) {
	double A = b/a, B = c/a, C = d/a, D = e/a;
	// Set k = - b/(4.*a)
	double p = B - 3.*A*A/8;
	double q = C - A*B/2. + A*A*A/8;
	double r = D - A*C/4. + A*A*B/16. - 3*A*A*A*A/256.;
	double t = b/(4.*a);
	complex_d_t* roots;
	
	if (D == 0) {
		x[0] = complex_d_t();
		roots = new complex_d_t[3];
		cardano(1, A, B, C, roots);
		x[1] = roots[0];
		x[2] = roots[1];
		x[3] = roots[2];
		delete[] roots;
		return;
	}

	// Biquadratic case
	if (q == 0) {
		roots = new complex_d_t [2];
		quadform(1., p, r, roots);
		for (int i = 0; i < 2; i++) {
			complex_d_t sqrt_root = sqrt(roots[i]);
			x[2*i] = sqrt_root - t;
			x[2*i+1] = -sqrt_root - t;
		}
		delete[] roots;
	}
	else {
		complex_d_t beta [3];
		cardano(1., 2.*p, p*p - 4.*r, -q*q, beta);
		complex_d_t beta0 = beta[1];
		complex_d_t sqrt_beta0 = sqrt(beta0);
		complex_d_t casePos = sqrt(beta0 - 2.*(p + beta0 + q/sqrt_beta0));
		complex_d_t caseNeg = sqrt(beta0 - 2.*(p + beta0 - q/sqrt_beta0));
		x[0] = 0.5*(sqrt_beta0 + casePos) - t;
		x[1] = 0.5*(-sqrt_beta0 + caseNeg) - t;
		x[2] = 0.5*(sqrt_beta0 - casePos) - t;
		x[3] = 0.5*(-sqrt_beta0 - caseNeg) - t;
	}
}

#endif
