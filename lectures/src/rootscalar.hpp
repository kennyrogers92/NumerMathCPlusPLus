// C++ header file for root-finding algorithms.
#ifndef ROOTSCALAR_HPP_INCLUDE
#define ROOTSCALAR_HPP_INCLUDE

// Standard library includes
#include <iostream>
#include <iomanip>
#include <limits>
#include <cmath>
#include <vector>

// Local include
#include "timer.hpp"
#include "max.hpp"
#include "min.hpp"
#include "sgn.hpp"
#include "linalg.hpp"
#include "rootspoly.hpp"

namespace root_scalar {

// Type for a function with input and output having the type double
using UniVarFunction = double(double);

// Struct for parameters in scalar-root finding algorithms
struct param
{
     double    tol   = std::numeric_limits<double>::epsilon();
     int       maxit = 1000;
     double    reftol = 1e-3;
     int       refmax = 100;
     bool      ref = 1;
     double    abstol = 0.9;
     double    funtol = 0.1;

};

// Constant for Inexact Newton-Raphson method
const double eps = std::numeric_limits<double>::epsilon();
const double root_eps = std::sqrt(eps);

// Struct for solution to scalar-root finding problems
struct Result
{
     int            numit;
     int            maxit;
     double         x;
     double         funval;
     double         error;
     double         tol;
     double         elapsed_time;
     std::string    method_name;
     std::string    termination_flag;

     // default constructor
     Result() {}

     // user-defined constructor
     Result (const int& numit, const int& maxit, const double& x,
          const double& funval, const double& error, const double& tol,
          const double& elapsed_time, const std::string& method_name,
          const std::string& termination_flag) {
          this->numit              = numit;
          this->maxit              = maxit;
          this->x                  = x;
          this->funval             = funval;
          this->error              = error;
          this->tol                = tol;
          this->elapsed_time       = elapsed_time;
          this->method_name        = method_name;
          this->termination_flag   = termination_flag;
     }

     void print() {
          std::cout << "ROOT FINDER:                     "
               << method_name << std::endl;
          std::cout << std::setprecision(16);
          std::cout << std::fixed;
          std::cout << "APPROXIMATE ROOT / LAST ITERATE: "
               << x << std::endl;
          std::cout << "TERMINATION:                     "
               << termination_flag << std::endl;
          std::cout << std::scientific;
          std::cout << "FUNCTION VALUE:                  "
               << funval << std::endl;
          std::cout << "ERROR:                           "
               << error << std::endl;
          std::cout << "TOLERANCE:                       "
               << tol << std::endl;
          std::cout << "NUM ITERATIONS:                  "
               << numit << std::endl;
          std::cout << "MAX ITERATIONS:                  "
               << maxit << std::endl;
          std::cout << "ELAPSED TIME:                    "
               << elapsed_time << " seconds" << std::endl;
          std::cout << std::defaultfloat;
     }

     bool operator<(const root_scalar::Result& that) const {
          return numit < that.numit;
     }
};

// Bisection method for approximating a solution of scalar equation f(x) = 0.
Result bisection(UniVarFunction &f, double a, double b, param &parameter) {
     timer stopwatch;
     stopwatch.start();
     std::string term_flag = "Success";
     double err = b - a;
     double fa = f(a), fb = f(b);
     int k = 0;
     double c;

     // Check if endpoints are roots or f has opposite signs on endpoints
     if (fa == 0) {
          c = a; err = 0;
     }
     if (fb == 0) {
          c = b; err = 0;
     }
     if (fa*fb > 0) {
          std::cerr << "Method Fails!" << std::endl;
          return Result();
     }

     // main while loop
     double fc;
     while ((err > parameter.tol) && (k < parameter.maxit)) {
          c = (a + b) / 2.;
          fc = f(c);
          if (fc*fa > 0) {
               a = c, fa = fc;
          }
          else {
               b = c;
          }
          err = abs(b - a);
          k++;
     }
     if ((err > parameter.tol) && (k == parameter.maxit)) {
          term_flag = "Fail";
     }
     stopwatch.stop();
     return Result(k, parameter.maxit, c, fc, err, parameter.tol,
          stopwatch.get_elapsed_time(), "BISECTION", term_flag);
}

// Chord method for approximating a solution of scalar equation f(x) = 0.
Result chord(UniVarFunction &f, double a, double b, double x, param &parameter) {
     timer stopwatch;
     stopwatch.start();
     std::string term_flag = "Success";
     double err = parameter.tol + 1;
     double q = (f(b) - f(a)) / (b - a);
     double fx = f(x);
     int k = 0;

     // main while loop
     while ((err > parameter.tol) && (k < parameter.maxit)) {
          double x_old = x;
          x = x - fx / q;
          fx = f(x);
          err = parameter.abstol * abs(x - x_old)
               + parameter.funtol * abs(fx);
          k++;
     }
     if ((err > parameter.tol) && (k == parameter.maxit)) {
          term_flag = "Fail";
     }
     stopwatch.stop();
     return Result(k, parameter.maxit, x, fx, err, parameter.tol,
          stopwatch.get_elapsed_time(), "Chord", term_flag);
}

// Secant method for approximating a solution of scalar equation f(x) = 0.
Result secant(UniVarFunction &f, double x0, double x1, param &parameter) {
     timer stopwatch;
     stopwatch.start();
     std::string term_flag = "Success";
     double err = parameter.tol + 1;
     double f0 = f(x0), f1 = f(x1);
     int k = 1;

     // main while loop
     while ((err > parameter.tol) && (k < parameter.maxit)) {
          double q = (f1 - f0) / (x1 - x0);
          double x_temp = x1;
          x1 = x1 - f1 / q;
          x0 = x_temp;
          f0 = f1;
          f1 = f(x1);
          err = parameter.abstol * abs(x1 - x0)
               + parameter.funtol * abs(f1);
          k++;
     }
     if ((err > parameter.tol) && (k == parameter.maxit)) {
          term_flag = "Fail";
     }
     stopwatch.stop();
     return Result(k, parameter.maxit, x1, f1, err, parameter.tol,
          stopwatch.get_elapsed_time(), "Secant", term_flag);
}

// Regula Falsi method for approximating a
// solution of scalar equation f(x) = 0.
Result regfalsi(UniVarFunction &f, double x0, double x1, param &parameter) {
     timer stopwatch;
     stopwatch.start();
     std::string term_flag = "Success";
     double err = parameter.tol + 1;
     std::vector<double> x_arr {x0, x1};
     std::vector<double> f_arr {f(x0), f(x1)};
     int k = 1;

     // main while loop
     while ((err > parameter.tol) && (k < parameter.maxit)) {
          double xc = x_arr.at(k), fc = f_arr.at(k);
          int j = k - 1;
          double xj = x_arr.at(j), fj = f_arr.at(j);
          while ((fj * fc >= 0) && (j > 1)) {
               j = j - 1;
               xj = x_arr.at(j);
               fj = f_arr.at(j);
          }
          double q = (fc - fj) / (xc - xj);
          double x = xc - fc / q;
          x_arr.push_back(x), f_arr.push_back(f(x));
          err = parameter.abstol * abs(x - xc)
               + parameter.funtol * abs(f(x));
          k++;
     }
     if ((err > parameter.tol) && (k == parameter.maxit)) {
          term_flag = "Fail";
     }
     stopwatch.stop();
     return Result(k, parameter.maxit, x_arr.back(), f_arr.back(),
          err, parameter.tol, stopwatch.get_elapsed_time(), "RegulaFalsi", term_flag);
}

// Newton-Raphson method for approximating
// a solution of scalar equation f(x) = 0.
Result newtonraphson(UniVarFunction &f, UniVarFunction &df, double x,
     param &parameter) {
     timer stopwatch;
     stopwatch.start();
     std::string term_flag = "Success";
     double err = parameter.tol + 1;
     double fx = f(x);
     int k = 0;

     // main while loop
     while ((err > parameter.tol) && (k < parameter.maxit)) {
          double x_old = x;
          x = x - fx / df(x);
          fx = f(x);
          err = parameter.abstol * abs(x - x_old)
               + parameter.funtol * abs(fx);
          k++;
     }
     if ((err > parameter.tol) && (k == parameter.maxit)) {
          term_flag = "Fail";
     }
     stopwatch.stop();
     return Result(k, parameter.maxit, x, fx, err, parameter.tol,
          stopwatch.get_elapsed_time(), "Newton-Raphson", term_flag);
}

// Inexact Newton-Raphson methods for approximating a
// solution of scalar equation f(x) = 0.
Result inexactforward(UniVarFunction &f, double x, param &parameter) {
     timer stopwatch;
     stopwatch.start();
     std::string term_flag = "Success";
     double err = parameter.tol + 1;
     double fx = f(x);
     int k = 0;

     // main while loop
     while ((err > parameter.tol) && (k < parameter.maxit)) {
          double x_old = x;
          double df = (f(x + root_eps) - f(x)) / root_eps;
          x = x - fx / df;
          fx = f(x);
          err = parameter.abstol * abs(x - x_old)
               + parameter.funtol * abs(fx);
          k++;
     }
     if ((err > parameter.tol) && (k == parameter.maxit)) {
          term_flag = "Fail";
     }
     stopwatch.stop();
     return Result(k, parameter.maxit, x, fx, err, parameter.tol,
          stopwatch.get_elapsed_time(), "Forward", term_flag);
}

Result inexactbackward(UniVarFunction &f, double x, param &parameter) {
     timer stopwatch;
     stopwatch.start();
     std::string term_flag = "Success";
     double err = parameter.tol + 1;
     double fx = f(x);
     int k = 0;

     // main while loop
     while ((err > parameter.tol) && (k < parameter.maxit)) {
          double x_old = x;
          double df = (f(x) - f(x - root_eps)) / root_eps;
          x = x - fx / df;
          fx = f(x);
          err = parameter.abstol * abs(x - x_old)
               + parameter.funtol * abs(fx);
          k++;
     }
     if ((err > parameter.tol) && (k == parameter.maxit)) {
          term_flag = "Fail";
     }
     stopwatch.stop();
     return Result(k, parameter.maxit, x, fx, err, parameter.tol,
          stopwatch.get_elapsed_time(), "Backward", term_flag);
}

Result inexactcenter(UniVarFunction &f, double x, param &parameter) {
     timer stopwatch;
     stopwatch.start();
     std::string term_flag = "Success";
     double err = parameter.tol + 1;
     double fx = f(x);
     int k = 0;

     // main while loop
     while ((err > parameter.tol) && (k < parameter.maxit)) {
          double x_old = x;
          double df = (f(x + root_eps) - f(x - root_eps)) / (2*root_eps);
          x = x - fx / df;
          fx = f(x);
          err = parameter.abstol * abs(x - x_old)
               + parameter.funtol * abs(fx);
          k++;
     }
     if ((err > parameter.tol) && (k == parameter.maxit)) {
          term_flag = "Fail";
     }
     stopwatch.stop();
     return Result(k, parameter.maxit, x, fx, err, parameter.tol,
          stopwatch.get_elapsed_time(), "Center", term_flag);
}

// Steffensen method for approximating a solution of scalar equation f(x) = 0.
Result steffensen(UniVarFunction &f, double x, param &parameter) {
     timer stopwatch;
     stopwatch.start();
     std::string term_flag = "Success";
     double err = parameter.tol + 1;
     double fx = f(x);
     int k = 0;

     // main while loop
     while ((err > parameter.tol) && (k < parameter.maxit)) {
          double x_old = x;
          double q = (f(x + fx) - fx) / fx;
          x = x - fx / q;
          fx = f(x);
          err = parameter.abstol * abs(x - x_old)
               + parameter.funtol * abs(fx);
          k++;
     }
     if ((err > parameter.tol) && (k == parameter.maxit)) {
          term_flag = "Fail";
     }
     stopwatch.stop();
     return Result(k, parameter.maxit, x, fx, err, parameter.tol,
          stopwatch.get_elapsed_time(), "Steffensen", term_flag);
}

// Fixpoint method for approximating a solution of scalar equation f(x) = 0.
Result fixpoint(UniVarFunction &g, double x, param &parameter) {
     timer stopwatch;
     stopwatch.start();
     std::string term_flag = "Success";
     double err = parameter.tol + 1;
     int k = 0;

     // main while loop
     while ((err > parameter.tol) && (k < parameter.maxit)) {
          double x_old = x;
          x = g(x);
          err = abs(x - x_old);
          k++;
     }
     if ((err > parameter.tol) && (k == parameter.maxit)) {
          term_flag = "Fail";
     }
     stopwatch.stop();
     return Result(k, parameter.maxit, x, g(x), err, parameter.tol,
          stopwatch.get_elapsed_time(), "FixPoint", term_flag);
}

// Muller method for approximating a solution of scalar equation f(x) = 0.
Result muller(UniVarFunction &f, double x0, double x1, double x2,
     param &parameter) {
     timer stopwatch;
     stopwatch.start();
     std::string term_flag = "Success";
     double err = parameter.tol + 1;
     double f0 = f(x0), f1 = f(x1), f2 = f(x2);
     int k = 2;

     // main while loop
     while ((err > parameter.tol) && (k < parameter.maxit)) {
          double f01 = (f1 - f0) / (x1 - x0);
          double f12 = (f2 - f1) / (x2 - x1);
          double f012 = (f12 - f01) / (x2 - x0);
          double w = f12 + f012*(x2 - x1);
          double alpha = w*w - 4.*f2*f012;
          double x;
          if (alpha >= 0.) {
               // Reference used max(w-sqrt(alpha), w+sqrt(alpha))
               double d = max(w - std::sqrt(alpha), w + std::sqrt(alpha));
               x = x2 - 2.*f2/d;
          }
          else {
               x = x2 - f2/f12;
          }
          x0 = x1, x1 = x2, x2 = x;
          f0 = f1; f1 = f2; f2 = f(x);
          err = parameter.abstol * abs(x2 - x1)
               + parameter.funtol * abs(f2);
          k++;
     }
     if ((err > parameter.tol) && (k == parameter.maxit)) {
          term_flag = "Fail";
     }
     stopwatch.stop();
     return Result(k, parameter.maxit, x2, f2, err, parameter.tol,
          stopwatch.get_elapsed_time(), "Muller", term_flag);
}

// Generalization of secant and muller method for approximating
// a solution of scalar equation f(x) = 0.
Result rootpolyinterp(UniVarFunction &f, std::vector<double>& x,
     param &parameter) {
     timer stopwatch;
     stopwatch.start();
     std::string term_flag = "Success";
     size_t n = x.size();
     double xc= x[n-1];
     double fx;
     double err = parameter.tol + 1.;
     int k = n - 1;
     // main while loop
     while ((err > parameter.tol) && (k < parameter.maxit)) {
          // Calculate the polynomial p passing thru x
          std::vector<double> p = linalg::polyinterp(f, x);
          std::vector<complex_d_t> q;
          for (size_t i = 0; i < n; i++) {
               q.push_back(p[i]);
          }
          // Compute the real roots of p
          std::vector<double> z;
          roots_poly::param roots_poly_parameter;
          std::vector<complex_d_t> temp =
               roots_poly::newtonHorner(q, xc, roots_poly_parameter).z_arr;
          // Filter the real roots of p
          for (size_t i = 0; i < n; i++) {
               if (abs(temp[i].imag()) <= eps) {
                    z.push_back(temp[i].real());
               }
          }
          // Choose the minimum of the real roots
          xc = min(z);
          fx = f(xc);
          err = parameter.abstol * abs(xc - x[n-1])
               + parameter.funtol * abs(fx);
          for (size_t j = 0; j < n - 1; j++) {
               x[j] = x[j+1];
          }
          x[n - 1] = xc;
          k++;
     }
     if ((err > parameter.tol) && (k == parameter.maxit)) {
          term_flag = "Fail";
     }
     stopwatch.stop();
     return Result(k, parameter.maxit, xc, fx, err, parameter.tol,
          stopwatch.get_elapsed_time(), "RootPolyInterp", term_flag);
}

// Generalization of newton method for approximating
// a solution of scalar equation f(x) = 0.
Result sidisecant(UniVarFunction &f, std::vector<double> x, param &parameter) {
     timer stopwatch;
     stopwatch.start();
     std::string term_flag = "Success";
     size_t n = x.size();
     complex_d_t xc= x[n-1];
     double fx;
     double err = parameter.tol + 1.;
     int k = n - 1;
     while ((err > parameter.tol) && (k < parameter.maxit)) {
          // Calculate the polynomial p passing thru x
          std::vector<double> p = linalg::polyinterp(f, x);
          std::vector<complex_d_t> q;
          for (size_t i = 0; i < n; i++) {
               q.push_back(p[i]);
          }
          // To hold function value of p under x[n-1] and the remainder q 
          complex_d_t qx;
          std::vector<complex_d_t> dq;
          std::tie(qx, dq) = roots_poly::horner(q, xc);
          complex_d_t dqx = roots_poly::horner(dq, xc).first;
          xc = xc - qx/dqx;
          fx = f(xc.real());
          err = parameter.abstol * abs(xc - x[n - 1])
               + parameter.funtol * abs(fx);
          for (size_t j = 0; j < n - 1; j++) {
               x[j] = x[j + 1];
          }
          x[n - 1] = xc.real();
          k++;
     }
     if ((err > parameter.tol) && (k == parameter.maxit)) {
          term_flag = "Fail";
     }
     stopwatch.stop();
     return Result(k, parameter.maxit, x[n-1], fx, err, parameter.tol,
          stopwatch.get_elapsed_time(), "SidiSecant", term_flag);
}

// Modified Dekker-Brent method that uses quadratic interpolation
// for approximating a solution of scalar equation f(x) = 0.
Result dekkerbrent(UniVarFunction &f, double a, double b, param &parameter) {
     timer stopwatch;
     stopwatch.start();
     std::string term_flag = "Success";
     double err = parameter.tol + 1;
     double delta = parameter.tol + 2.*eps*abs(b);
     double fa = f(a), fb = f(b);
     int k = 0;

     // Check endpoints if they are zeros
     if (fa == 0) {
          b = a;
          err = 0;
     }
     if (fb == 0) {
          err = 0;
     }
     // Interchange a and b if f(a) is closer to zero
     if (abs(fa) < abs(fb)) {
          double temp = a, f_temp = fa;
          a = b, fa = fb;
          b = temp, fb = f_temp;
     }
     // Initialize third point c
     double c = a, fc = fa;

     // main while loop
     while ((err > delta) && (k < parameter.maxit)) {
          double b_old = b;
          double z;
          
          // linear interpolation
          if ((a == b) || (b == c) || (a == c)) {
               z = b - fb*(b - a)/(fb - fa);
          }
          // quadratic interpolation
          else {
               double f01 = (fa - fc)/(a - c);
               double f12 = (fb - fa)/(b - a);
               double f012 = (f12 - f01)/(b - c);
               double w = f12 + f012*(b - a);
               double alpha = w*w - 4.*fb*f012;
               if (alpha >= 0) {
                    double d = w + sgn(w)*std::sqrt(alpha);
                    z = b - 2.*fb/d;
               }
               else {
                    z = b - fb / f12;
               }
          }

          double m = (c - b)/2.;   // Half the distance between b and c
          double b_temp;           // Hold the candidate for b
          // If z lies in between b and b+m, take z as candidate
          // Otherwise, take the midpoint of b and c as candidate
          if ( (min(b, b + m) < z) && (z < max(b, b + m)) ) {
               b_temp = z;
          }
          else {
               b_temp = b + m;
          }

          // If candidate is far enough from previous b, update b
          // Otherwise, perturbed slightly from b
          if (abs(b - b_temp) > delta) {
               b = b_temp;
          }
          else {
               b = b + delta*sgn(m);
          }
          
          // Update a, b, c
          a = b_old, fa = fb;
          fb = f(b);
          if (fa*fb < 0) {
               c = a, fc = fa;
          }
          if (abs(fc) < abs(fb)) {
               a = b, fa = fb;
               b = c, fb = fc;
               c = a, fc = fa;
          }
          delta = parameter.tol + 2.*eps*abs(b);
          err = abs(m);
          k++;
          if (abs(fb) <= parameter.tol) {
               break;
          }
     }
     if ((err > delta) && (k == parameter.maxit)) {
          term_flag = "Fail";
     }
     stopwatch.stop();
     return Result(k, parameter.maxit, b, fb, err, parameter.tol,
          stopwatch.get_elapsed_time(), "Modified Dekker-Brent", term_flag);
}

// Aitken method for accelerating fix point method approximating a solution of
// the scalar equation f(x) = 0.
Result aitken(UniVarFunction &g, double x, param &parameter) {
     timer stopwatch;
     stopwatch.start();
     std::string term_flag = "Success";
     double err = parameter.tol + 1;
     int k = 0;

     // main while loop
     while ((err > parameter.tol) && (k < parameter.maxit)) {
          double x_old = x;
          double x1 = g(x);
          double x2 = g(x1);
          x = x2 - std::pow((x2 - x1), 2)/(x2 - 2.*x1 + x);
          err = abs(x - x_old);
          k++;
     }
     if ((err > parameter.tol) && (k == parameter.maxit)) {
          term_flag = "Fail";
     }
     stopwatch.stop();
     return Result(k, parameter.maxit, x, g(x), err, parameter.tol,
          stopwatch.get_elapsed_time(), "Aitken", term_flag);
}

// Accelerating Newton-Raphson method for
// approximating a solution of scalar equation f(x) = 0. 
Result modnewton(UniVarFunction &f, UniVarFunction &df, double x,
     long m, param &parameter) {
     timer stopwatch;
     stopwatch.start();
     std::string term_flag = "Success";
     double err = parameter.tol + 1;
     double fx = f(x);
     int k = 0;

     // main while loop
     while ((err > parameter.tol) && (k < parameter.maxit)) {
          double x_old = x;
          x = x - m*fx/df(x);
          fx = f(x);
          err = parameter.abstol * abs(x - x_old) + parameter.funtol * abs(fx);
          k++;
     }
     if ((err > parameter.tol) && (k == parameter.maxit)) {
          term_flag = "Fail";
     }
     stopwatch.stop();
     return Result(k, parameter.maxit, x, fx, err, parameter.tol,
          stopwatch.get_elapsed_time(), "Modified Newton", term_flag);
}

// Adaptive newton-raphson method for accelerating fix point method approximating
// a solution of the scalar equation f(x) = 0. */
 Result adaptivenewton(UniVarFunction &f, UniVarFunction &df, double x,
     param &parameter, double m0, double alpha, double _lambda) {
     timer stopwatch;
     stopwatch.start();
     std::string term_flag = "Success";
     double
          err = parameter.tol + 1,
          err_abs = err,
          fx = f(x),
          m = m0,
          lambda = 1;
     int k = 0;

     // main while loop
     while ((err > parameter.tol) && (k < parameter.maxit)) {
          double x_old = x;
          x = x - m*fx/df(x);
          fx = f(x);
          double err_temp = err_abs;
          err_abs = abs(x - x_old);
          err = parameter.abstol*err_abs + parameter.funtol*abs(fx);
          double lambda_old = lambda;
          lambda = err_abs / err_temp;
          if ((abs(lambda - lambda_old) < alpha) && (lambda > _lambda)) {
               double m_temp = 1.0 / abs(1. - lambda);
               if (m < m_temp) {
                    m = m_temp;
               }
          }
          k++;
     }
     if ((err > parameter.tol) && (k == parameter.maxit)) {
          term_flag = "Fail";
     }
     stopwatch.stop();
     return Result(k, parameter.maxit, x, fx, err, parameter.tol,
          stopwatch.get_elapsed_time(), "Adaptive Newton", term_flag);
}

}    // end of namespace rootscalar

#endif