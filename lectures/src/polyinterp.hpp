/**
 * @file polyinterp.cpp
 * @author kennyrogers (kenraposas92@gmail.com)
 * @brief C++ Header file for polynomial interpolation
 * @version 0.1
 * @date 2024-03-14
 * 
 * @copyright Copyright (c) 2024
 */
#ifndef POLYINTERP
#define POLYINTERP

// Standard Library Includes
#include <vector>
#include <iostream>
#include <iomanip>
#include <cmath>

// Local Includes
#include "max.hpp"
#include "rootspoly.hpp"

namespace polyinterp {

    // Vector_t are representation of polynomials.
    // Dont use it in any other context as addition on vectors of different
    // dimensions are permissible when viewed as polynomials.
    using vector_t = std::vector<double>;
    using matrix_t = std::vector<vector_t>;
    using UniVarFunction = double(double);

    // Output stream for vector_t and matrix_t
    std::ostream& operator<<(std::ostream& os, const vector_t& X) {
        int n = X.size();
        for (int k = n - 1; k >= 0; k--) {
            double coeff = X[k];
            if (coeff == 0.) {
                continue;
            }
            if (coeff < -1.) {
                coeff *= -1.;
                if (k == n - 1) {
                    os << "-";
                }
                else {
                    os << " - ";
                }
            }
            else {
                if (k < n - 1) {
                    os << " + ";
                }
            }
            if (coeff != 1. || k == 0) {
                os << coeff;
            }
            if (k == 1) {
                os << "x";
            }
            if (k > 1) {
                os << "x^" << k;
            }
        }
        return os;
    }

    // Arithmetic Operations on Vectors
    vector_t operator+(const vector_t& A, const vector_t& B) {
        size_t n = A.size();
        size_t m = B.size();
        size_t deg = max(n, m);

        vector_t C;
        C.reserve(deg);

        for (size_t i = 0; i < deg; i++) {
            double a = (i < n) ? A[i] : 0.0;
            double b = (i < m) ? B[i] : 0.0;
            C.push_back(a + b);
        }
        return C;
    }

    vector_t operator+=(vector_t& A, const vector_t& B) {
        size_t m = B.size();
        if (A.size() < m) {
            A.resize(m, 0.);
        }
        A = A + B;
        return A;
    }

    vector_t operator-(const vector_t& A) {
        size_t n = A.size();
        vector_t B (n, 0);
        for (size_t i = 0; i < n; i++) {
            B[i] = -A[i];
        }
        return B;
    }

    vector_t operator-(const vector_t& A, const vector_t& B) {
        return A + (-B);   
    }

    // Scalar multiplication with a vector
    vector_t operator*(const vector_t& A, const double& x) {
        size_t n = A.size();
        vector_t C (n, 0);
        for (size_t i = 0; i < n; i++) {
            C[i] = x*A[i];
        }
        return C;
    }

    vector_t operator*(const double& x, const vector_t& A) {
        return A*x;
    }

    vector_t operator/(const vector_t& A, const double x) {
        return A*(1./x);
    }

    // Definition of multiplication of polynomials stored as a vector
    vector_t operator*(const vector_t& p, const vector_t& q) {
        vector_t ans;
        size_t n = p.size() + q.size() - 2;
        for (size_t k = 0; k <= n; k++) {
            double c = 0.;
            for (size_t j = 0; j <= k; j++) {
                if ((j < p.size()) && (k - j < q.size())) {
                    c += p[j]*q[k - j];
                }
            }
            ans.push_back(c);
        }
        return ans;
    }

    // Lagrange Interpolating polynomial of f given x nodes
    vector_t lagrange(const UniVarFunction& f, const vector_t& x) {
        size_t n = x.size();
        vector_t p {0.};
        for (size_t k = 0; k < n; k++) {
            vector_t l {1.};
            for (size_t i = 0; i < n; i++) {
                if (i == k) {
                    continue;
                }
                vector_t temp {-x[i], 1.};
                temp = temp/(x[k] - x[i]);
                l = l*temp;
            }
            p += f(x[k])*l;
        }
        return p;
    }

    // Newton Form of the Lagrange Interpolating Polynomial
    vector_t newtonlagrange(const UniVarFunction& f, const vector_t& x) {
        size_t n = x.size() - 1;
        matrix_t d (n + 1, vector_t (n + 1, 0));
        for (size_t k = 0; k <= n; k++) {
            d[k][0] = f(x[k]);
        }
        for (size_t k = 1; k <= n; k++) {
            for (size_t l = k; l <= n; l++) {
                d[l][k] = (d[l][k - 1] - d[l - 1][k - 1])/(x[l] - x[l - k]);
            }
        }
        vector_t p {d[n][n]};
        for (size_t k = 1; k <= n; k++) {
            vector_t temp {-x[n - k], 1.};
            p = p*temp + vector_t {d[n - k][n - k]};
        }
        return p;
    }

    // Hermite Interpolating Polynomial
    vector_t hermite(const UniVarFunction& f, const UniVarFunction& df,
        const vector_t& x) {
        size_t n = x.size();
        vector_t p {0.};
        for (size_t k = 0; k < n; k++) {
            double delta = 1.;
            vector_t nu {1.};
            double ell1 = 0.;
            for (size_t i = 0; i < n; i++) {
                if (i == k) {
                    continue;
                }
                nu = nu*vector_t{-x[i], 1.};
                delta = delta * (x[k] - x[i]);
                ell1 += 1./(x[k] - x[i]);
            }
            vector_t ell2 = (nu/delta)*(nu/delta);
            vector_t eta = vector_t{-x[k], 1.}*ell2;
            vector_t h = (vector_t{1.} - 2.*ell1*vector_t{-x[k], 1.})*ell2;
            p += f(x[k])*h + df(x[k])*eta;
        }
        return p;
    }

    // Lagrange Interpolating Polynomial at affine-transformed Chebyshev Points
    vector_t chebyshev(const UniVarFunction& f, size_t n,
        double a=0., double b=1.) {
        // Construct n+1 Chebyshev points
        vector_t x;
        for (size_t k = 0; k <=n; k++) {
            double temp = cos((k+0.5)*M_PI/(n+1.));
            temp = 0.5*((b - a)*temp + a + b);
            x.push_back(temp);
        }
        return newtonlagrange(f, x);
    }

    // Construct the nth Legendre Polynomial
    vector_t legendre(int n) {
        if (n == 0) {
            return vector_t{1.};
        }
        if (n == 1) {
            return vector_t{0., 1.};
        }
        return vector_t{0., (2. - 1./n)} * legendre(n - 1)
            + vector_t{-1. + 1./n} * legendre(n - 2);
    }


    // Lagrange Interpolating Polynomial at affine-transformed Legendge Points
    vector_t legendre(const UniVarFunction&f, size_t n,
        double a=0., double b=1.) {
        // Construct the legendre polynomial
        vector_t L = legendre(n+1);
        std::vector<complex_d_t> p;
        for (size_t k = 0; k < L.size(); k++) {
            p.push_back(L[k]);
        }
        // get the zeros of legendre polynomial
        roots_poly::param parameter;
        std::vector<complex_d_t> z;
        z = roots_poly::newtonHorner(p, a+b/2., parameter).z_arr;
        vector_t x;
        for (size_t k = 0; k < z.size(); k++) {
            x.push_back(0.5*(b - a)*z[k].real() + a + b);
        }
        return newtonlagrange(f, x);
    }

}

#endif