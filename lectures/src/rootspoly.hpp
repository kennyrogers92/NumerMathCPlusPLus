/**
 * @file rootspoly.cpp
 * @author kennyrogers (kenraposas92@gmail.com)
 * @brief C++ Header file for finding roots of real polynomial of degree
 * greater than 5. For real polynomials of degree less than 5, use quadform,
 * cardano and quarform for degree 2, 3, and 4, respectively. Real polynomials
 * of the form a0 + a1x + ... + anx^n where are stored in a vector in the form
 * [a0, a1, ..., an].
 * @version 0.1
 * @date 2024-03-10
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef ROOTSPOLY_HPP
#define ROOTSPOLY_HPP

// Standard Library Includes
#include <limits>
#include <vector>
#include <utility>
#include <iostream>
#include <iomanip>

// Local includes
#include "timer.hpp"
#include "complex.hpp"
#include "sgn.hpp"
#include "max.hpp"

namespace roots_poly {

// Constant for machine epsilon
const double eps = std::numeric_limits<double>::epsilon();

// Struct for parameters of rootspoly
struct param 
{
    int     maxit = 1000;
    int     refmax = 100;
    double  tol = eps;
    double  reftol = 1e-3;
    bool    ref = true;
    double  abstol = 0.9;
    double  funtol = 0.1;
};

// Struct for roots of real polynomials
struct Result
{
    std::vector<int>            numit_arr;
    int                         maxit;
    std::vector<int>            refnum_arr;
    int                         refmax;
    std::vector<complex_d_t>    z_arr;
    std::vector<complex_d_t>    pz_arr;
    std::vector<double>         error_arr;
    std::vector<double>         error_ref_arr;
    double                      tol;
    double                      reftol;
    bool                        ref;
    std::string                 method_name;

    // Constructors
    Result() {}

    Result(const std::vector<int>& numit_arr, const int& maxit,
        const std::vector<int>& refnum_arr, const int& refmax,
        const std::vector<complex_d_t>& z_arr,
        const std::vector<complex_d_t>& pz_arr,
        const std::vector<double>& error_arr,
        const std::vector<double>& error_ref_arr,
        const double& tol, const double& reftol, const bool& ref,
        const std::string& method_name) {
        this->numit_arr         = numit_arr;
        this->maxit             = maxit;
        this->refnum_arr        = refnum_arr;
        this->refmax            = refmax;
        this->z_arr             = z_arr;
        this->pz_arr            = pz_arr;
        this->error_arr         = error_arr;
        this->error_ref_arr     = error_ref_arr;
        this->tol               = tol;
        this->reftol            = reftol;
        this->ref               = ref;
        this->method_name       = method_name;
    }

    void print() const {
        if (ref) {
            std::cout << "POLYNOMIAL ROOT FINDER:   "
                << method_name << std::endl;
            std::cout << "MAX ITERATIONS:           " 
                << maxit << std::endl;
            std::cout << "TOLERANCE:                "
                << tol << std::endl;
            std::cout << "REFINEMENT:               True" << std::endl;
            std::cout << "MAX REFINEMENT ITERATIONS:"
                << refmax << std::endl;
            std::cout << "REF TOLERANCE:            "
                << reftol << std::endl;
        }
        else 
            {std::cout << "POLYNOMIAL ROOT FINDER:   "
                << method_name << std::endl;
            std::cout << "MAX ITERATIONS:           " 
                << maxit << std::endl;
            std::cout << "TOLERANCE:                "
                << tol << std::endl;
            std::cout << "REFINEMENT:               False" << std::endl;
        }
        print_table();
    }

    void print_table() const {
        if (ref) {
            std::cout << std::string(129, '-') << std::endl;
            std::cout << "REAL PART\t\tIMAG PART\t\t|FUN VAL|\t\t"
                << "ERROR\t\t\tITER\tERROR_REF\tREF" << std::endl;
            std::cout << std::string(129, '-') << std::endl;
            size_t n = z_arr.size();
            for (size_t k = 0; k < n; k++) {
                std::cout << std::showpos << std::scientific
                << std::setprecision(15) << z_arr[k].real() << "\t"
                << z_arr[k].imag() << "\t" << std::noshowpos
                << std::norm(pz_arr[k]) << "\t" << error_arr[k] << "\t"
                << numit_arr[k] << "\t" << error_ref_arr[k] << "\t" <<
                refnum_arr[k] << std::endl;
            } 
            std::cout << std::string(129, '-') << std::endl;
            return;
        }
        std::cout << std::string(100, '-') << std::endl;
        std::cout << "REAL PART\t\tIMAG PART\t\t|FUN VAL|\t\t"
            << "ERROR\t\t\tITER";
        std::cout << "\n" << std::string(100, '-') << std::endl;
        size_t n = z_arr.size();
        for (size_t k = 0; k < n; k++) {
            std::cout << std::showpos << std::scientific
            << std::setprecision(15) << z_arr[k].real() << "\t"
            << z_arr[k].imag() << "\t" << std::noshowpos << std::norm(pz_arr[k])
            << "\t" << error_arr[k] << "\t" << numit_arr[k] << std::endl;
        } 
        std::cout << std::string(100, '-') << std::endl;
    }
};

/**
 * @brief Evaluates the polynomial using the usual polynomial evaluation
 * @param p polynomial
 * @param z point to evaluate p with
 * @return complex_d_t value of p at z
 */
complex_d_t usualPolyEval(const std::vector<complex_d_t>& p,
    const complex_d_t& z) {
    complex_d_t v = p.at(0), x = 1.;
    for (size_t j = 1; j <= p.size(); j++) {
        x *= z;
        v += p.at(j)*x;
    }
    return v;
}

/**
 * @brief Evaluates the polynomial using nested multiplication and also returns
 * the remainder when p is divided by (x-z)
 * @param p polynomial
 * @param z point to evaluate p at
 * @return std::pair<complex_d_t, std::vector<double> >  value of p at z and
 * remainder of p when divided by x-z
 */
std::pair<complex_d_t, std::vector<complex_d_t> >
    horner(const std::vector<complex_d_t>& p,
    complex_d_t& z) {
    size_t n = p.size() - 1;
    std::vector<complex_d_t> q (n+1, 0.);
    q[n] = p[n];
    for (int k = n-1; k > -1; k--) {
        q[k] = p[k] + q[k+1]*z;
    }
    complex_d_t b0 = q[0];
    q.erase(q.begin());
    return std::make_pair(b0, q);
}

/**
 * @brief Newton-Horner Variation of the Polynomial Deflation Method
 * @param p polynomial
 * @param z initial iterate
 * @return struct Result for the rootspoly via the method
 */
Result newtonHorner(const std::vector<complex_d_t>& p, complex_d_t z,
    param& parameter) {
    timer stopwatch;
    stopwatch.start();
    size_t n = p.size() - 1;
    std::vector<complex_d_t> p_temp = p;

    // Hold roots, function value, iter, iter_ref, err, err_ref
    std::vector<complex_d_t> z_arr;
    std::vector<complex_d_t> pz_arr;
    std::vector<int> numit_arr;
    std::vector<int> refnum_arr;
    std::vector<double> error_arr;
    std::vector<double> error_ref_arr;

    // Polynomial Deflation 
    for (size_t m = 0; m < n; m++) {
        double err_ref = parameter.tol*parameter.reftol;
        int k = 0;
        z = z*complex_d_t(1., 1.);
        double err = parameter.tol + 1.;
        // Linear Polynomial
        if (m == n - 1) {
            k++;
            z = -p_temp[0]/p_temp[1];
            err = std::norm(horner(p_temp, z).first);
        }
        // Polynomial of higher degree
        else {
            while ((err > parameter.tol) && (k < parameter.maxit)) {
                k++;
                complex_d_t z_old = z;
                complex_d_t pz;
                std::vector<complex_d_t> q;
                std::tie(pz, q) = horner(p_temp, z);
                complex_d_t qz = horner(q, z).first;
                if (std::norm(qz) > eps) {
                    z = z - pz/qz;
                    err = max(std::norm(z - z_old), std::norm(pz));
                }
                else {
                    err = 0.;
                }
            }
        }
        numit_arr.push_back(k);
        error_arr.push_back(err);

        complex_d_t pz;
        // Refinement Step
        if (parameter.ref) {
            // Reset iteration counter for refinement
            k = 0;
            complex_d_t z_ref = z;
            err = parameter.tol + 1.;
            while ((err > err_ref) && (k < parameter.refmax)) {
                k++;
                std::vector<complex_d_t> q;
                std::tie(pz, q) = horner(p, z_ref);
                complex_d_t qz = horner(q, z_ref).first;
                if (std::norm(qz) > eps) {
                    complex_d_t z_temp = z_ref;
                    z_ref = z_ref - pz/qz;
                    err = max(std::norm(z_ref - z_temp), std::norm(pz));
                }
                else {
                    err = 0;
                }
            }
            z = z_ref;
            error_ref_arr.push_back(err);
            refnum_arr.push_back(k);
        }
        // Otherwise, compute for pz and update array of funvals
        else {
            pz = horner(p, z).first;
        }
        pz_arr.push_back(pz);
        p_temp = horner(p_temp, z).second;
        z_arr.push_back(z);
    }
    stopwatch.stop();
    return Result(numit_arr, parameter.maxit, refnum_arr, parameter.refmax,
        z_arr, pz_arr, error_arr, error_ref_arr, parameter.tol, parameter.reftol,
        parameter.ref, "NEWTONHORNER");
}

/**
 * @brief Muller-Horner Variation of the Polynomial Deflation Method
 * @param p polynomial
 * @param x0 first point
 * @param x1 second point
 * @param x2 third point
 * @return struct Result for the rootspoly via the method
 */
Result mullerHorner(const std::vector<complex_d_t>& p, const double& x0,
    const double& x1, const double& x2, param& parameter) {
    timer stopwatch;
    stopwatch.start();
    size_t n = p.size() - 1;
    std::vector<complex_d_t> p_temp = p;

    // Hold roots, function value, iter, iter_ref, err, err_ref
    std::vector<complex_d_t> z_arr;
    std::vector<complex_d_t> pz_arr;
    std::vector<int> numit_arr;
    std::vector<int> refnum_arr;
    std::vector<double> error_arr;
    std::vector<double> error_ref_arr;

    complex_d_t z;
    // Polynomial Deflation 
    for (size_t m = 0; m < n; m++) {
        double err_ref = parameter.tol*parameter.reftol;
        int k = 0;
        double err = parameter.tol + 1;
        complex_d_t z0 = x0, z1 = x1, z2 = x2;
        // Linear Polynomial
        if (m == n- 1) {
            k++;
            z = -p_temp[0]/p_temp[1];
        }
        // Polynomial of Higher Degree
        else {
            while ((err > parameter.tol) && (k < parameter.maxit)) {
                k++;
                complex_d_t f0 = horner(p_temp, z0).first;
                complex_d_t f1 = horner(p_temp, z1).first;
                complex_d_t f2 = horner(p_temp, z2).first;
                complex_d_t f01 = (f1 - f0)/(z1 - z0);
                complex_d_t f12 = (f2 - f1)/(z2 - z1);
                complex_d_t f012 = (f12 - f01)/(z2 - z0);
                complex_d_t w = f12 + (z2 - z1)*f012;
                complex_d_t alpha = w*w - 4.*f2*f012;
                complex_d_t d = max(std::norm(w - sqrt(alpha)),
                    std::norm(w + sqrt(alpha)));
                if (std::norm(d) > eps) {
                    z = z2 - 2.*f2/d;
                    err = max(std::norm(z - z2), std::norm(f2));
                }
                else {
                    err = 0.;
                }
                z0 = z1, z1 = z2, z2 = z;
            }
        }
        numit_arr.push_back(k);
        error_arr.push_back(err);

        complex_d_t pz;
        // Refinement Step
        if (parameter.ref) {
            // Reset iteration counter for refinement
            int k = 0;
            complex_d_t z_ref = z;
            double err = parameter.tol + 1.;
            while ((err > err_ref) && (k < parameter.refmax)) {
                k++;
                std::vector<complex_d_t> q;
                std::tie(pz, q) = horner(p, z_ref);
                complex_d_t qz = horner(q, z_ref).first;
                if (std::norm(qz) > eps) {
                    complex_d_t z_temp = z_ref;
                    z_ref = z_ref - pz/qz;
                    err = max(std::norm(z_ref - z_temp), std::norm(pz));
                }
                else {
                    err = 0.;
                }
            }
            z = z_ref;
            error_ref_arr.push_back(err);
            refnum_arr.push_back(k);
        }
        // Otherwise, compute for pz and update array of funvals
        else {
            pz = horner(p, z).first;
        }
        pz_arr.push_back(pz);
        p_temp = horner(p_temp, z).second;
        z_arr.push_back(z);
    }
    stopwatch.stop();
    return Result(numit_arr, parameter.maxit, refnum_arr, parameter.refmax,
        z_arr, pz_arr, error_arr, error_ref_arr, parameter.tol, parameter.reftol,
        parameter.ref, "MULLERHORNER");
}

}

#endif