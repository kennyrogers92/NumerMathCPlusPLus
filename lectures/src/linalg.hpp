/**
 * @file linalg.hpp
 * @author kennyrogers (kenraposas92@gmail.com)
 * @brief C++ Header file for linear algebra system support for C++
 * Real vectors and matrices are typed as vector_t and matrix_t respectively 
 * @version 0.1
 * @date 2024-03-09
 * 
 * @copyright Copyright (c) 2024
 * 
 */
#ifndef LINALG_HPP
#define LINALG_HPP

#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <limits>
#include <string>
#include <tuple>
#include <functional>
#include "timer.hpp"
#include "sgn.hpp"
#include "max.hpp"

namespace linalg {

     // Create new types from the standard
     using vector_t = std::vector<double>;
     using matrix_t = std::vector<vector_t>;	

     // Constant for machine epsilon
     const double eps = std::numeric_limits<double>::epsilon();

     // Output stream for vector_t and matrix_t
     std::ostream& operator<<(std::ostream& os, const vector_t& X) {
          os << "[\t";
          for (double elt : X) {
               os << std::fixed << std::scientific << std::setprecision(5)
                    << elt << "\t";
          }
	     os << "]";
          return os;
     }

     std::ostream& operator<<(std::ostream& os, const matrix_t& X) {
          for (vector_t row : X) os << row << std::endl;
          return os;
     }

     // Arithmetic Operations on Vectors
     vector_t operator+(const vector_t& A, const vector_t& B) {
          size_t n = A.size();
          vector_t C (n, 0);
          for (size_t i = 0; i < n; i++) {
               C[i] = A[i] + B[i];
          }
          return C;
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

     // Arithmetic Operations on Matrices
     matrix_t operator+(const matrix_t& A, const matrix_t& B) {
          size_t n = A.size(), m = A.at(0).size();
          matrix_t C (n, vector_t (m, 0));
          for (size_t i = 0; i < n; i++) {
               for (size_t j = 0; j < m; j++) {
                    C[i][j] = A[i][j] + B[i][j];
               }
          }
          return C;
     }

     matrix_t operator-(const matrix_t& A) {
          size_t n = A.size();
          matrix_t B (n, vector_t (n, 0));
          for (size_t i = 0; i < n; i++) {
               B[i] = -A[i];
          }
          return B;
     }

     matrix_t operator-(const matrix_t& A, const matrix_t& B) {
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

     // Scalar multiplication with a matrix
     matrix_t operator*(const matrix_t& A, const double& x) {
          size_t n = A.size(), m = A.at(0).size();
          matrix_t C (n, vector_t (m, 0));
          for (size_t i = 0; i < n; i++) {
               C[i] = x*A[i];
          }
          return C;
     }

     matrix_t operator*(const double& x, const matrix_t& A) {
          return A*x;
     }

     // Matrix-Vector product. We use SAXPY for this procedure
     vector_t operator*(const matrix_t& A, const vector_t& x) {
          size_t m = A.size(), n = A.at(0).size();
          vector_t y (m, 0);
          for (size_t j = 0; j < n; j++) {
               for (size_t i = 0; i < m; i++) {
                    y[i] += A[i][j]*x[j];
               }
          }
          return y;
     }

     // Product of two matrices A and B with the appropriate sizes
     matrix_t operator*(const matrix_t& A, const matrix_t& B) {
          size_t m = A.size(), n = A.at(0).size(), p = B.at(0).size();
          if (n != B.size()) {
               std::cerr << "Warning: Matrices are not in right dimensions";
               return matrix_t();
          }
          matrix_t C (m, vector_t (p, 0));
          for (size_t k = 0; k < p; k++) {
               for (size_t j = 0; j < n; j++) {
                    for (size_t i = 0; i < m; i++) {
                         C[i].at(k) += A[i][j] * B[j].at(k);
                    }
               }
          }
          return C;
     }

     // Inner product of n-vector x and y
     double dot(const vector_t& x, const vector_t& y) {
          size_t n = x.size();
          double s = 0;
          for (size_t k = 0; k < n; k++) {
               s += x.at(k) * y.at(k);
          }
          return s;
     }

     // Outer product of n-vector x and y
     matrix_t cross(const vector_t& x, const vector_t& y) {
          size_t n = x.size();
          matrix_t A (n, vector_t (n, 0));
          for (size_t i = 0; i < n; i++) {
               for (size_t j = 0; j < n; j++) {
                    A[i][j] = x[i] * y[j];
               }
          }
          return A;
     }

     // Transpose of a n by m matrix A
     matrix_t transpose(const matrix_t& A) {
          size_t n = A.size(), m = A.at(0).size();
          matrix_t T (n, vector_t (m, 0));
          for (size_t i = 0; i < n; i++) {
               for (size_t j = 0; j < m; j++) {
                    T[j][i] = A[i][j];
               }
          }
          return T;
     }

     // Euclidean Norm of a vector
     double Euclid_norm(const vector_t& X) {
          size_t n = X.size();
          double sum = 0;
          for (size_t i = 0; i < n; i++) {
               sum += X[i] * X[i];
          }
          return std::sqrt(sum);
     }

     // Max norm
     double Max_norm(const vector_t& X) {
          size_t n = X.size();
          double t = abs(X[0]);
          for (size_t i = 1; i < n; i++) {
               if (t < abs(X[i])) {
                    t = abs(X[i]);
               }
          }
          return t;
     }

     // Solving Linear Systems

     // Direct Methods
     namespace direct {

     // struct for solution to linear systems by direct methods
     struct Result
     {
          vector_t       x;
          double         residual_max;
          double         elapsed_time;
          std::string    method_name;

          // default constructor
          Result() {}

          // user-defined constructor
          Result (const vector_t &x, const double &residual_max,
               const double &elapsed_time, const std::string &method_name)
          {
               this->x                  = x;
               this->residual_max       = residual_max;
               this->elapsed_time       = elapsed_time;
               this->method_name        = method_name;
          }

          void print_result()
          {
               std::cout << "METHOD:                        "
                    << method_name << std::endl;
               std::cout << std::setprecision(16);
               std::cout << std::fixed;
               std::cout << "APPROXIMATE SOLUTION:          ";
               std::cout << x << std::endl;
               std::cout << std::scientific;
               std::cout << "RESIDUAL MAX:                  "
                    << residual_max << std::endl;
               std::cout << "ELAPSED TIME:                  "
                    << elapsed_time << " seconds" << std::endl;
               std::cout << std::defaultfloat;
          }
     };

     /// Triangular Systems

     // Forward substitution by rows of a lower triangular matrix that solves
     // the equation Lx = b.
     Result forwardsubrow(const matrix_t& L, const vector_t& b) {
          timer stopwatch;
          stopwatch.start();
          size_t n = b.size();
          vector_t x;
          x.push_back(b[0] / L[0][0]);
          for (size_t i = 1; i < n; i++) {
               double s = 0;
               for (size_t j = 0; j < i; j++) {
                    s += L[i][j] * x[j];
               }
               x.push_back((b[i] - s) / L[i][i]);
          }
          stopwatch.stop();
          double residual_max = Max_norm(L*x - b);
          return Result(x, residual_max, stopwatch.get_elapsed_time(),
               "FORWARDSUBROW");
     }

     // Forward substitution by columns of a lower triangular matrix that solves
     // the equation Lx = b. 
     Result forwardsubcol(const matrix_t& L, const vector_t& b) {
          timer stopwatch;
          stopwatch.start();
          vector_t b_old = b;
          size_t n = b.size();
          for (size_t j = 0; j < n - 1; j++) {
               b_old[j] = b_old[j]/L[j][j];
                    for (size_t i = j + 1; i < n; i++) {
                         b_old[i] = b_old[i] - L[i][j]*b_old[j];
                    }
               }
          b_old[n - 1] = b_old[n - 1] / L[n - 1][n - 1];
          stopwatch.stop();
          double residual_max = Max_norm(L*b_old - b);
          return Result(b_old, residual_max, stopwatch.get_elapsed_time(),
               "FORWARDSUBCOL");
     }

     // Backward substitution by rows of an upper triangular matrix that solves
     // the equation Ux = b. 
     Result backwardsubrow(const matrix_t& U, const vector_t& b) {
          timer stopwatch;
          stopwatch.start();
          size_t n = b.size();
          vector_t x (n, 0);
          x[n - 1] = b[n - 1] / U[n - 1][n - 1];
          for (int i = n - 2; i >= 0; i--) {
               double s = 0;
               for (size_t j = i + 1; j < n; j++) {
                    s += U[i][j] * x[j];
               }
               x[i] = (b[i] - s) / U[i][i];
          }
          stopwatch.stop();
          double residual_max = Max_norm(U*x - b);
          return Result(x, residual_max, stopwatch.get_elapsed_time(),
               "BACKWARDSUBROW");
     }

     // Backward substitution by columns of an upper triangular matrix that
     // solves the equation Ux = b. */
     Result backwardsubcol(const matrix_t& U, const vector_t& b) {
          timer stopwatch;
          stopwatch.start();
          vector_t b_old = b;
          size_t n = b.size();
          for (size_t j = n - 1; j > 0; j--) {
               b_old[j] = b_old[j] / U[j][j];
               for (size_t i = 0; i < j; i++) {
                    b_old[i] = b_old[i] - U[i][j] * b_old[j];
               }
          }
          b_old[0] = b_old[0] / U[0][0];
          stopwatch.stop();
          double residual_max = Max_norm(U*b_old - b);
          return Result(b_old, residual_max, stopwatch.get_elapsed_time(),
               "BACKWARDSUBCOL");
     }

     /// LU Factorization

     // Single-precision version of LU factorization of A.
     matrix_t LUkji(const matrix_t& A) {
          size_t n = A.size();
          matrix_t B = A;
          for (size_t k = 0; k < n; k++) {
               for (size_t j = k + 1; j < n; j++) {
                    B[j][k] = B[j][k]/B[k][k];
               }
               for (size_t j = k + 1; j < n; j++) {
                    for (size_t i = k + 1; i < n; i++) {
                         B[i][j] = B[i][j]  - B[i][k]*A[k][j];
                    }
               }
          }
          return B;
     }

     // Generalized single-precision version of LU factorization of A
     matrix_t LUjki(const matrix_t& A) {
          size_t n = A.size();
          matrix_t B = A;
          for (size_t j = 0; j < n; j++) {
               for (size_t k = 0; k < j; k++) {
                    for (size_t i = k + 1; i < n; i++) {
                         B[i][j] = B[i][j] - B[i][k] * B[k][j];
                    }
               }
               for (size_t k = j + 1; k < n; k++) {
                    B[k][j] = B[k][j] / B[j][j];
               }
          }
          return B;
     }

     // Doolittle method of LU factorization of A 
     matrix_t LUijk(const matrix_t& A) {
          size_t n = A.size();
          matrix_t B = A;
          for (size_t j = 1; j < n; j++) {
               B[j][0] = B[j][0] / B[0][0];
          }
          for (size_t i = 1; i < n; i++) {
               for (size_t j = i; j < n; j++) {
                    double s = 0;
                    for (size_t k = 0; k < i; k++) {
                         s += B[i][k] * B[k][j];
                    }
                    B[i][j] = B[i][j] - s;
               }
               for (size_t j = i + 1; j < n; j++) {
                    double s = 0;
                    for (size_t k = 0; k < i; k++) {
                         s += B[j][k] * B[k][i];
                    }
                    B[j][i] = (B[j][i] - s) / B[i][i];
               }
          }
          return B;
     }

     // LU factorization with partial pivoting followed Doolittle Method
     std::tuple<matrix_t, std::vector<int> > LUpartial(const matrix_t& A,
          double eps=1) {
          size_t n = A.size();
          matrix_t B = A;
          std::vector<int> p (n);
          for (size_t i = 0; i < n; i++) {
               p[i] = i+1;
          }
          for (size_t i = 0; i < n; i++) {
               if (abs(B[i][i]) < eps) {
                    // Find an index idx where |A[idx][i]| = max_1<=j<=n(|A[i][j]|)
                    vector_t pivot;
                    for (size_t j = i; j < n; j++) {
                         pivot.push_back(abs(B[i][j]));
                    }
                    long idx;
                    max(pivot, idx);
                    // Interchange idxth and ith rows of A as well as entries of p
                    std::swap(B[i], B[idx]);
                    std::swap(p[i], p[idx]);
               }
               if (i == 0) {
                    for (size_t j = 1; j < n; j++) {
                         B[j][0] = B[j][0] / B[0][0];
                    }
               }
               else {
                    for (size_t j = i; j < n; j++) {
                         double s = 0;
                         for (size_t k = 0; k < i; k++) {
                              s += B[i][k] * B[k][j];
                         }
                         B[i][j] = B[i][j] - s;
                    }
                    for (size_t j = i + 1; j < n; j++) {
                         double s = 0;
                         for (size_t k = 0; k < i; k++) {
                              s += B[j][k] * B[k][i];
                         }
                         B[j][i] = (B[j][i] - s) / B[i][i];
                    }
               }
          }
          return std::make_tuple(B, p);
     }

     // LU factorization with complete pivoting followed Doolittle Method
     std::tuple<matrix_t, std::vector<int>, std::vector<int> > LUcomplete(
          const matrix_t& A, double eps=1.) {     
          size_t n = A.size();
          matrix_t B = A;
          std::vector<int> p (n), q (n);
          for (size_t i = 0; i < n; i++) {
               p[i] = i+1;
               q[i] = i+1;
          }
          for (size_t i = 0; i < n; i++) {
               if (abs(B[i][i]) < eps) {
                    // Find an index r,s where |A[r][s]| = max_1<=k,j<=n(|A[k][j]|)
                    vector_t pivot;
                    for (size_t k = i; k < n; k++) {
                         for (size_t j = i; j < n; j++) {
                              pivot.push_back(abs(B[k][j]));
                         }
                    }
                    long idx;
                    max(pivot, idx);
                    int r = idx / (n-i);
                    int s = idx % (n-i);
                    // Interchange rows and cols of A as well as entries of p
                    std::swap(B[i], B[r]);
                    B = transpose(B);
                    std::swap(B[i], B[s]);
                    B = transpose(B);
                    std::swap(p[i], p[r]);
                    std::swap(q[i], q[s]);
               }
               if (i == 0) {
                    for (size_t j = 1; j < n; j++) {
                         B[j][0] = B[j][0] / B[0][0];
                    }
               }
               else {
                    for (size_t j = i; j < n; j++) {
                         double s = 0;
                         for (size_t k = 0; k < i; k++) {
                              s += B[i][k] * B[k][j];
                         }
                         B[i][j] = B[i][j] - s;
                    }
                    for (size_t j = i + 1; j < n; j++) {
                         double s = 0;
                         for (size_t k = 0; k < i; k++) {
                              s += B[j][k] * B[k][i];
                         }
                         B[j][i] = (B[j][i] - s) / B[i][i];
                    }
               }
          }
          return std::make_tuple(B, p, q);
     }

     // Returns the LU factorization of A, where L and U are stored
     // as separate matrices. Note: This procudure takes as argument
     // the return of one of the three LU Factorization procedures 
     std::pair<matrix_t, matrix_t> getLU(const matrix_t& A) {
          size_t n = A.size();
          matrix_t L (n, vector_t (n, 0)), U (n, vector_t (n, 0));
          for (size_t i = 0; i < n; i++) {
               L[i][i] = 1;
               for (size_t j = 0; j < i; j++) {
                    L[i][j] = A[i][j];
               }
               for (size_t j = i; j < n; j++) {
                    U[i][j] = A[i][j];
               }
          }
          return std::pair<matrix_t, matrix_t>(L, U);
     }

     // Solves the equation Ax = b through LU factorization
     Result solveLU(const matrix_t& A, const vector_t& b,
          std::string method="ijk") {
          timer stopwatch;
          stopwatch.start();
          matrix_t B = A;
          std::string method_name;
          if (method == "kji") {
               B = LUkji(A);
               method_name = "LUSolve - KJI";
          }
          else if (method == "jki") {
               B = LUjki(A);
               method_name = "LUSolve - JKI";
          }
          else if (method == "ijk") {
               B = LUijk(A);
               method_name = "LUSolve - IJK";
          }
          matrix_t L, U;
          std::tie(L, U) = getLU(B);
          vector_t y = forwardsubcol(L, b).x;
          vector_t x = backwardsubcol(U, y).x;
          stopwatch.stop();
          double residual_max = Max_norm(A*x - b);
          return Result(x, residual_max, stopwatch.get_elapsed_time(),
               method_name); 
     }

     // Solves the equation Ax = b through LU factorization with pivoting
     Result solveLUPivot(const matrix_t& A, const vector_t& b, double eps=1.,
          std::string method="partial") {
          timer stopwatch;
          stopwatch.start();
          matrix_t B = A;
          // Containers for p, q
          std::vector<int> p, q;
          std::string method_name;
          if (method == "partial") {
               std::tie(B, p) = LUpartial(A, eps);
               method_name = "LUPartialPivotSolve";
          }
          if (method == "complete") {
               std::tie(B, p, q) = LUcomplete(A, eps);
               method_name = "LUCompletePivotSolve";
          }
          
          matrix_t L, U;
          std::tie(L, U) = getLU(B);
          // Permute b according to p

          // Proceed as before
          vector_t y = forwardsubcol(L, b).x;
          vector_t x = backwardsubcol(U, y).x;
          // Permute x if pivoting is complete
          stopwatch.stop();
          double residual_max = Max_norm(A*x - b);
          return Result(x, residual_max, stopwatch.get_elapsed_time(),
               method_name); 
     }

     // Refinement of solution x` obtained by direct solvers
     // by solving the system Az = b - Ax` and taking x = x` + z
     // By default we use LU factorization with partial pivoting 
     // for the direct solvers. One can choose which method to use
     // by modifying lines 5 and 8 internally
     Result iterativerefinement(const matrix_t& A, const vector_t& b,
          double tol=2.*eps, int maxit=10) {
          timer stopwatch;
          stopwatch.start();
          double err = tol + 1;
          int k = 0;
          vector_t x = solveLU(A, b).x;
          while ((err > tol) && (k < maxit)) {
               vector_t r = b - A*x;
               vector_t z = solveLU(A, r).x;
               x = x + z;
               err = Euclid_norm(z)/Euclid_norm(x);
          }
          stopwatch.stop();
          double residual_max = Max_norm(A*x - b);
          return Result(x, residual_max, stopwatch.get_elapsed_time(),
               "IterativeRef");
     }

     // QR factorization of A by the Gram-Schidmt Orthogonalization 
     std::pair<matrix_t, matrix_t> ModifiedGramSchmidt(const matrix_t& A) {
          size_t m = A.size();
          size_t n = A.at(0).size();
          matrix_t B = A;
          matrix_t Q (m, vector_t (n, 0)), R = B;
          double q = Euclid_norm(transpose(B)[0]);
          for (size_t i = 0; i < m; i++) {
               Q[i][0] = B[i][0] / q;
          }
          for (size_t k = 1; k < n; k++) {
               for (size_t j = 0; j < k; j++) {
                    double s = 0;
                    for (size_t i = 0; i < m; i++) {
                         s += Q[i][j] * B[i][k];
                    }
                    for (size_t i = 0; i < m; i++) {
                         B[i][k] = B[i][k] - s * Q[i][j];
                    }
               }
               q = Euclid_norm(transpose(B)[k]);
               for (size_t i = 0; i < m; i++) {
                    Q[i][k] = B[i][k] / q;
               }
          }
          R = transpose(Q)*R;
          return std::pair(Q, R);
     }

     // QR factorization of A by the Gram-Schidmt Orthogonalization that
     // guarantees linear independece 
     std::pair<matrix_t, matrix_t> ModifiedGramSchmidt2(const matrix_t& A) {
          size_t m = A.size();
          size_t n = A.at(0).size();
          matrix_t B = A;
          matrix_t Q (m, vector_t (n, 0)), R (n, vector_t (n, 0));
          for (size_t k = 0; k < n; k++) {
               R[k][k] = Euclid_norm(transpose(B)[k]);
               for (size_t j = 0; j < m; j++) {
                    Q[j][k] = B[j][k] / R[k][k];
               }
               for (size_t j = k + 1; j < n; j++) {
                    for (size_t i = 0; i < m; i++) {
                         R[k][j] = R[k][j] + Q[i][k]*B[i][j];
                    }
                    for (size_t i = 0; i < m; i++) {
                         B[i][j] = B[i][j] - R[k][j] * Q[i][k];
                    }
               }     
          }
          return std::pair<matrix_t, matrix_t>(Q, R);
     }

     // Solves the equation Ax = b through QR factorization. By default
     // uses ModifiedGSO2
     Result solveQR(const matrix_t& A, const vector_t& b,
          int method=1) {
          timer stopwatch;
          stopwatch.start();
          matrix_t Q, R;
          std::string method_name;
          if (method == 0) {
               std::tie(Q, R) = ModifiedGramSchmidt(A);
               method_name = "QRSolve - GSO";
          }
          else {
               std::tie(Q, R)  = ModifiedGramSchmidt2(A);
               method_name = "QRSolve - GSO2";
          }
          vector_t x = backwardsubcol(R, transpose(Q)*b).x;
          stopwatch.stop();
          double residual_max = Max_norm(A*x - b);
          return Result(x, residual_max, stopwatch.get_elapsed_time(),
               method_name);
     }

     // LU factorization of a symmetric positive definite (SPD) A where U = L^T
     matrix_t Cholesky(const matrix_t& A) {
          size_t n = A.size();
          matrix_t L (n, vector_t (n, 0));
          L[0][0] = std::sqrt(A[0][0]);
          for (size_t i = 1; i < n; i++) {
               for (size_t j = 0; j < i; j++) {
                    double s = 0;
                    for (size_t k = 0; k < j; k++) {
                         s += L[i][k] * L[j][k];
                    }
                    L[i][j] = (A[i][j] - s) / L[j][j];
               }
               double s = 0;
               for (size_t j = 0; j < n; j++) {
                    s += L[i][j] * L[i][j];
               }
               L[i][i] = std::sqrt(A[i][i] - s);
          }
          return L;
     }

     // Solves the equation Ax = b through Cholesky factorization
     Result solveChol(const matrix_t& A, const vector_t& b) {
          timer stopwatch;
          stopwatch.start();
          matrix_t L = Cholesky(A);
          matrix_t L_T = transpose(L);
          vector_t y = forwardsubcol(L, b).x;
          vector_t x = backwardsubcol(L_T, y).x;
          stopwatch.stop();
          double residual_max = Max_norm(A*x - b);
          return Result(x, residual_max, stopwatch.get_elapsed_time(),
               "CholSolve");
     }

     }

     // Iterative Methods
     namespace iterative {

     // Struct for solution to linear systems by iterative methods
     struct Result
     {
          vector_t       x;
          int            numit;
          int            maxit;
          double         rel_err;
          double         tol;
          double         elapsed_time;
          std::string    method_name;
          std::string    termination_flag;

          // default constructor
          Result() {}

          // user-defined constructor
          Result (const vector_t &x, const int &numit, const int &maxit, 
               const double &rel_err, const double tol,
               const double &elapsed_time, const std::string &method_name,
               const std::string &termination_flag)
          {
               this->x                  = x;
               this->numit              = numit;
               this->maxit              = maxit;
               this->rel_err       = rel_err;
               this->tol                = tol;
               this->elapsed_time       = elapsed_time;
               this->method_name        = method_name;
               this->termination_flag   = termination_flag;
          }

          void print_result()
          {
               std::cout << "METHOD:                        "
                    << method_name << std::endl;
               std::cout << std::setprecision(16);
               std::cout << std::fixed;
               std::cout << "APPROXIMATE SOLUTION:          ";
               std::cout << x << std::endl;
               std::cout << "TERMINATION:                   "
                    << termination_flag << std::endl;
               std::cout << std::scientific;
               std::cout << "RESIDUAL MAX:                  "
                    << rel_err << std::endl;
               std::cout << "TOLERANCE:                     "
                    << tol << std::endl;
               std::cout << "NUM ITERATIONS:                "
                    << numit << std::endl;
               std::cout << "MAX ITERATIONS:                "
                    << maxit << std::endl;
               std::cout << "ELAPSED TIME:                  "
                    << elapsed_time << " seconds" << std::endl;
               std::cout << std::defaultfloat;
          }
     };

     // Struct for parameters in Iterative methods
     struct param
     {
          double    tol = 1e-12;
          int       maxit = 1e4;
          double    omega = 0.5;
     };

     // Jacobi Over-relaxation Method  
     Result JOR(const matrix_t& A, const vector_t& b, const vector_t& x,
          param &parameter) {
          timer stopwatch;
          stopwatch.start();
          std::string term_flag = "Success";
          size_t n = b.size();

          // Pivoting step

          vector_t y = x;
          double b_norm = Euclid_norm(b);
          double err = Euclid_norm(b - A*y) / b_norm;
          int k = 0;

          // main while loop
          while ((err > parameter.tol) && (k < parameter.maxit)) {
               vector_t y_old;
               double s;
               for (size_t i = 0; i < n; i++) {
                    y_old = y, s = 0;
                    for (size_t j = 0; j < n; j++) {
                         if (j != i) {
                              s += A[i][j] * y_old[j];
                         }
                    }
                    y[i] = parameter.omega*(b[i] - s) / A[i][i]
                         + (1 - parameter.omega)*y_old[i];
               }
               err = Euclid_norm(b - A*y) / b_norm;
               k++;
          }
          if ((err > parameter.tol) && (k == parameter.maxit)) {
               term_flag = "Fail";
          }
          stopwatch.stop();
          double rel_err = Max_norm(A*y - b);
          return Result(y, k, parameter.maxit, rel_err, parameter.tol,
               stopwatch.get_elapsed_time(), "JOR", term_flag);
     }

     // Successive Over-relaxation Method
     Result SOR(const matrix_t& A, const vector_t& b, const vector_t x,
          param &parameter) {
          timer stopwatch;
          stopwatch.start();
          std::string term_flag = "Success";
          size_t n = b.size();

          // Pivoting step

          vector_t y = x;
          double b_norm = Euclid_norm(b);
          double err = Euclid_norm(b - A*y) / b_norm;
          int k = 0;

          // main while loop
          while ((err > parameter.tol) && (k < parameter.maxit)) {
               vector_t y_old;
               double s;
               for (size_t i = 0; i < n; i++) {
                    y_old = y, s = 0;
                    for (size_t j = 0; j < i; j++) {
                         s = s + A[i][j] * y[j];
                    }
                    for (size_t j = i + 1; j < n; j++) {
                         s = s + A[i][j] * y_old[j];
                    }
                    y[i] = parameter.omega * (b[i] - s) / A[i][i]
                         + (1 - parameter.omega) * y_old[i];
               }
               err = Euclid_norm(b - A*y) / b_norm;
               k++;
          }
          if ((err > parameter.tol) && (k == parameter.maxit)) {
               term_flag = "Fail";
          }
          stopwatch.stop();
          double rel_err = Max_norm(A*y - b);
          return Result(y, k, parameter.maxit, rel_err, parameter.tol,
               stopwatch.get_elapsed_time(), "SOR", term_flag);
     }

     // Backward Symmetric Over-relaxation method
     Result BSOR(const matrix_t& A, const vector_t& b, const vector_t x,
          param &parameter) {
          timer stopwatch;
          stopwatch.start();
          std::string term_flag = "Success";
          size_t n = b.size();

          // Pivoting step

          vector_t y = x;
          double b_norm = Euclid_norm(b);
          double err = Euclid_norm(b - A*y) / b_norm;
          int k = 0;

          // main while loop
          while ((err > parameter.tol) && (k < parameter.maxit)) {
               vector_t y_old;
               double s;
               for (size_t i = 0; i < n; i++) {
                    y_old = y, s = 0;
                    for (size_t j = 0; j < i; j++) {
                         s = s + A[i][j] * y_old[j];
                    }
                    for (size_t j = i + 1; j < n; j++) {
                         s = s + A[i][j] * y[j];
                    }
                    y[i] = parameter.omega * (b[i] - s) / A[i][i]
                         + (1 - parameter.omega) * y_old[i];
               }
               err = Euclid_norm(b - A*y) / b_norm;
               k++;
          }
          if ((err > parameter.tol) && (k == parameter.maxit)) {
               term_flag = "Fail";
          }
          stopwatch.stop();
          double rel_err = Max_norm(A*y - b);
          return Result(y, k, parameter.maxit, rel_err, parameter.tol,
               stopwatch.get_elapsed_time(), "BSOR", term_flag);
     }

     // Symmetric Successive Over-relaxation Method
     Result SSOR(const matrix_t& A, const vector_t& b, const vector_t x,
          param &parameter) {
          timer stopwatch;
          stopwatch.start();
          std::string term_flag = "Success";
          size_t n = b.size();

          // Pivoting step

          vector_t y = x;
          double b_norm = Euclid_norm(b);
          double err = Euclid_norm(b - A*y) / b_norm;
          int k = 0;

          // main while loop
          while ((err > parameter.tol) && (k < parameter.maxit)) {
               vector_t z = y, y_old = y;
               double r, s;
               for (size_t i = 0; i < n; i++) {
                    r = 0, s = 0;
                    for (size_t j = 0; j < i; j++) {
                         r = r + A[i][j] * z[j];
                         s = s + A[i][j] * z[j];
                    }
                    for (size_t j = i + 1; j < n; j++) {
                         r = r + A[i][j] * y_old[j];
                         s = s + A[i][j] * y[j];
                    }
                    z[i] = parameter.omega * (b[i] - r) / A[i][i]
                         + (1 - parameter.omega) * y_old[i];
                    y[i] = parameter.omega * (b[i] - s) / A[i][i]
                         + (1 - parameter.omega) * z[i];
               }
               err = Euclid_norm(b - A*y) / b_norm;
               k++;
          }
          if ((err > parameter.tol) && (k == parameter.maxit)) {
               term_flag = "Fail";
          }
          stopwatch.stop();
          double rel_err = Max_norm(A*y - b);
          return Result(y, k, parameter.maxit, rel_err, parameter.tol,
               stopwatch.get_elapsed_time(), "SSOR", term_flag);
     }

     // Optimization-based methods for symmetric positive definite matrices

     // Steepest Descent Method 
     Result steepestdescent(const matrix_t& A, const vector_t& b,
          const vector_t& x, param &parameter) {
          timer stopwatch;
          stopwatch.start();
          std::string term_flag = "Success";
          vector_t y = x;
          vector_t r = b - A*y;
          double rho = dot(r, r);
          int k = 0;

          // main while loop
          while (k < parameter.maxit) {
               double s = dot(r, A*r);
               double alpha = rho / s;
               y = y + alpha*r;
               r = b - A*y;
               rho = dot(r, r);
               k++;
               if (std::sqrt(rho) <= parameter.tol) {
                    break;
               }
          }
          if ((std::sqrt(rho) > parameter.tol) && (k == parameter.maxit)) {
               term_flag = "Fail";
          }
          stopwatch.stop();
          double rel_err = Max_norm(A*y - b);
          return Result(y, k, parameter.maxit, rel_err, parameter.tol,
               stopwatch.get_elapsed_time(), "SteepestDescent", term_flag);
     }

     // Conjugate Gradient Method
     Result conjugategradient(const matrix_t& A, const vector_t& b,
          const vector_t& x, param &parameter) {
          timer stopwatch;
          stopwatch.start();
          std::string term_flag = "Success";
          vector_t y = x;
          vector_t r = b - A*y;
          vector_t d = r;
          double rho = dot(r, r);
          int k = 0;

          // main while loop
          while (k < parameter.maxit) {
               vector_t w = A*d;
               double alpha = rho / dot(d, w);
               y = y + alpha*d;
               r = r - alpha*w;
               double rho_old = rho;
               rho = dot(r, r);
               k++;
               if (std::sqrt(rho) <= parameter.tol) {
                    break;
               }
               double beta = rho_old / rho;
               d = r + beta*d;
          }
          if ((std::sqrt(rho) > parameter.tol) && (k == parameter.maxit)) {
               term_flag = "Fail";
          }
          stopwatch.stop();
          double rel_err = Max_norm(A*y - b);
          return Result(y, k, parameter.maxit, rel_err, parameter.tol,
               stopwatch.get_elapsed_time(), "ConjugateGradient", term_flag);
     }

     // Conjugate Gradient Normal Residual Method
     Result cgnormalresidual(const matrix_t& A, const vector_t& b,
          const vector_t& x, param &parameter) {
          timer stopwatch;
          stopwatch.start();
          std::string term_flag = "Success";
          vector_t y = x;
          vector_t r = b - A*y;
          vector_t d = transpose(A)*r;
          vector_t q = d;
          double sigma = dot(q, q);
          int k = 0;

          // main while loop
          while (k < parameter.maxit) {
               vector_t w = A*d;
               double alpha = sigma / dot(w, w);
               y = y + alpha*d;
               r = r - alpha*w;
               k++;
               if (std::sqrt(dot(r, r)) <= parameter.tol) {
                    break;
               }
               double sigma_old = sigma;
               q = transpose(A)*r;
               sigma = dot(q, q);
               double beta = sigma_old / sigma;
               d = q + beta*d;
          }
          if ((std::sqrt(dot(r, r)) > parameter.tol) && (k == parameter.maxit)) {
               term_flag = "Fail";
          }
          stopwatch.stop();
          double rel_err = Max_norm(A*y - b);
          return Result(y, k, parameter.maxit, rel_err, parameter.tol,
               stopwatch.get_elapsed_time(), "CGNormalResidual", term_flag);
     }

     // CG Residual Method
     Result cgresidual(const matrix_t& A, const vector_t& b, const vector_t& x,
          param &parameter) {
          timer stopwatch;
          stopwatch.start();
          std::string term_flag = "Success";
          vector_t y = x;
          vector_t r = b - A*y;
          vector_t d = r;
          int k = 0;

          // main while loop
          while (k < parameter.maxit) {
               vector_t w = A*d;
               double omega = dot(w, w);
               double alpha = dot(r, w) / omega;
               y = y + alpha*d;
               r = r - alpha*w;
               k++;
               if (std::sqrt(dot(r, r)) <= parameter.tol) {
                    break;
               }
               double beta = -dot(w, A*r) / omega;
               d = r + beta*d;
          }
          if (std::sqrt(dot(r, r) > parameter.tol) && (k == parameter.maxit)) {
               term_flag = "Fail";
          }
          stopwatch.stop();
          double rel_err = Max_norm(A*y - b);
          return Result(y, k, parameter.maxit, rel_err, parameter.tol,
               stopwatch.get_elapsed_time(), "CGResidual", term_flag);
     }

     }

     // Nonlinear solvers

     namespace nonlinear {

     // Type for a function with input and output having the type vector and matrix
     using VectVarFunction = std::function<vector_t(vector_t)>;
     using MatVarFunction = std::function<matrix_t(vector_t)>;

     // Struct for parameters in solvers of nonlinear system
     struct param
     {
          double    tol=eps;
          int       maxit=100;
     };

     // Struct for Result of nonlinear system solvers
     struct Result
     {
          int       numit;
          int       maxit;
          vector_t  x;
          double    funval_norm;
          double    error;
          double    tol;
          double    elapsed_time;
          std::string    method_name;
          std::string    termination_flag;
          
          // Default Constructor
          Result() {}

          // user-defined constructor
          Result(const int& numit, const int& maxit, const vector_t& x,
               const double& funval_norm, const double& error, const double& tol,
               const double& elapsed_time, const std::string& method_name,
               const std::string termination_flag) {
               this->numit              = numit;
               this->maxit              = maxit;
               this->x                  = x;
               this->funval_norm        = funval_norm;
               this->error              = error;
               this->tol                = tol;
               this->elapsed_time       = elapsed_time;
               this->method_name        = method_name;
               this->termination_flag   = termination_flag;
          }

          void print() {
               std::cout << "METHOD:                        "
                    << method_name << std::endl;
               std::cout << std::setprecision(16);
               std::cout << std::fixed;
               std::cout << "TERMINATION:                   "
                    << termination_flag << std::endl;
               std::cout << std::scientific;
               std::cout << "FUNVAL NORM:                   "
                    << funval_norm << std::endl;
               std::cout << "ERROR:                         "
                    << error << std::endl;
               std::cout << "TOLERANCE:                     "
                    << tol << std::endl;
               std::cout << "NUM ITERATIONS:                "
                    << numit << std::endl;
               std::cout << "MAX ITERATIONS:                "
                    << maxit << std::endl;
               std::cout << "ELAPSED TIME:                  "
                    << elapsed_time << " seconds" << std::endl;
               std::cout << "APPROXIMATE SOLUTION:          \n\t";
               std::cout << x << std::endl;
          }
     };

     // Generalization of the Newton Method to solve linear systems of eqns
     Result newton(VectVarFunction f, MatVarFunction Df, const vector_t& x,
          param &parameter) {
          timer stopwatch;
          stopwatch.start();
          std::string term_flag = "Success";
          vector_t y = x;
          double err = parameter.tol + 1.;
          int k = 0;
          while ((err > parameter.tol) && (k < parameter.maxit)) {
               // Initialize Jacobian Matrix J
               vector_t delta = direct::solveLU(Df(y), -f(y)).x;
               y = y + delta;
               err = Max_norm(delta);
               k++;
          }
          if ((err > parameter.tol) && (k < parameter.maxit)) {
               term_flag = "Fail";
          }
          stopwatch.stop();
          double fx_norm = Euclid_norm(f(y));
          return Result(k, parameter.maxit, y, fx_norm, err, parameter.tol,
               stopwatch.get_elapsed_time(), "NEWTON", term_flag);
     }

     // Broyden Method that approximates Jacobian
     Result broyden(VectVarFunction f, const matrix_t& Q, const vector_t& x,
          param &parameter) {
          timer stopwatch;
          stopwatch.start();
          std::string term_flag = "Success";
          matrix_t R = Q;
          vector_t y = x;
          double err = parameter.tol + 1.;
          int k = 0;
          while (k < parameter.maxit) {
               vector_t delta = direct::solveLU(R, -f(y)).x;
               y = y + delta;
               err = Max_norm(delta);
               if (err > parameter.tol) {
                    R = R + cross(f(y), delta)*(1./dot(delta, delta));
               }
               else {
                    break;
               }
               k++;
          }
          if ((err > parameter.tol) && (k < parameter.maxit)) {
               term_flag = "Fail";
          }
          stopwatch.stop();
          double fx_norm = Euclid_norm(f(y));
          return Result(k, parameter.maxit, y, fx_norm, err, parameter.tol,
               stopwatch.get_elapsed_time(), "BROYDEN", term_flag);
     }

     // FixPoint Method adoption to nonlinear system
     Result fixpoint(VectVarFunction g, const vector_t &x, param &parameter) {
          timer stopwatch;
          stopwatch.start();
          std::string term_flag = "Success";
          vector_t y = x;
          double err = parameter.tol + 1.;
          int k = 0;
          while ((err > parameter.tol) && (k < parameter.maxit)) {
               vector_t y_old = y;
               y = g(y);
               err = Euclid_norm(y - y_old);
               k++;
          }
          if ((err > parameter.tol) && (k == parameter.maxit)) {
               term_flag = "Fail";
          }
          stopwatch.stop();
          double fx_norm = Euclid_norm(g(y));
          return Result(k, parameter.maxit, y, fx_norm, err, parameter.tol,
               stopwatch.get_elapsed_time(), "FIXPOINT", term_flag);
     }

     }

     using UniVarFunction = double(double);
     // Returns the coefficients of interpolating polynomial of f given
     // interpolating nodes. Uses the method of Undetermined coefficients
     // This method tries several linear system solver until one successfully
     // solves for a solution
     vector_t polyinterp(UniVarFunction& f, const vector_t& x) {
          size_t n = x.size();
          // Construct the Vandermonde Matrix
          matrix_t V (n, vector_t (n));
          for (size_t i = 0; i < n; i++) {
               double v = 1.;
               for (size_t j = 0; j < n; j++) {
                    V[i][j] = v;
                    v = v*x[i];
               }
          }
          // Construct the coefficient vector b of images of f under x
          vector_t b (n);
          for (size_t i = 0; i < n; i++) {
               b[i] = f(x[i]);
          }
          vector_t p = direct::solveLU(V, b).x;
          return p;
     }
}

#endif
