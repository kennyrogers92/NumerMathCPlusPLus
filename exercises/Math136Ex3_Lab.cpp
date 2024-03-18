#include <iostream>
#include <iomanip>
#include <string>
#include <cmath>
#include <vector>
#include <utility>
#include "../lectures/src/linalg.hpp"

using namespace linalg;

// Define An(n), bn(n)
matrix_t A(int n) {
    matrix_t B (n, vector_t (n, 0));
    for (int i = 0; i < n; i++) {
        B[i][i] = 2.;
    }
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if ((i == j + 1) || (j == i + 1)) {
                B[i][j] = -1;
            }
        }
    }
    return B*(n + 1.)*(n + 1.);
}

vector_t b(int n) {
    std::vector<double> c;
    for (int i = 1; i < n + 1; i++) {
        c.push_back(sin(i*M_PI/(n + 1.)));
    }
    return c*M_PI*M_PI;
}

// In this program, we test all direct, iterative and gradiet methods 
// in solving linear systems
int main() {
    timer stopwatch;
    stopwatch.start();
    int n = 250;
    matrix_t An = A(n);
    vector_t bn = b(n);

    // Direct Methods
    std::vector<direct::Result> result0;
    result0.push_back(direct::solveLU(An, bn, "kji"));
    result0.push_back(direct::solveLU(An, bn, "jki"));
    result0.push_back(direct::solveLU(An, bn, "ijk"));
    result0.push_back(direct::solveLUPivot(An, bn, 1., "partial"));
    // result0.push_back(direct::solveLUPivot(An, bn, 1., "complete"));
    // result0.push_back(direct::iterativerefinement(An, bn));
    result0.push_back(direct::solveQR(An, bn, 0));
    result0.push_back(direct::solveQR(An, bn));
    result0.push_back(direct::solveChol(An, bn));

    int m = result0.size();

    // Print results
    std::cout << "> DIRECT METHODS" << std::endl;
    std::cout << std::string(87, '-') << std::endl;
    std::cout << "METHOD\t\t\tRESIDUAL MAX NORM\t" <<
        "ERROR NORM\t\tELAPSED TIME" << std::endl;
    std::cout << std::string(87, '-') << std::endl;
    for (int i = 0; i < m; i++) {
        double err = linalg::Max_norm(result0[i].x - std::pow(M_PI, -2)*bn)
            /(n+1);
        std::cout << result0[i].method_name << "\t";
        if (result0[i].method_name.size() < 15) {
            std::cout << "\t";
        }
        std::cout << std::fixed << std::scientific << std::setprecision(15)
            << result0[i].residual_max << "\t" << err << "\t"
            << std::setprecision(5) << result0[i].elapsed_time << " sec\n";
    }
    std::cout << std::string(87, '-') << std::endl;

    // Iterative methods
    iterative::param parameters;
    std::vector<iterative::Result> result1;
    vector_t xn (n);
    for (int i = 0; i < 3; i++) {
        parameters.omega = 0.5*(i+1);
        result1.push_back(iterative::JOR(An, bn, xn, parameters));
        result1.push_back(iterative::SOR(An, bn, xn, parameters));
        result1.push_back(iterative::BSOR(An, bn, xn, parameters));
        result1.push_back(iterative::SSOR(An, bn, xn, parameters));
    }
    
    m = result1.size();

    // Print results
    std::cout << "\n> ITERATIVE METHODS" << std::endl;
    std::cout << std::string(63, '-') << std::endl;
    std::cout << "METHOD\tOMEGA\tRELATIVE ERROR\t\t" <<
        "NUMIT\tELAPSED TIME" << std::endl;
    std::cout << std::string(63, '-') << std::endl;
    for (int i = 0; i < m; i++) {
        std::cout << result1[i].method_name << "\t"
            << std::defaultfloat << std::setprecision(2) << ((i)/4+1)*0.5 << "\t"
            << std::fixed << std::scientific << std::setprecision(15)
            << result1[i].rel_err << "\t" << result1[i].numit << "\t"
            << std::setprecision(5) << result1[i].elapsed_time << " sec\n";
    }
    std::cout << std::string(63, '-') << std::endl;

    // Gradient methods
    result1.clear();
    parameters.tol = 2.*eps;
    parameters.maxit = 500;
    result1.push_back(iterative::steepestdescent(An, bn, xn, parameters));
    result1.push_back(iterative::conjugategradient(An, bn, xn, parameters));
    result1.push_back(iterative::cgresidual(An, bn, xn, parameters));
    result1.push_back(iterative::cgnormalresidual(An, bn, xn, parameters));
    
    m = result1.size();

    // Print results
    std::cout << "\n> GRADIENT METHODS" << std::endl;
    std::cout << std::string(71, '-') << std::endl;
    std::cout << "METHOD\t\t\tRELATIVE ERROR\t\tNUMIT\tELAPSED TIME" << std::endl;
    std::cout << std::string(71, '-') << std::endl;
    for (int i = 0; i < m; i++) {
        std::cout << result1[i].method_name << "\t";
        if (result1[i].method_name.size() < 16) {
            std::cout << "\t";
        }
        std::cout << std::fixed << std::scientific << std::setprecision(15)
            << result1[i].rel_err << "\t" << result1[i].numit << "\t"
            << std::setprecision(5) << result1[i].elapsed_time << " sec\n";
    }
    std::cout << std::string(71, '-') << std::endl;

    stopwatch.stop();
    std::cout << "\nThis program ran for " << std::defaultfloat
        << stopwatch.get_elapsed_time() << " seconds." << std::endl;
    return 0;
}