// Sample implementations of LU factorization to solve linear systems.

#include "../src/linalg.hpp"

using namespace linalg;

// Construct SPD matrix from the discretization of the Dirichlet Laplacian
// on the unit square
matrix_t discretize(int n);


// Print Results
int main() {
    int n = 900;
    vector_t x (n, 1.);
    matrix_t A = discretize(n);
    vector_t b = A*x;

    // Compute x for each method
    std::vector<direct::Result> res;
    res.push_back(direct::solveLU(A, b, "kji"));
    res.push_back(direct::solveLU(A, b, "jki"));
    res.push_back(direct::solveLU(A, b, "ijk"));
    
    std::cout << std::string(79, '-') << std::endl;
    std::cout << "METHOD\t\tMAX ABS ERROR\t\t" <<
        "RESIDUAL MAX NORM\tELAPSED TIME" << std::endl;
    std::cout << std::string(79, '-') << std::endl;
    for (int i = 0 ; i < 3; i++) {
        double error = Max_norm(x - res[i].x);
        std::cout << res[i].method_name << "\t" << std::fixed
            << std::scientific << std::setprecision(15)
            << error << "\t" << res[i].residual_max << "\t"
            << std::setprecision(5) << res[i].elapsed_time << " sec\n";
    }
    std::cout << std::string(79, '-') << std::endl;
    return 0;
}