// Sample implementations of solving linear systems numerically.

#include "../src/linalg.hpp"

using namespace linalg;

// Construct Lower Triangular Matrix L and Upper Triangular Matrix U
matrix_t benchmarkL(int n) {
    int k = 1;
    matrix_t A (n, vector_t (n, 0.));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= i; j++) {
            A[i][j] = std::pow((k+1.)/n, 2);
            k++;
        }
    }
    return A; 
}

matrix_t benchmarkU(int n) {
    int k = 1;
    matrix_t A (n, vector_t (n, 0.));
    for (int i = n-1; i > -1; i--) {
        for (int j = n-1; j >= i; j--) {
            A[i][j] = std::pow((k+1.)/n, 2);
            k++;
        }
    }
    return A;
}

// Print Results
int main() {
    int n = 1000;
    matrix_t L = benchmarkL(n);
    matrix_t U = benchmarkU(n);
    vector_t x (n, 1.);

    // Compute x for each method
    std::vector<direct::Result> res;
    vector_t b = L*x;
    res.push_back(direct::forwardsubrow(L, b));
    res.push_back(direct::forwardsubcol(L, b));
    b = U*x;
    res.push_back(direct::backwardsubrow(U, b));
    res.push_back(direct::backwardsubcol(U, b));
    
    std::cout << std::string(79, '-') << std::endl;
    std::cout << "METHOD\t\tMAX ABS ERROR\t\t" <<
        "RESIDUAL MAX NORM\tELAPSED TIME" << std::endl;
    std::cout << std::string(79, '-') << std::endl;
    for (int i = 0 ; i < 4; i++) {
        double error = Max_norm(x - res[i].x);
        std::cout << res[i].method_name << "\t" << std::fixed
            << std::scientific << std::setprecision(15)
            << error << "\t" << res[i].residual_max << "\t"
            << std::setprecision(5) << res[i].elapsed_time << " sec\n";
    }
    std::cout << std::string(79, '-') << std::endl;
    return 0;
}