// Sample implementations of solving nonlinear systems numerically.

#include "../src/linalg.hpp"

using namespace linalg;

// Define vector-valued function f
vector_t f(const vector_t& x) {
    size_t n = x.size();
    vector_t ans (n);
    ans[0] = std::pow(x[0], 2) + std::pow(x[1], 2) - 4.;
    ans[1] = std::pow((x[0] - 1.), 2) + std::pow((x[1] + 1.), 2) - 1.;
    return ans;
}

// Define Jacobian of f
matrix_t Df(const vector_t& x) {
    size_t n = x.size();
    matrix_t ans (n, vector_t(n));
    ans[0][0] = 2.*x[0];
    ans[0][1] = 2.*x[1];
    ans[1][0] = 2.*x[0] - 2.;
    ans[1][1] = 2.*x[1] + 2.;
    return ans;
}

int main() {    
    nonlinear::param parameters;
    vector_t x {3., 1.5};
    nonlinear::newton(f, Df, x, parameters).print();
    x = vector_t {-1.5, -3.};
    nonlinear::newton(f, Df, x, parameters).print();
    x = vector_t {3., 1.5};
    matrix_t Q = Df(x);
    nonlinear::broyden(f, Q, x, parameters).print();
    x = vector_t {-1.5, -3.};
    Q = Df(x);
    nonlinear::broyden(f, Q, x, parameters).print();

    return 0;
}