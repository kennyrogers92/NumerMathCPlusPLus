// Sample implementation of muller method method.
#include <string>
#include "../src/rootspoly.hpp"

std::vector<complex_d_t> pn(int n) {
     std::vector<complex_d_t> pn;
     for (int i = 0; i < n; i++) {
          pn.push_back(-1);
     }
     pn.push_back(1);
     return pn;
}

int main() {
     // parameter object
     roots_poly::param parameter;
     parameter.maxit = 100;
     parameter.refmax = 100;
     parameter.tol = std::pow(10, 3) * std::numeric_limits<double>::epsilon();
     parameter.reftol = 1e-3;

     // Muller-Horner without refinement
     parameter.ref = 0;
     roots_poly::Result
          deg5_result = roots_poly::newtonHorner(pn(5), 0.5, parameter);
     deg5_result.print_table();

     // Muller-Horner with refinement
     parameter.ref = 1;
     deg5_result = roots_poly::newtonHorner(pn(5), 0.5, parameter);
     deg5_result.print_table();
     return 0;
}