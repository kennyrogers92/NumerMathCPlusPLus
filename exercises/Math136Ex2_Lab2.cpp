#include <string>
#include "../lectures/src/rootspoly.hpp"

std::vector<complex_d_t> pn(int n) {
     std::vector<complex_d_t> pn;
     for (int i = 0; i < n; i++) {
          pn.push_back(-1);
     }
     pn.push_back(1);
     return pn;
}

int main() {
     roots_poly::param parameter;
     parameter.maxit = 1000;
     parameter.refmax = 100;
     parameter.tol = std::pow(10, 3) * std::numeric_limits<double>::epsilon();
     parameter.reftol = 1e-3;
     parameter.ref = true;

     for (int j = 1; j < 5; j++) {
          int n = 5 * j;
          std::complex<double> z (1, 1);
          std::vector<complex_d_t> degn_result
               = roots_poly::newtonHorner(pn(n), z, parameter).z_arr;
          std::cout << "\n> Degree " << n << std::endl;
          std::cout << std::string(75, '-') << std::endl;
          std::cout << "REAL PART\t\t\tIMAG PART\t\t\t|FUNVAL|\n";
          std::cout << std::string(75, '-') << std::endl;
          for (int i = 0; i < n; i++) {
               double funval = std::abs(roots_poly::horner(pn(n),
                    degn_result.at(i)).first)/std::numeric_limits<double>::epsilon();
               std::cout << std::showpos << std::scientific
                    << std::setprecision(15) << degn_result.at(i).real()
                    << "\t\t" << std::showpos << std::scientific
                    << std::setprecision(15) << degn_result.at(i).imag()
                    << "\t\t" << std::noshowpos << std::fixed
                    << std::setprecision(1) << funval << "eps\n";
          }
          std::cout << std::string(75, '-') << std::endl;
     }
     return 0;
}
