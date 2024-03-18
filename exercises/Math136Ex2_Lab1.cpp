#include <string>
#include "../lectures/src/rootscalar.hpp"

double f(double x) {
     return std::pow(x, 5) - std::pow(x, 4) - std::pow(x, 3) - std::pow(x, 2)
          - x - 1;
}

double df(double x) {
     return 5*std::pow(x, 4) - 4*std::pow(x, 3) - 3*std::pow(x, 2)
          - 2*x - 1;
}

int main() {
     double eps = std::numeric_limits<double>::epsilon();
     root_scalar::param parameter;
	parameter.tol = eps;
	parameter.maxit = 1000;

     std::vector<std::vector<root_scalar::Result>> method_result;
     std::vector<root_scalar::Result> newton_result;
     std::vector<root_scalar::Result> secant_result;
     std::vector<root_scalar::Result> muller_result;
     for (size_t k = 0; k < 10; k++) {
          newton_result.push_back(root_scalar::
               newtonraphson(f, df, 0.2*k, parameter));
          secant_result.push_back(root_scalar::
               secant(f, 0.2*k, 0.2*k + 0.5, parameter));
          muller_result.push_back(root_scalar::
               muller(f, 0.2*k, 0.2*k + 0.5, 0.2*k + 1.0, parameter));
     }
     method_result.push_back(newton_result);
     method_result.push_back(secant_result);
     method_result.push_back(muller_result);

     for (size_t i = 0; i < 3; i++) {
          std::cout << "\n> " << method_result[i][0].method_name << " Method\n";
          std::cout << std::string(63, '-') << std::endl;
          std::cout << "K\tAPPROXIMATE ROOT\tNUMIT\t|FUNVAL|\tERROR\n";
          std::cout << std::string(63, '-') << std::endl;
          double method_ave_root = 0;
          long method_ave_numit = 0;
          double method_ave_funval = 0;
          double method_ave_error = 0;
          long method_success_iterate = 0;
          for (size_t k = 0; k < 10; k++) {
               if (method_result[i][k].termination_flag == "Success") {
                    method_ave_root += method_result[i][k].x;
                    method_ave_numit += method_result[i][k].numit;
                    method_ave_funval += std::abs(method_result[i][k].funval)/eps;
                    method_ave_error += method_result[i][k].error / eps;
                    method_success_iterate++;
               }
               std::cout << k << "\t" << std::fixed << std::setprecision(15)
               << method_result[i][k].x << "\t" << method_result[i][k].numit
               << "\t" << std::fixed << std::setprecision(1)
               << std::abs(method_result[i][k].funval)/eps << "eps\t\t"
               << std::fixed << std::setprecision(1)
               << method_result[i][k].error/eps << "eps\n";
          }
          std::cout << std::string(63, '-') << std::endl;
          std::cout << "AVE\t"
               << std::fixed << std::setprecision(15)
               << method_ave_root/method_success_iterate << "\t"
               << method_ave_numit / method_success_iterate << "\t"
               << std::fixed << std::setprecision(1)
               << method_ave_funval/method_success_iterate << "eps\t\t"
               << std::fixed << std::setprecision(1)
               << method_ave_error / method_success_iterate << "eps\n";
          std::cout << std::string(63, '-') << std::endl;
     }
     std::cout << std::endl;
     return 0;
}
