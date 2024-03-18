#include "../src/rootscalar.hpp"
#include <string>

double f(double x) {
	return std::pow((x*x - 1), 7) * std::log(x);
}

double df(double x) {
     return std::pow((x*x - 1), 6) *
     (14.0 * x * std::log(x + std::numeric_limits<double>::epsilon())
     + (x * x - 1) / x);
}

int main() {
     // Parameter object
     root_scalar::param parameter;
	parameter.tol = 1e-10;
	parameter.maxit = 1e4;
     // Approximate a root by modified newton-raphson methods
     root_scalar::Result
          Newton = newtonraphson(f, df, 0.8, parameter),
          Modnewton = modnewton(f, df, 0.8, 8, parameter),
          Adaptivenewton = adaptivenewton(f, df, 0.8, parameter, 1, 1e-3, 1e-2);
     std::vector<root_scalar::Result> method
          {Newton, Modnewton, Adaptivenewton};
     std::vector<std::string> name
          {"Newton", "ModifiedNewton", "AdaptiveNewton"};
     for (size_t i = 0; i < 3; i++) {
          method[i].method_name = name[i];
     }
     std::cout << std::string(94, '-') << std::endl;
     std::cout << "METHOD\t\tAPPROXIMATE ROOT\tFUNCTION VALUE"
          << "\t\tERROR\t\t\tNITERS\n";
     std::cout << std::string(94, '-') << std::endl;
     for (size_t i = 0; i < 3; i++) {
          std::cout << method.at(i).method_name;
          if (method.at(i).method_name.size() < 7) {
               std::cout << "\t";
          }
          std::cout << "\t" << std::fixed
          << std::setprecision(17) << method.at(i).x << "\t" << std::showpos
          << std::scientific << std::setprecision(12) << method.at(i).funval
          << std::noshowpos << "\t" << method.at(i).error << "\t"
          << std::fixed << method.at(i).numit << std::endl;
     }
     std::cout << std::string(94, '-') << std::endl;

     return 0;
}
