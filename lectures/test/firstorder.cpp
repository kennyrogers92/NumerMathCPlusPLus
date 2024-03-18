#include "../src/rootscalar.hpp"
#include <string>
#include <vector>
#include <algorithm>

// Define f(x) = 0.25*cos^2(2x) - x^2
double f(double x) {
	return 0.25*std::pow(std::cos(2.0*x), 2) - std::pow(x, 2);
}
double g(double x) {
     return -std::sin(2.0*x)*std::cos(2.0*x) - 2*x;
}

int main() {
     // parameter object
     root_scalar::param parameter;
	parameter.tol = 1e-15;
	parameter.maxit = 100;
     // approximate a root by inexact newton-raphson methods
     root_scalar::Result
          Chord = chord(f, 0.0, 0.5, 1.0, parameter),
          Secant = secant(f, 0.5, 0.0, parameter),
          RegulaFalsi = regfalsi(f, 0.5, 0.0, parameter),
          NewtonRaphson = newtonraphson(f, g, 0.5, parameter),
	  	Inexactcenter = inexactcenter(f, 0.5, parameter),
          Steffensen = steffensen(f, 0.5, parameter);
	Inexactcenter.method_name= "Inexact Center";
     std::vector<root_scalar::Result> method = {Chord, Secant,
               RegulaFalsi, NewtonRaphson, Inexactcenter, Steffensen};
     std::sort(method.begin(), method.end());
     std::cout << std::string(94, '-') << std::endl;

     std::cout << "METHOD\t\tAPPROXIMATE ROOT\t"
          << "FUNCTION VALUE\t\tERROR\t\t\tNITERS\n";
     std::cout << std::string(94, '-') << std::endl;
     for (size_t i = 0; i < 6; i++) {
          std::cout << method.at(i).method_name;
          if (method.at(i).method_name.size() < 7) {
               std::cout << "\t\t";
          }
          else {
               std::cout << "\t";
          }
          std::cout << std::fixed << std::setprecision(17)
               << method.at(i).x << "\t" << std::scientific
               << std::setprecision(12)<< method.at(i).funval << "\t"
               << method.at(i).error << "\t" << std::fixed
               << method.at(i).numit << std::endl;
     }
     std::cout << std::string(94, '-') << std::endl;
     return 0;
}
