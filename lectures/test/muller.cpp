#include "../src/rootscalar.hpp"
#include <string>

// Define f(x) = 0.25*cos^2(2x) - x^2
double f(double x) {
	return 0.25*std::pow(std::cos(2.0*x), 2) - std::pow(x, 2);
}

int main() {
    // parameter object
    root_scalar::param parameter;
	parameter.tol = 1e-16;
	parameter.maxit = 100;
    std::vector<double> x (4);
    for (int i = 0; i < 4; i++) {
        x[i] = 0.1*i;
    }
    // approximate a root by Muller's method and its generalizations
    std::vector<root_scalar::Result> method;
    method.push_back(muller(f, 0., -0.25, -0.5, parameter));
    method.push_back(rootpolyinterp(f, x, parameter));
    for (int i = 0; i < 4; i++) {
        x[i] = - 0.1*i;
    }
    method.push_back(sidisecant(f, x, parameter));
    std::cout << std::string(94, '-') << std::endl;
    std::cout << "METHOD\t\tAPPROXIMATE ROOT\tFUNCTION VALUE"
        << "\t\tERROR\t\t\tNITERS\n";
    std::cout << std::string(94, '-') << std::endl;
    for (size_t i = 0; i < method.size(); i++) {
        std::cout << method.at(i).method_name;
        if (method.at(i).method_name.size() < 8) {
            std::cout << "\t";
        }
        std::cout << "\t" << std::fixed
        << std::setprecision(17) << method.at(i).x << "\t"
        << std::showpos << std::scientific << std::setprecision(12)
        << method.at(i).funval << std::noshowpos << "\t"
        << method.at(i).error << "\t" << std::fixed
        << method.at(i).numit << std::endl;
    }
    std::cout << std::string(94, '-') << std::endl;
    return 0;
}
