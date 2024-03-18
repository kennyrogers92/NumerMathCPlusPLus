// Sample Implementation of Polynomial Interpolation Algorithms

#include "../src/polyinterp.hpp"

using namespace polyinterp;

double f(double x) {
    return std::pow(M_E, 3.*x)*sin(2.*M_PI*x);
}

double df(double x) {
    return std::pow(M_E, 3.*x)*(3.*sin(2.*M_PI*x) + (2.*M_PI*cos(2.*M_PI*x)));
}

double g(double x) {
    return 1./(1. + x*x);
}

int main() {
    vector_t x (4);
    for (int i = 0; i < 4; i++) {
        x[i] = 1./3.*i;
    }
    std::cout << "Interpolation of exp(3x)sin(2pix) in the interval [0, 1] " <<
        "with 4 equally-spaced nodes\n"; 
    std::cout << "Via Lagrange:\t" << lagrange(f, x) << std::endl;
    std::cout << "Via Newton-Lagrange:\t" << newtonlagrange(f,x) << std::endl;
    std::cout << "Via Hermite with 2 nodes:\t" <<
        hermite(f, df, vector_t{0.,1.}) << std::endl;
    std::cout << "Via Hermite with 3 nodes:\t"<<
        hermite(f, df, vector_t{0., 0.5, 1.}) << std::endl;
    x.resize(5, 0.);
    for (int i = 0; i < 5; i++) {
        x[i] = 1./4.*i;
    }

    std::cout << "Interpolation of exp(3x)sin(2pix) in the interval [0, 1] " <<
        "with 5 equally-spaced nodes\n"; 
    std::cout << "Via Lagrange:\t" << lagrange(f, x) << std::endl;
    std::cout << "Via Newton-Lagrange:\t" << newtonlagrange(f, x) << std::endl;

    std::cout << "Interpolation of 1/(1+x^2) in the interval [-5, 5] " <<
        "with nodes:\n"; 
    std::cout << "Via Chebyschev:\t" << chebyshev(g, 15, -5., 5.) << std::endl;
    std::cout << "Via Legendre:\t" << legendre(g, 15, -5., 5.) << std::endl;
    
}