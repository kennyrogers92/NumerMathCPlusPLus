#include <iostream>
#include <ctime>
#include "../src/max.hpp"
#include "../src/min.hpp"

int main() {
    size_t N;
    std::cout << "Enter the size of the vector: "; std::cin >> N;
    std::vector<double> v(N);
    //reset the random generator based on the system clock
    srand(unsigned(time(0)));
    // generate a random vector of size N with elements in [0, 1]
    for (size_t k = 0; k < N; k++) {
        v[k] = rand() / static_cast<double>(RAND_MAX);
    }
    // print vector
    std::cout << "\nVector: (";
    for (size_t k = 0; k < N - 1; k++) {
        std::cout << v[k] << ", ";
    }
    std::cout << v[N - 1] << ")" << std::endl;

    // get the maximum entry and the corresponding index
    long i, j;
    double m = max(v, i);
    //print result
    std::cout << "\nmax is " << m << " at " << i << "th index" << std::endl;

    // get the minimum entry and the corresponding index
    m = min(v, j);
    //print result
    std::cout << "\nmin is " << m << " at " << j << "th index" << std::endl;
    return 0;
}