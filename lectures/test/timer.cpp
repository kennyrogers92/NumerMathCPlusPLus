#include <iostream>
#include "../src/timer.hpp"

int main() {
     // create a timer object
     timer stopwatch;

     // start TIMER_HPP_INCLUDE
     stopwatch.start();
     // Compute sum of first ten million positive integers.
     unsigned long sum = 0;
     for (unsigned long k = 1; k < 100000001; k++) {
          sum += k;
     }
     // end TIMER_HPP_INCLUDE
     stopwatch.stop();

     std::cout << "Elapsed time for computing the sum of first\n"
          << "ten million positive integers is " << std::scientific
          << stopwatch.get_elapsed_time() << "seconds.\n" << std::endl;
     return 0;
}
