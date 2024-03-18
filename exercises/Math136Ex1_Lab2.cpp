#include <iostream>
#include <cmath>
#include <iomanip>
#include <climits>
#include "../lectures/src/binary64.hpp"
#include "../lectures/src/sgn.hpp"

/**
 * @brief Chops the given float up to p-digit precision.
 * fmod in the standard library is used to get the
 * fractional part of x after multiplying 10^p. 
 * @param x float to be chopped
 * @param p number of digit precision
 * @return double chopped float
 */
double dec_chop(const Binary64 &x, int p) {
   double t = x.data.float64;
   double fracpart;
   fracpart = std::fmod(t, std::pow(10.0, -p+1));
   t -= fracpart;
   return t;
}

/**
 * @brief Rounds the given float up to p-digit precision.
 * fmod in the standard library is used to get the
 * fractional part of x after multiplying 10^p. 
 * @param x float to be rounded
 * @param p number of digit precision
 * @return double rounded float
 */
double dec_round(const Binary64 &x, int p) {
   double t = x.data.float64;
   double fracpart;
   fracpart = std::fmod(t, std::pow(10.0, -p+1));
   if (std::abs(fracpart) < (5*std::pow(10.0, -p)) ||
      std::abs(fracpart) >= (9*std::pow(10.0, -p))) {
      t = t - fracpart;
   }
   else {
      t = t - fracpart + (std::pow(10.0, -p+1)*sgn(t)); 
   }
   return t;
}

/**
 * @brief Chops the given float in its binary representation
 * to p-digit precision. The method modify_bits is available
 * in the binary64 header file.
 * @param x float to be chopped
 * @param p number of digit precision
 * @return double chopped float
 */
double bin_chop(const Binary64 &x, int p) {
   Binary64 t = x;
   for (int k = x.NBIT_MANTISSA-p-1; k > -1; k--) {
      t.modify_bits(k, 0);
   }
   return t.data.float64;
}

/**
 * @brief Rounds the given float in its binary representation
 * to p-digit precision. The method modify_bits is available
 * in the binary64 header file.
 * @param x float to be rounded
 * @param p number of digit precision
 * @return double rounded float
 */
double bin_round(const Binary64 &x, int p) {
   Binary64 t = x;
   if (x.bits[x.NBIT_MANTISSA-p-1] == 1) {
      int j = 0;
      while (x.bits[x.NBIT_MANTISSA-p+j] == 1) {
         t.modify_bits(x.NBIT_MANTISSA-p+j, 0);
         j++;
      }
      t.modify_bits(x.NBIT_MANTISSA-p+j, 1);
   }
   for (int k = x.NBIT_MANTISSA-p-1; k > -1; k--) {
      t.modify_bits(k, 0);
   }
   return t.data.float64;
}


double relative_error(double actual, double estimate) {
   return std::abs((estimate-actual)/std::abs(actual));
}

int main () {
   double x = -10*M_PI;
   int i, p = 15, p1 = 52;

   std::cout << "\nUSING DECIMAL CONVERSION: \n" << std::endl;

   std::cout << "p\t" << "CHOPPED FLOAT\t" << "RELATIVE ERROR"
      << "RATIO: REL_ERROR/PRED_ERROR" << std::endl;
   for (i = 1; i <= p; i++){
      double estimate = dec_chop(x, i);
      double rel_error = relative_error(x, estimate);
      std::cout << i << "\t" << std::fixed << std::setprecision(14)
         << estimate << "\t" << std::scientific << std::setprecision(10)
         << rel_error  << "\t"
         << rel_error/std::pow(10.0, -i + 1) << std::endl;
   }
   std::cout << "\n" << std::endl;

   std::cout << "p\t" << "CHOPPED FLOAT\t" << "RELATIVE ERROR"
      << "RATIO: REL_ERROR/PRED_ERROR" << std::endl;
   for (i = 1; i <= p; i++){
      double estimate = dec_round(x, i);
      double rel_error = relative_error(x, estimate);
      std::cout << i << "\t" << std::fixed << std::setprecision(14)
         << estimate << "\t" << std::scientific << std::setprecision(10)
         << rel_error  << "\t"
         << rel_error/std::pow(10.0, -i + 1) << std::endl;
   }

   std::cout << "\nUSING BINARY CONVERSION: \n" << std::endl;

   std::cout << "p\t" << "CHOPPED FLOAT\t" << "RELATIVE ERROR"
      << "RATIO: REL_ERROR/PRED_ERROR" << std::endl;
   for (i = 1; i <= p1; i++){
      double estimate = bin_chop(x, i);
      double rel_error = relative_error(x, estimate);
      std::cout << i << "\t" << std::fixed << std::setprecision(14)
         << estimate << "\t" << std::scientific << std::setprecision(10)
         << rel_error  << "\t"
         << rel_error/std::pow(10.0, -i + 1) << std::endl;
   }
   std::cout << "\n" << std::endl;

   std::cout << "p\t" << "CHOPPED FLOAT\t\t" << "RELATIVE ERROR"
      << "RATIO: REL_ERROR/PRED_ERROR" << std::endl;
   for (i = 1; i <= p1; i++){
      double estimate = bin_round(x, i);
      double rel_error = relative_error(x, estimate);
      std::cout << i << "\t" << std::fixed << std::setprecision(14)
         << estimate << "\t" << std::scientific << std::setprecision(10)
         << rel_error  << "\t"
         << rel_error/std::pow(10.0, -i + 1) << std::endl;
   }

   return 0;
}
