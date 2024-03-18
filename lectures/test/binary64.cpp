// Implementation of the Binary64 struct.
#include <cmath>
#include <limits>
#include <utility>
#include "../src/binary64.hpp"

int main() {
	/* Descriptions:
	 * 	M_PI		: constant pi(ratio of circle's perimeter to diamter
	 * 	M_E			: constant e (base of natural log)
	 * 	ZERO		: zero
	 * 	N_ZERO		: negative zero
	 * 	DBL_EPS		: floating-point epsilon
	 * 	DBL_MIN		: smallest positive floating-point number
	 * 	DBL_MAX		: largest positive floating-point number
	 * 	DENORM_MIN	: smallest denormalize numbers
	 * 	INFINITY	: infinity
	 * 	N_INFINITY	: negative infinity
	 * 	NAN			: NaN (Not a Number)
	 * 	ANOTHER_NAN	: Another NaN
	 */
	std::pair<std::string, Binary64>
	numbers[] = {
		{"M_PI",			Binary64(M_PI)},
		{"M_E",				Binary64(M_E)},
		{"ZERO",			Binary64(0.0)},
		{"N_ZERO",			Binary64(-0.0)},
		{"DBL_EPS",			Binary64(std::numeric_limits<double>::epsilon())},
		{"DBL_MIN",			Binary64(std::numeric_limits<double>::min())},
		{"DBL_MAX",			Binary64(std::numeric_limits<double>::max())},
		{"DENORM_MIN",		Binary64(0.0).modify_bits(0,1)},
		{"INFINITY",		Binary64(INFINITY)},
		{"N_INFINITY",		Binary64(-INFINITY)},
		{"NAN",				Binary64(NAN)},
		{"ANOTHER_NAN",		Binary64(NAN).modify_bits({0, 51, 63}, {0, 1, 0})}
	};
	for(std::pair<std::string, Binary64> &number : numbers) {
		std::cout << number.first << std::endl;
		std::cout << number.second << std::endl;
	}
	return 0;
}
