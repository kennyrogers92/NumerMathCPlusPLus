/* Ways of computing the machine epsilon for double-precision floating-point
   number system. */
#include <iostream>
#include <iomanip>
#include <limits>
#include <cmath>

int main() {
	double MACH_EPS[4];
	// Value using definition of machine epsilon.
	MACH_EPS[0] = std::pow(2, -52);
	// Value obtained from numeric_limits in <climits>.
	MACH_EPS[1] = std::numeric_limits<double>::epsilon();
	// Value using a for-loop.
	MACH_EPS[2] = 1.0;
	while (1.0 + 0.5*MACH_EPS[2] > 1.0)
		MACH_EPS[2] *= 0.5;
	// Value using a simple arithmetic computation.
	MACH_EPS[3] = 7.0/3.0 - 4.0/3.0 - 1.0;

	// Print elements of array MACH_EPS.
	for (int k = 0; k < 4; k++) {
		std::cout << "MACH_EPS[" << k << "] = " << std::setprecision(16)
			<< MACH_EPS[k] << std::endl;
	}

	// Check if all elements in array MACH_EPS are equal.
	bool allEqual = true;
	for (int k = 1; k < 4; k++) {
		if (MACH_EPS[0] != MACH_EPS[k]) {
			allEqual = false;
			break;
		}
	}
	if (allEqual)
		std::cout << "All elements of MACH_EPS[] are equal!" << std::endl;

	return 0;
}
