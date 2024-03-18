/* Calculates the IEEE 754 standard 64-bit binary floating-point
   representation of a real number. */
#ifndef BINARY64_HPP_INCLUDE
#define BINARY64_HPP_INCLUDE

#include <iostream>
#include <iomanip>
#include <bitset>

struct Binary64 {
	static const int NBIT = 63, NBIT_EXPONENT = 11, NBIT_MANTISSA = 52;
	union { double float64; unsigned long ulong64; } data;
	std::bitset<NBIT> bits;

	// Initialization with a double precision floating-point number input.
	Binary64(double number) {
		this->data.float64 = number;
		this->bits = std::bitset<NBIT>(data.ulong64);
	}

	// Initialization with an unsigned long integer.
	Binary64(unsigned long number) {
		this->data.ulong64 = number;
		this->bits = std::bitset<NBIT>(number);
	}

	// Initialization with the bits as input.
	Binary64(std::bitset<NBIT> &bits) {
		this->bits = bits;
		this->data.ulong64 = bits.to_ulong();
	}

	// Change the stored bit at idx with val.
	Binary64 modify_bits(size_t idx, int val) {
		this->bits[idx] = val;
		this->data.ulong64 = this->bits.to_ulong();
	return *this;
	}

	// Change the stored bits at idxs wth vals.
	Binary64 modify_bits(std::initializer_list<size_t> &&idxs,
		std::initializer_list<int> &&vals) {
		for (size_t k = 0; k < idxs.size(); k++)
			this->bits[*(idxs.begin() + k)] = *(vals.begin() + k);
		this->data.ulong64 = this->bits.to_ulong();
		return *this;
	}
};

int getNumZeroBitsExponent(const Binary64 &number) {
	int numZeroBitsExponent = 0;
	for (int k = number.NBIT - 2; k > number.NBIT_MANTISSA; k--) {
		if (number.bits[static_cast<size_t>(k)] == 0)
			numZeroBitsExponent++;
	}
	return numZeroBitsExponent;
}

int getNumZeroBitsMantissa(const Binary64 &number) {
	int numZeroBitsMantissa = 0;
	for (int k = number.NBIT_MANTISSA - 1; k >= 0; k--) {
		if (number.bits[static_cast<size_t>(k)] == 0)
			numZeroBitsMantissa++;
	}
	return numZeroBitsMantissa;
}

void printSignBit(const Binary64 &number) {
	std::cout << number.bits[number.NBIT-1] << " ";
}

void printExponentBits(const Binary64 &number) {
	for (int k = number.NBIT - 2; k > number.NBIT_MANTISSA - 1; k--)
		std::cout << number.bits[static_cast<size_t>(k)];
	std::cout << " ";
}

void printMantissaBits(const Binary64 &number) {
	for (int k = number.NBIT_MANTISSA - 1; k >= 0; k--)
		std::cout << number.bits[static_cast<size_t>(k)];
}

void printStatus(const Binary64 &number, int &numZeroBitsExponent, int &numZeroBitsMantissa) {
	std::cout << "\nStatus\t: ";
	//NaNs and infinities
	if (numZeroBitsExponent == 0)
		std::cout << "Exceptional";
	//De-normalized numbers
	if (numZeroBitsExponent == number.NBIT_EXPONENT
			&& numZeroBitsMantissa < number.NBIT_MANTISSA)
		std::cout << "De-normal";
	//Positive and negative 0
	if (numZeroBitsExponent == number.NBIT_EXPONENT
			&& numZeroBitsMantissa == number.NBIT_MANTISSA) {
		if (number.bits[number.NBIT-1] == 0)
			std::cout << "Normal";
		else
			std::cout << "Exceptional";
	}
	// Normalized numbers
	if (numZeroBitsExponent > 0
			&& numZeroBitsExponent < number.NBIT_EXPONENT)
		std::cout << "Normal";
}

//Output stream of a Binary64 struct.
std::ostream& operator<<(std::ostream &output, const Binary64 &number) {
	int numZeroBitsExponent = getNumZeroBitsExponent(number);
	int numZeroBitsMantissa = getNumZeroBitsMantissa(number);
	std::cout.precision(16);
	output << "Float64 : " << number.data.float64 << std::endl;
	output << "Binary64: ";
	printSignBit(number);
	printExponentBits(number);
	printMantissaBits(number);
	printStatus(number, numZeroBitsExponent, numZeroBitsMantissa);
	output << "\nHex Form: Ox" << std::hex << number.data.ulong64
		<< std::dec << "\n";
	return output;
}

#endif
