#ifndef TIMER_HPP_INCLUDE
#define TIMER_HPP_INCLUDE

#include <ctime>

/*
 * Class implementation for the Matlab functions tic and toc.
 */
class timer {
private:
	clock_t 	tic; // begin
	clock_t 	toc; // end
	double elapsed_time;	// time in seconds

public:
	// Default constructor.
	timer() {}

	// Begin timer.
	void start() {
		tic = clock();
	}

	// End timer and store elapsed time in seconds.
	void stop() {
		toc = clock();
		elapsed_time = static_cast<double>(toc - tic)/CLOCKS_PER_SEC;
	}

	// Return the elapsed time in seconds.
	double get_elapsed_time() {
		return elapsed_time;
	}

	// Default destructor.
	~timer() {}

};	// End of timer declaration

#endif
