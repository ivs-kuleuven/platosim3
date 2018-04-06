#ifndef INCLUDE_MATHEMATICS_H_
#define INCLUDE_MATHEMATICS_H_

#include <cmath>
#include "Exceptions.h"

namespace Mathematics {

	double expint(double x);

	const double EULER = exp(1.0);			// Euler's constant gamma
	constexpr int MAX_NUM_TERMS = 100;		// Maximum number of iterations allowed
	constexpr double FPMIN = 1.0e-30; 		// Close to smallest representable floating-point number.
	constexpr double EPS = 6.0e-8;			// Relative error, or absolute error near the zero of Ei at x = 0.3725.

}
#endif /* INCLUDE_MATHEMATICS_H_ */
