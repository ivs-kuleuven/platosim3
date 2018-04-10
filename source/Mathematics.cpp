#include "Mathematics.h"

namespace Mathematics
{

/**
 * Computes the exponential integral Ei(x) for x > 0;
 *
 * The implementation was adopted from "Numerical Recipes in C"
 * (2nd Edition, Sect. 6.3).
 */
double expint(double x)
{
	float fact, prev, sum, term;
	unsigned int k;

	int test = Mathematics::MAX_NUM_TERMS;

	if (x <= 0.0)
		throw UnsupportedException("Bad argument in Ei");

	// Limit the number of terms in the Taylor series of Eq. (6.3.10) to a minimum
	// if x is very small (to avoid underflow in convergence test)

	if (x <  Mathematics::FPMIN)
	{
		return Mathematics::EULER + log(x);
	}

	// Use Taylor series of Eq. (6.3.10)

	if (x <= -log(Mathematics::EPS)) {

		sum = 0.0;
		fact = 1.0;

		for (k = 1; k <= Mathematics::MAX_NUM_TERMS; k++) {
			fact *= x / k;
			term = fact / k;
			sum += term;

			if (term < Mathematics::EPS * sum)
				break;
		}

		if (k > Mathematics::MAX_NUM_TERMS)
			throw UnsupportedException("Series failed in Ei");

		return Mathematics::EULER + log(x) + sum;
	}

	// Use asymptotic series of Eq. (6.3.11)
	// (start with the 2nd term)

	else {

		sum = 0.0;
		term = 1.0;

		for (k = 1; k <= Mathematics::MAX_NUM_TERMS; k++) {
			prev = term;
			term *= k / x;

			// Since the final sum is greater than one, term itself approximates
			// the relative error

			if (term < Mathematics::EPS)
				// Still converging: add new term

				break;

			if (term < prev)
				// Diverging: subtract previous term and exit

				sum += term;

			else
			{
				sum -= prev;
				break;
			}
		}

		return exp(x) * (1.0 + sum) / x;
	}
}

}
