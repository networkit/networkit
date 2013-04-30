/*
 * MissingMath.h
 *
 *  Created on: 21.03.2013
 *      Author: cls
 */

#ifndef MISSINGMATH_H_
#define MISSINGMATH_H_

#include <cinttypes>

namespace Aux {

class MissingMath {

public:

	static int64_t binomial(int64_t n, int64_t k) {
		if (k == 0) return 1;
		if (2 * k > n) return binomial(n, n - k);

		int64_t b = n - k + 1;
		for (int64_t i = 2; i <= k; ++i) {
			b = b * (n - k + i);
			b = b / i;
		}
		return b;
	}

};

} /* namespace Aux */
#endif /* MISSINGMATH_H_ */
