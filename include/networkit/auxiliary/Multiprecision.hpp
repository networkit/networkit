/*
 * Multiprecision.hpp
 *
 *  Created on: 21.03.2013
 *      Author: cls
 */

#ifndef NETWORKIT_AUXILIARY_MUTLIPRECISION_HPP
#define NETWORKIT_AUXILIARY_MUTLIPRECISION_HPP

#include <ttmath/ttmath.hpp>

namespace NetworKit {
    using bigfloat = ttmath::Big<TTMATH_BITS(64),TTMATH_BITS(64)>;	///< big floating point number
}

#endif // NETWORKIT_AUXILIARY_MUTLIPRECISION_HPP