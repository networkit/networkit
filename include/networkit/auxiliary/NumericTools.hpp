/*
 * NumericTools.hpp
 *
 *  Created on: 15.02.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef NETWORKIT_AUXILIARY_NUMERIC_TOOLS_HPP_
#define NETWORKIT_AUXILIARY_NUMERIC_TOOLS_HPP_

#include <algorithm>
#include <limits>

namespace Aux {

/**
 * Tools to deal with limited precision in numeric computations.
 */
namespace NumericTools {

static constexpr double machineEpsilon = std::numeric_limits<double>::epsilon();

static constexpr double acceptableError = 1e-12;

template <typename T>
bool willOverflow(const T &pX, const T &pValue, const T &pMax = std::numeric_limits<T>::max()) {
    return pMax - pValue < pX;
}

template <typename T>
bool willUnderflow(const T &pX, const T &pValue, const T &pMin = std::numeric_limits<T>::min()) {
    return pMin + pValue > pX;
}

/**
 * @brief Robust strict less-than for doubles with tolerance.
 *
 * Returns true iff x is strictly smaller than y by more than the adaptive
 * @ref tolerance(x, y, absTol, relTol). In other words, tiny numerical noise
 * around equality will NOT count as x < y.
 *
 * Mathematically: x < y - tolerance(x, y, absTol, relTol)
 *
 * @param x Left-hand value.
 * @param y Right-hand value.
 * @param error
 * @return True if x is meaningfully less than y; false otherwise.
 */
inline bool lt(double x, double y, double error = acceptableError);

/**
 * @brief Robust strict greater-than for doubles with tolerance.
 *
 * Returns true iff x is strictly larger than y by more than the adaptive
 * @ref tolerance(x, y, absTol, relTol). This ignores tiny roundoff so values
 * that are effectively equal do not count as x > y.
 *
 * Mathematically: x > y + tolerance(x, y, absTol, relTol)
 *
 * @param x Left-hand value.
 * @param y Right-hand value.
 * @param absTol Absolute tolerance (default 1e-12).
 * @param relTol Relative tolerance (default 1e-12).
 * @return True if x is meaningfully greater than y; false otherwise.
 */
inline bool gt(double x, double y, double error = acceptableError);

/**
 * Test doubles for equality within a given error.
 */
bool equal(double x, double y, double error = acceptableError);

/**
 * Test doubles for equality within a given error.
 */
bool le(double x, double y, double error = acceptableError);

/**
 * Test doubles for equality within a given error.
 */
bool ge(double x, double y, double error = acceptableError);

/**
 * Test doubles for equality within the smallest possible error.
 */
bool logically_equal(double a, double b, double error_factor = 1.0);

} /* namespace NumericTools */

} // namespace Aux
#endif // NETWORKIT_AUXILIARY_NUMERIC_TOOLS_HPP_
