/*
 * NumericTools.h
 *
 *  Created on: 15.02.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef NUMERICTOOLS_H_
#define NUMERICTOOLS_H_

#include <limits>



namespace Aux {

class NumericTools {
public:

	static constexpr double machineEpsilon = std::numeric_limits<double>::epsilon();

	static constexpr double acceptableError = 1e-12;

	NumericTools();
	virtual ~NumericTools();

	template <typename T> static bool willOverflow(const T& pX, const T& pValue, const T& pMax = std::numeric_limits<T>::max()) {
		return pMax - pValue < pX;
	}

	template <typename T> static bool willUnderflow(const T& pX, const T& pValue, const T& pMin = std::numeric_limits<T>::min()) {
		return pMin + pValue > pX;
	}


	/**
	 * Test doubles for equality within a given error.
	 */
	static bool equal(const double x, const double y, const double error = acceptableError);

	/**
	 * Test doubles for equality within a given error.
	 */
	static bool le(const double x, const double y, const double error = acceptableError);

	/**
	 * Test doubles for equality within a given error.
	 */
	static bool ge(const double x, const double y, const double error = acceptableError);




};

} /* namespace NetworKit */
#endif /* NUMERICTOOLS_H_ */
