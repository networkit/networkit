/*
 * NumericTools.cpp
 *
 *  Created on: 15.02.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#include "NumericTools.h"
#include <cmath>
#include <algorithm>

namespace Aux {

namespace NumericTools{

bool equal(const double x, const double y, const double error) {
	return (x <= (y + error)) && (x >= (y - error));
}

bool le(const double x, const double y, const double error) {
	return (x <= (y + error));
}

bool ge(const double x, const double y, const double error) {
	return (x >= (y - error));
}

bool logically_equal(double a, double b, double error_factor) {
	return a==b || std::abs(a-b)<std::abs(std::min(a,b))*std::numeric_limits<double>::epsilon()*error_factor;
}

} /* namespace NumericTools */

} /* namespace NetworKit */
