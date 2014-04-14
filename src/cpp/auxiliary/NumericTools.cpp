/*
 * NumericTools.cpp
 *
 *  Created on: 15.02.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#include "NumericTools.h"

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

} /* namespace NumericTools */

} /* namespace NetworKit */
