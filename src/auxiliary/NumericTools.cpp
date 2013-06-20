/*
 * NumericTools.cpp
 *
 *  Created on: 15.02.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#include "NumericTools.h"

namespace NetworKit {

NumericTools::NumericTools() {
	// TODO Auto-generated constructor stub

}

NumericTools::~NumericTools() {
	// TODO Auto-generated destructor stub
}

bool NumericTools::equal(const double x, const double y, const double error) {
	// TODO: needs testing
	double eps;
	if (error == 0.0) {
		eps = std::numeric_limits<double>::epsilon();
	} else {
		eps = error;
	}
	return (x <= (y + eps)) && (x >= (y - eps));
}

} /* namespace NetworKit */
