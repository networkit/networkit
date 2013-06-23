/*
 * NumericTools.cpp
 *
 *  Created on: 15.02.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#include "NumericTools.h"

namespace Aux {

NumericTools::NumericTools() {
	// TODO Auto-generated constructor stub

}

NumericTools::~NumericTools() {
	// TODO Auto-generated destructor stub
}

bool NumericTools::equal(const double x, const double y, const double error) {
	return (x <= (y + error)) && (x >= (y - error));
}

bool NumericTools::le(const double x, const double y, const double error) {
	return (x <= (y + error));
}

bool NumericTools::ge(const double x, const double y, const double error) {
	return (x >= (y - error));
}

} /* namespace NetworKit */
