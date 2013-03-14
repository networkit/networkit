/*
 * NumericTools.cpp
 *
 *  Created on: 15.02.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#include "NumericTools.h"

namespace EnsembleClustering {

NumericTools::NumericTools() {
	// TODO Auto-generated constructor stub

}

NumericTools::~NumericTools() {
	// TODO Auto-generated destructor stub
}

bool NumericTools::equal(const double x, const double y) {
	// TODO: needs testing
	double eps = std::numeric_limits<double>::epsilon();
	return (x <= (y + eps)) && (x >= (y - eps));
}

} /* namespace EnsembleClustering */
