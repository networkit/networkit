/*
 * Conductance2.h
 *
 *  Created on: 06.06.2013
 *      Author: cls
 */

#ifndef CONDUCTANCE_H
#define CONDUCTANCE_H

#include <algorithm>

#include "LocalQualityMeasure.h"

namespace NetworKit {

class Conductance: public NetworKit::LocalQualityMeasure {
public:
	Conductance();
	virtual ~Conductance();

	virtual double getQuality(std::unordered_set<node>& C, Graph& G);
};

} /* namespace NetworKit */
#endif /* CONDUCTANCE_H */
