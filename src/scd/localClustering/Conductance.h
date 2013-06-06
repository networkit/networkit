/*
 * Conductance.h
 *
 *  Created on: 06.06.2013
 *      Author: Yassine Marrakchi
 */

#ifndef CONDUCTANCE_H_
#define CONDUCTANCE_H_

#include "LocalQualityMeasure.h"
#include <unordered_set>

namespace NetworKit {

class Conductance: public NetworKit::LocalQualityMeasure {

public:

	Conductance();

	virtual ~Conductance();

	virtual double getQuality(std::unordered_set<node>& zeta, Graph& G);
};

} /* namespace NetworKit */
#endif /* CONDUCTANCE_H_ */
