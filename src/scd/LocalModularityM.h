/*
 * LocalModularityM.h
 *
 *  Created on: 06.06.2013
 *      Author: Yassine Marrakchi
 */

#ifndef LOCALMODULARITYM_H_
#define LOCALMODULARITYM_H_

#include "LocalQualityMeasure.h"
#include <unordered_set>

namespace NetworKit {

class LocalModularityM: public NetworKit::LocalQualityMeasure {

public:

	LocalModularityM();

	virtual ~LocalModularityM();

	virtual double getQuality(std::unordered_set<node>& zeta, Graph& G);
};

} /* namespace NetworKit */
#endif /* LOCALMODULARITYM_H_ */
