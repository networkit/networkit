/*
 * QualityMeasure.h
 *
 *  Created on: 10.12.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef QUALITYMEASURE_H_
#define QUALITYMEASURE_H_

#include <networkit/structures/Partition.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

/**
 * @ingroup community
 * Abstract base class for all clustering quality measures.
 */
class QualityMeasure {


public:
	virtual double getQuality(const Partition& zeta, const Graph& G) = 0;
};

} /* namespace NetworKit */
#endif /* QUALITYMEASURE_H_ */
