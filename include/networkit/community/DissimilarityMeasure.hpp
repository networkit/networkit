/*
 * DissimilarityMeasure.h
 *
 *  Created on: 19.01.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef DISSIMILARITYMEASURE_H_
#define DISSIMILARITYMEASURE_H_

#include <networkit/graph/Graph.hpp>
#include <networkit/structures/Partition.hpp>
#include <networkit/structures/Cover.hpp>

namespace NetworKit {


/**
 * @ingroup community
 * Base class for all clustering dissimilarity measures.
 */
class DissimilarityMeasure {

public:

	virtual double getDissimilarity(const Graph& G, const Partition& first, const Partition& second) = 0;


	virtual double getDissimilarity(const Graph &G, const Cover &first, const Cover &second);
};

} /* namespace NetworKit */
#endif /* DISSIMILARITYMEASURE_H_ */
