/*
 * EdgeScoring.h
 *
 *  Created on: 15.10.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef EDGESCORING_H_
#define EDGESCORING_H_

#include "../graph/Graph.h"

namespace NetworKit {


/**
 * @ingroup scoring
 * Abstract base class for algorithms associating a score with an edge.
 */
template<typename T>
class EdgeScoring {

protected:

	Graph* G;	//!< pointer to the graph

public:

	EdgeScoring(Graph& G);

	virtual ~EdgeScoring();

	virtual void scoreEdges(int attrId) = 0;

	virtual T edgeScore(node u, node v) const = 0;
};


template<typename T>
EdgeScoring<T>::EdgeScoring(Graph& G) {
	this->G = &G;
}

template<typename T>
EdgeScoring<T>::~EdgeScoring() {

}


} /* namespace NetworKit */
#endif /* EDGESCORING_H_ */
