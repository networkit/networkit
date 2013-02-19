/*
 * EdgeScoring.h
 *
 *  Created on: 15.10.2012
 *      Author: cls
 */

#ifndef EDGESCORING_H_
#define EDGESCORING_H_

#include "../graph/Graph.h"

namespace EnsembleClustering {


template<typename T>
class EdgeScoring {

protected:

	Graph* G;	//!< pointer to the graph

public:

	EdgeScoring(Graph& G);

	virtual ~EdgeScoring();

	virtual void scoreEdges() = 0;

	virtual T edgeScore(node u, node v) const = 0;
};


template<typename T>
EdgeScoring<T>::EdgeScoring(Graph& G) {
	this->G = &G;
}

template<typename T>
EdgeScoring<T>::~EdgeScoring() {
	// TODO Auto-generated destructor stub
}


} /* namespace EnsembleClustering */
#endif /* EDGESCORING_H_ */
