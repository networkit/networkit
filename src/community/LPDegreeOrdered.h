/*
 * LPDegreeOrdered.h
 *
 *  Created on: 24.09.2013
 *      Author: cls
 */

#ifndef LPDEGREEORDERED_H_
#define LPDEGREEORDERED_H_

#include "Clusterer.h"

namespace NetworKit {

typedef cluster label; // a label is the same as a cluster id

/**
 * Label propagation-based community detection algorithm which
 * processes nodes in increasing order of node degree.
 */
class LPDegreeOrdered: public NetworKit::Clusterer {
public:
	LPDegreeOrdered();
	virtual ~LPDegreeOrdered();
	virtual Clustering run(Graph& G);
};

} /* namespace NetworKit */
#endif /* LPDEGREEORDERED_H_ */
