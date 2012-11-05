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

class EdgeScoring {

class Node;

public:

	EdgeScoring();

	virtual ~EdgeScoring();

	virtual double scoreEdge(Node u, Node v) = 0;
};

} /* namespace EnsembleClustering */
#endif /* EDGESCORING_H_ */
