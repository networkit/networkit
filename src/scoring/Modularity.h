/*
 * Modularity.h
 *
 *  Created on: 15.10.2012
 *      Author: cls
 */

#ifndef MODULARITY_H_
#define MODULARITY_H_

#include "EdgeScoring.h"
#include "../graph/Graph.h"

class Cluster; // TODO: implement in src/clustering
class Clustering; // TODO: implement in src/clustering

namespace EnsembleClustering {

class Modularity: public EnsembleClustering::EdgeScoring {

public:

	Modularity();

	virtual ~Modularity();

	virtual double scoreEdge(id u, id v);

	double mod(Cluster* C, Clustering* clustering);

	double cutweight(Cluster* c, Cluster* d);

	double weight(Cluster* c);
};

} /* namespace EnsembleClustering */
#endif /* MODULARITY_H_ */
