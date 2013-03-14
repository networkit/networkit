/*
 * RandomClusterer.h
 *
 *  Created on: 02.02.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef RANDOMCLUSTERER_H_
#define RANDOMCLUSTERER_H_

#include "Clusterer.h"

#include "../base/ClusteringGenerator.h"

namespace EnsembleClustering {

class RandomClusterer: public EnsembleClustering::Clusterer {
public:
	RandomClusterer();
	virtual ~RandomClusterer();

	virtual Clustering run(Graph& G);
};

} /* namespace EnsembleClustering */
#endif /* RANDOMCLUSTERER_H_ */
