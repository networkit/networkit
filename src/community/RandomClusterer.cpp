/*
 * RandomClusterer.cpp
 *
 *  Created on: 02.02.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#include "RandomClusterer.h"
#include "../clustering/ClusteringGenerator.h"

namespace NetworKit {

RandomClusterer::RandomClusterer() {
	// TODO Auto-generated constructor stub

}

RandomClusterer::~RandomClusterer() {
	// TODO Auto-generated destructor stub
}

Clustering RandomClusterer::run(Graph& G) {
	ClusteringGenerator gen;
	return gen.makeRandomClustering(G, 42); // TODO: does this need a command line parameter?
}

} /* namespace NetworKit */
