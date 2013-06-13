/*
 * Surprise.cpp
 *
 *  Created on: 21.03.2013
 *      Author: cls
 */

#include "Surprise.h"

namespace NetworKit {

Surprise::Surprise() {
	// TODO Auto-generated constructor stub

}

Surprise::~Surprise() {
	// TODO Auto-generated destructor stub
}

double Surprise::getQuality(const Clustering& zeta, const Graph& G) {


	count n = G.numberOfNodes();
	count m = G.numberOfEdges();

	count mMax = (n * (n - 1)) / 2;	// maximum number of edges


	count k = 0; 					// number of clusters
	// TODO: get number of clusters
	// TODO: get map (or list?) of cluster sizes


	count muStar; // TODO:
	count mu; // TODO:

	double sur;

	double sum = 0.0;
	count upper = std::min(muStar, m);
	for (count j = mu; j < upper; j += 1) {
		Aux::MissingMath::binomial(muStar, j);
	}
	sur = -1 * std::log(sum);

	return sur;

}

} /* namespace NetworKit */
