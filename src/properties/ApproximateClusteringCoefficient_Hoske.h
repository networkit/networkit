/*
 * ApproximateClusteringCoefficient.h
 *
 *  Created on: 08.04.2013
 *      Author: dhoske, 14.11.2013
 */

#ifndef APPROXIMATECLUSTERINGCOEFFICIENT_H
#define APPROXIMATECLUSTERINGCOEFFICIENT_H

#include "../graph/Graph.h"

namespace NetworKit {
	class ApproximateClusteringCoefficient_Hoske {
	public:
		/**
		 * Calculates the approximate clustering coefficient with
		 * a probabilistic approach.
		 *
		 *   global: global or average local coefficient?
		 */
		double calculate(bool global, const Graph& G, count k);

		/**
		 * Calculates the number of iterations necessary for
		 * error probability error and variance variance.
		 */
		count niters(double variance, double error);
	};
} /* namespace NetworKit */
#endif /* APPROXIMATECLUSTERINGCOEFFICIENT_H */
