/*
 * ClusteringCoefficient_Hoske.h
 *
 *  Created on: 08.04.2013
 *      Author: dhoske, 14.11.2013
 */

#ifndef EXACTCLUSTERINGCOEFFICIENT_HOSKE_H
#define EXACTCLUSTERINGCOEFFICIENT_HOSKE_H

#include "../graph/Graph.h"

namespace NetworKit {
	namespace ExactClusteringCoefficient {
		/**
		 * Calculates the exact clustering coefficient.
		 *   global: global or average local coefficient?
		 */
		double calculate(bool global, const Graph& G);
	};
} /* namespace NetworKit */
#endif /* EXACTCLUSTERINGCOEFFICIENT_HOSKE_H */
