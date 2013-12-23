/*
 * ExactClusteringCoefficient_Hoske.cpp
 *
 *  Created on: 14.11.2013
 *      Author: dhoske
 */

#include "ExactClusteringCoefficient_Hoske.h"
#include "../properties/GraphProperties.h"
#include <random>
#include <cmath>

namespace NetworKit { namespace ExactClusteringCoefficient {

double calculate(bool global, const Graph& G) {
    std::vector<double> coefficients = GraphProperties::localClusteringCoefficients(G);

	double sum = 0.0;
	count normalise = 0;
	G.forNodes([&] (node u) {
		if (global) {
			sum += coefficients[u] * (G.degree(u) * (G.degree(u) - 1));
			normalise += std::max(G.degree(u) * (G.degree(u) - 1), count(0));
		} else {
			sum += coefficients[u];
			if (G.degree(u) >= 2) {
				normalise++;
			}
		}

	});

	return normalise != 0 ? sum / normalise : 0.0;
}

}} /* namespace NetworKit */
