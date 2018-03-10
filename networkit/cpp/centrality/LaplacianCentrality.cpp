/*
 * LaplacianCentrality.cpp
 *
 *  Created on: 08.03.2018
 *      Author: Kolja Esders
 */

#include "LaplacianCentrality.h"

namespace NetworKit {

LaplacianCentrality::LaplacianCentrality(const Graph& G): Centrality(G) {
}

void LaplacianCentrality::run() {
	scoreData = std::vector<double>(G.upperNodeIdBound(), 0.0);
	double totalLaplacianEnergy = 0.0;

	G.parallelForNodes([&](node u) {
		count degreeU = G.degree(u);
		totalLaplacianEnergy += degreeU * degreeU;
		double energyLossOnNodeDrop = degreeU * degreeU;

		G.forNeighborsOf(u, [&](node v, edgeweight ew) {
			energyLossOnNodeDrop += ew * (ew + 2 * G.degree(v));
			totalLaplacianEnergy += ew * ew;
		});

		scoreData[u] = energyLossOnNodeDrop;
	});

	G.parallelForNodes([&](node u) {
		scoreData[u] /= totalLaplacianEnergy;
	});

	hasRun = true;
}

} /* namespace NetworKit */
