/*
 * EnsembleClusterer.cpp
 *
 *  Created on: 17.12.2012
 *      Author: cls
 */

#include "EnsembleClusterer.h"

namespace EnsembleClustering {

EnsembleClusterer::EnsembleClusterer() {
	// TODO Auto-generated constructor stub

}

EnsembleClusterer::~EnsembleClusterer() {
	// TODO Auto-generated destructor stub
}

Clustering& EnsembleClusterer::run(Graph& G) {

	// base clusterers calculate base clusterings
	for (Clusterer& clusterer : baseClusterers) {
		Clustering zeta = clusterer.run(G);
		baseClusterings.push_back(zeta);
	}

	// calculate overlap of base clusterings
	RegionGrowingOverlapper over;
	Clustering core = over.run(G, baseClusterings);

	// contract graph according to overlap clustering

	// submit contracted graph to base clusterers again

	// repeat until ?

	// TODO: implement

}

} /* namespace EnsembleClustering */
