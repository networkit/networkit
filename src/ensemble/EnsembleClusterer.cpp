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

	for (Clusterer& clusterer : baseClusterers) {
		Clustering zeta = clusterer.run(G);
		baseClusterings.push_back(zeta);
	}

	// TODO: implement

}

} /* namespace EnsembleClustering */
