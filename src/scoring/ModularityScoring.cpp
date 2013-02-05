/*
 * Modularity.cpp
 *
 *  Created on: 15.10.2012
 *      Author: cls
 */

#include "ModularityScoring.h"

namespace EnsembleClustering {


ModularityScoring::ModularityScoring(Graph& G) : EdgeScoring(G) {
	this->omegaE = this->G->totalEdgeWeight();
}

ModularityScoring::~ModularityScoring() {
	// TODO Auto-generated destructor stub
}

double ModularityScoring::scoreEdge(node u, node v) {

	// calculate $$\Delta mod(c, d) := \frac{1}{2 \omega(E)} \left ( 2 \omega(E) \omega(c,d) - \omega(c) \omega(d) \right ) $$

	double deltaMod = (1 / (2 * omegaE)) * (2 * omegaE * G->weight(u, v) - G->weight(u) * G->weight(v));
	return deltaMod;

}




} /* namespace EnsembleClustering */
