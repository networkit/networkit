/*
 * STINGERGTest.cpp
 *
 *  Created on: 10.01.2013
 *      Author: cls
 */

#include "STINGERGTest.h"

namespace EnsembleClustering {

TEST_F(STINGERGTest, testNodes) {
	Graph G;
	stinger* S = G.asSTINGER();

	// what is the weight of a node unknown to stinger

	G.setWeight(1, 42.0);

	double w = stinger_vweight(S, 17);
	DEBUG("weight of node 17: " << w);
}


} /* namespace EnsembleClustering */
