/*
 * MultiscaleBackboneGTest.cpp
 *
 *  Created on: 20.06.2014
 *      Author: Gerd Lindner
 */

#ifndef NOGTEST

#include "MultiscaleBackboneGTest.h"

#include "../../backbones/MultiscaleBackbone.h"

namespace NetworKit {

TEST_F(MultiscaleBackboneGTest, testSimple) {
	Graph g(6, true, false);

	g.setWeight(0, 1, 0.39);
	g.setWeight(0, 2, 0.01);
	g.setWeight(0, 3, 0.5);
	g.setWeight(0, 4, 0.1);
	g.setWeight(4, 5, 1);
	g.setWeight(3, 5, 0.1);

}

}
/* namespace NetworKit */

#endif /*NOGTEST */
