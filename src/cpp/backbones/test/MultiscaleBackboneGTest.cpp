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

	g.setWeight(0, 1, 10.0);
	g.setWeight(0, 2, 1.0);
	g.setWeight(0, 3, 5.0);
	g.setWeight(0, 4, 1.0);
	g.setWeight(4, 5, 1.0);
	g.setWeight(3, 5, 2.0);



}

}
/* namespace NetworKit */

#endif /*NOGTEST */
